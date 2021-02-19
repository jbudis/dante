import pysam
import gzip
import logging
import report
import sys
import io
import collections
import numpy as np

HIGHEST_QUALITY = 'I'


def is_seq(sequence):
    """
    Checks if this is a DNA sequence.
    :param sequence: str - string to check
    :return: bool - true if this is a DNA sequence
    """
    return len(sequence) > 3 and sequence[0] in 'ACTGN-' and sequence[1] in 'ACTGN-' and sequence[2] in 'ACTGN-'


def extract_pair(read_id, previous_left_pair=None):
    """
    Parse the read_id and find pair information and shorten it.
    :param read_id: str - read id with pair info
    :param previous_left_pair: bool/None - is this left pair? for consistency checking
    :return: tuple(str, bool/None) - read_id without pair info, extracted pair info
    """
    split = read_id.split()
    left_pair = None
    if split[-1][0] == '1':
        left_pair = True
    elif split[-1][0] == '2':
        left_pair = False

    if left_pair is not None:
        if previous_left_pair is not None:
            assert previous_left_pair == left_pair

    return read_id, left_pair


class Read:
    """
    Class for encapsulating the read object.
    """
    
    translate_table = str.maketrans('ACGT', 'TGCA')    
    
    def __init__(self, name, sequence, quality=None, map_qual=None, chromosome=None, ref_start=None, ref_end=None, left_pair=None, complement=False):
        """
        Initialize the Read object
        :param name: str - read id
        :param sequence: str - sequence
        :param quality: str - sequencing quality
        :param map_qual: int - mapping quality
        """
        self.name = name
        self.complement = complement
        self.sort_name = name.split()[0]
        self.sequence = sequence.upper()
        self.quality = quality if quality is not None else HIGHEST_QUALITY * len(self.sequence)
        self.map_qual = map_qual
        self.chromosome = chromosome
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.left_pair = left_pair

    def __len__(self):
        """
        Return the lengths of the sequence
        :return: int - length of the sequence
        """
        return len(self.sequence)

    def __str__(self):
        """
        Return string representation of this read.
        :return: str
        """
        if self.chromosome is not None:
            rs = str(self.ref_start if self.ref_start is not None else 'NA')
            re = str(self.ref_end if self.ref_end is not None else 'NA')
            return '%s %s:%s-%s' % (self.sequence, self.chromosome, rs, re)
        else:
            return self.sequence

    def is_left_current(self):
        """
        Return if the read is left-to-right mapped.
        :return: bool/None - read is left-to-right mapped/unknown mapping side
        """
        # return left_pair xor complement
        if self.left_pair is None:
            return None
        return self.left_pair != self.complement

    def reverse_complement(self):
        """
        Return a reverse complement of this read.
        :return: Read - reversed read
        """

        reversed_sequence = self.sequence[::-1].translate(self.translate_table)
        reversed_quality = self.quality[::-1]
        return Read(self.name, reversed_sequence, reversed_quality, self.map_qual, chromosome=self.chromosome, ref_start=self.ref_start, ref_end=self.ref_end,
                    left_pair=self.left_pair, complement=True)


class ReadFile:
    """
    Class for file-reading of reads (sam/bam/fastq/fasta).
    """

    SUPPORTED_FORMATS = ('sam', 'bam', 'fasta', 'fastq', 'fasta.gz', 'fastq.gz', 'fa', 'fq', 'fa.gz', 'fq.gz', 'txt')
    STRANDED_TYPES = ('yes', 'both', 'reverse')

    def __init__(self, file_path, stranded, maximal_reads=None, file_type=None, verbosity=0, chromosome=None, ref_start=None, ref_end=None, unmapped=False):
        """
        Initialize ReadFile.
        :param file_path: str/None - file path or None if reads come from stdin
        :param maximal_reads: int - read maximal of this number of reads
        :param stranded: 'yes'/'both'/'reverse' - whether to use stranded
        :param file_type: str - file type
        :param verbosity: int - how verbose are we
        :param chromosome: str - chromosome
        :param ref_start: int - reference start
        :param ref_end: int - reference end
        :param unmapped: boolean - include unmapped reads?
        """
        self.file_path = file_path
        self.maximal_reads = maximal_reads
        self.file_type = file_type if file_type else self.estimate_file_type()
        self.n_parsed_reads = 0
        assert stranded in ReadFile.STRANDED_TYPES
        self.iter_standard = stranded in ['yes', 'both']
        self.iter_reversed = stranded in ['reverse', 'both']

        self.distribution = collections.defaultdict(int)
        
        self.reader = self.load_reader(verbosity, chromosome, ref_start, ref_end, unmapped)

    def load_reader(self, verbosity, chromosome=None, ref_start=None, ref_end=None, unmapped=False):
        """
        Return reader according to file type and file path.
        :param verbosity: int - how verbose are we
        :param chromosome: str - chromosome
        :param ref_start: int - reference start
        :param ref_end: int - reference end
        :param unmapped: boolean - include unmapped reads?
        :return: reader object
        """
        readers = {
            'sam': lambda fn: ReadFile.iter_seqs_bam(fn, verbosity),
            'bam': lambda fn: ReadFile.iter_seqs_bam(fn, verbosity, chromosome, ref_start, ref_end, unmapped),
            'fasta': lambda fn: ReadFile.iter_seqs_fasta(fn, False),
            'fastq': lambda fn: ReadFile.iter_seqs_fastq(fn, False),
            'fasta.gz': lambda fn: ReadFile.iter_seqs_fasta(fn, True),
            'fastq.gz': lambda fn: ReadFile.iter_seqs_fastq(fn, True),
            'fa': lambda fn: ReadFile.iter_seqs_fasta(fn, False),
            'fq': lambda fn: ReadFile.iter_seqs_fastq(fn, False),
            'fa.gz': lambda fn: ReadFile.iter_seqs_fasta(fn, True),
            'fq.gz': lambda fn: ReadFile.iter_seqs_fastq(fn, True),
            'txt': lambda fn: ReadFile.iter_seqs_dante(fn)
        }
        return readers[self.file_type](self.file_path)

    def __iter__(self):
        return self.iter_reads()

    def iter_reads(self):
        """
        Iterate over the reads from reader.
        :return: iterator - over the reads
        """
        self.n_parsed_reads = 0
        for read in self.reader:
            if self.maximal_reads is not None and self.n_parsed_reads >= self.maximal_reads:
                raise StopIteration
            self.distribution[len(read)] += 1
            if self.iter_standard:
                self.n_parsed_reads += 1
                yield read
            if self.iter_reversed:
                self.n_parsed_reads += 1
                yield read.reverse_complement()

    def estimate_file_type(self):
        """
        Estimate file type from the suffix.
        :return: str - suffix (file type)
        """
        if self.file_path is None:
            error = 'Unknown file type of read file, use one of supported suffix or specify file type with "filetype" argument'
            report.log_str(error, priority=logging.ERROR)
            raise NotImplementedError(error)
        for suffix in ReadFile.SUPPORTED_FORMATS:
            if self.file_path.endswith('.%s' % suffix):
                return suffix
        error = 'Unknown file type of read file, use one of supported suffix or specify file type with "filetype" argument'
        report.log_str(error, priority=logging.ERROR)
        raise NotImplementedError(error)

    def get_distribution(self, up_to=None):
        """
        Returns the distribution of reads.
        :param up_to: int - the length of the return array, if None, then it goes to maximal
        :return: ndarray - distribution of read lengths
        """
        max_length = max(self.distribution.keys()) + 1
        if up_to is None:
            up_to = max_length

        distrib_array = np.zeros(max_length, dtype=int)
        for k, v in self.distribution.values():
            distrib_array[k] = v

        return distrib_array[:up_to]

    @staticmethod
    def iter_seqs_bam(file_name, verbosity=0, chromosome=None, pos_start=None, pos_end=None, unmapped=False):
        """
        Reader for bam files.
        :param file_name: str - file path to open
        :param verbosity: int - how verbose are we 0-3
        :param chromosome: str - chromosome
        :param pos_start: int - reference start
        :param pos_end: int - reference end
        :param unmapped: boolean - include unmapped reads?
        :return: iterator - reads a bam file iteratively
        """

        if file_name is not None:
            bam = pysam.AlignmentFile(file_name, "rb")
        else:
            bam = pysam.AlignmentFile(sys.stdin, "rb") 

        if chromosome is not None and pos_start is not None and pos_end is not None:
            # read reads from specific region
            try:
                region = bam.fetch(chromosome, pos_start, pos_end)
            except ValueError:
                print(f"Detected BAM file {file_name} without index, please sort and index with 'samtools sort' and 'samtools index'.")
                sys.exit()
        else:
            region = bam
            
        for read in region:
            if unmapped and not read.is_unmapped:
                continue
            left_pair = None
            if read.is_paired:
                left_pair = read.is_read1
            try:
                ref_start = read.reference_start
                ref_end = read.reference_end
                ref_id = str(read.reference_name)
                mapq = read.mapping_quality
                ref_id, left_pair_from_name = extract_pair(ref_id, None)
                if left_pair_from_name is not None and left_pair is not None and left_pair_from_name != left_pair and verbosity > 0:
                    warn = "WARNING: read inconsistency (left pair-end should end with '1', right with '2'): %s" % (str(read))
                    report.log_str(warn, priority=logging.WARNING)
                if left_pair is None:
                    left_pair = left_pair_from_name
            except ValueError:
                # unaligned read?
                ref_start = ref_end = ref_id = mapq = None
                
            yield Read(read.qname, str(read.seq), read.qual, mapq, ref_id, ref_start, ref_end, left_pair=left_pair)
        bam.close()

    @staticmethod
    def iter_seqs_fastq(file_name, is_gzipped):
        """
        Reader for fastq files.
        :param file_name: str - file path to open
        :param is_gzipped: bool - true if the input file is gzipped
        :return: iterator - reads a fastq file iteratively
        """
        left_pair = None
        if file_name is not None:
            if '_R1' in file_name:
                left_pair = True
            if '_R2' in file_name:
                left_pair = False

        # log the pair info
        report.log_str('Reading: %s (left_pair=%s)' % (str(file_name), str(left_pair)))

        if file_name is not None:
            reads = io.BufferedReader(gzip.open(file_name)) if is_gzipped else open(file_name)
        else:
            reads = sys.stdin
        rid = "NO_ID"
        sequence = ""
        left_pair_cur = None
        for i, line in enumerate(reads):
            if not isinstance(line, str):
                line = line.decode("utf-8")
            if i % 4 == 0:
                rid = line.strip()[1:]
                rid, left_pair_cur = extract_pair(rid, left_pair)
            elif i % 4 == 1:
                sequence = line.strip()
            elif i % 4 == 3:
                quality = line.strip()
                yield Read(rid, sequence, quality, left_pair=left_pair_cur)
        if reads != sys.stdin:
            reads.close()

    @staticmethod
    def iter_seqs_fasta(file_name, is_gzipped):
        """
        Reader for fasta files.
        :param file_name: str - file path to open
        :param is_gzipped: bool - true if the input file is gzipped
        :return: iterator - reads a fasta file iteratively
        """
        left_pair = None
        if type(file_name) is str:
            if 'R1' in file_name:
                left_pair = True
            if 'R2' in file_name:
                left_pair = False

        if file_name is not None:
            reads = io.BufferedReader(gzip.open(file_name)) if is_gzipped else open(file_name)
        else:
            reads = sys.stdin
        seq, rid = '', None
        left_pair_cur = None
        for line in reads:
            if not isinstance(line, str):
                line = line.decode("utf-8")
            if line[0] == '>':
                if seq:
                    yield Read(rid, seq, left_pair=left_pair_cur)
                seq, rid = '', line[1:-1]
                rid, left_pair_cur = extract_pair(rid, left_pair)
            else:
                seq += line[:-1]
        if rid is not None:
            yield Read(rid, seq, left_pair=left_pair_cur)
        if reads != sys.stdin:
            reads.close()

    @staticmethod
    def iter_seqs_dante(file_name):
        """
        Reader for dante files.
        :param file_name: str - file path to open
        :return: iterator - reads a dante annotation file iteratively
        """
        if file_name is not None:
            reads = open(file_name)
        else:
            reads = sys.stdin
        rid = "NO_ID"
        left_pair = None
        complement = False
        for i, line in enumerate(reads):
            if not isinstance(line, str):
                line = line.decode("utf-8")
            if line[0] == '>':
                rid = line[1:].strip()
                rid, left_pair = extract_pair(rid)
                complement = False
            elif rid != "NO_ID":
                if line.startswith('line=R'):
                    left_pair = None
                    if line[6] == '1':
                        left_pair = True
                    if line[6] == '2':
                        left_pair = False
                if line.startswith('complement='):
                    complement = line.strip().split('=')[1] == 'True'
                if is_seq(line):
                    sequence = line.strip().replace('-', '')
                    yield Read(rid, sequence, left_pair=left_pair, complement=complement)
                    rid = "NO_ID"
                    left_pair = None
                    complement = False
        if reads != sys.stdin:
            reads.close()
