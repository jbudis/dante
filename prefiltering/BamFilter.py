from prefiltering.DummyFilter import DummyFilter


class BamFilter(DummyFilter):
    """
    Filter class for BAM filtering according to the chromosome and reference starts/ends
    """

    def __init__(self, chromosome, ref_start, ref_end, min_mapq=None, overlap=1):
        """
        Initialize BamFilter object
        :param chromosome: int - chromosome of interest in bam files
        :param ref_start: int - reference start
        :param ref_end:  int - reference end
        :param overlap: int - needed overlap
        """
        super(self.__class__, self).__init__()
        self.chromosome = chromosome
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.overlap = overlap
        self.min_mapq = min_mapq

    @staticmethod
    def get_overlap(a_start, a_end, b_start, b_end):
        return max(0, min(a_end, b_end) - max(a_start, b_start))

    def filter_read(self, read):
        """
        Test whether the read has passsed the filter.
        :param read: Read - the specified read.
        :return: bool - True if read has passed the filter
        """
        try:
            # if self.min_mapq is not None and read.map_qual < self.min_mapq:
            #    return False
            return read.chromosome == self.chromosome and self.get_overlap(self.ref_start, self.ref_end, read.ref_start, read.ref_end) >= self.overlap
        except (TypeError, AttributeError):
            # print('Read w/o chromosome:', read)
            return True

    def encapsulate_reader(self, reader):
        """
        Encapsulates reader object with dummy filter.
        :param reader: iterator over sequences
        :return: reader with dummy filter
        """
        for read in reader:
            if self.filter_read(read):
                yield read

    def __str__(self):
        return 'BamFilter: %s:%d-%d' % (str(self.chromosome), self.ref_start, self.ref_end)
