from __future__ import print_function, division

import matplotlib

matplotlib.use('Agg')

import os
import cProfile
import pstats
import multiprocess
import sys
from datetime import datetime
import traceback
import logging

import arguments
from parser import ReadFile
import prefiltering
import templates
import report
import annotation
import postfilter
from annotation import Annotator
import all_call


def run_annot_iter(reader, annotators, filters, threads, annotated_read_prev, filtered_reads_prev, filter_on=True, report_every=100000):
    """
    Runs the annotator iteratively and in parallel on all reads in reader.
    :param reader: iterator - iterative reader
    :param annotators: list(Annotator) - annotator classes for annotating reads
    :param filters: list(filters) - filters for pre-filtering reads
    :param threads: int - number of threads to use
    :param annotated_read_prev: list(int) - previously annotated reads
    :param filtered_reads_prev: int - previously filtered reads
    :param filter_on: bool - whether the filtering in subprocesses is on
    :param report_every: int - report progress every this number of reads
    :return: list(2D), int, int - list of annotations for filtered reads, number of annotated reads, number of filtered reads
    """

    # construct shared value and lock
    annotated_reads = multiprocess.Array('i', annotated_read_prev)
    filtered_reads = multiprocess.Value('i', filtered_reads_prev)
    lock = multiprocess.Lock()

    def init_pool(l, ar, fr, an, flt, f_on, report_num):
        """
        Initializes a pool with variables
        :param l: Lock - shared lock object for all processes
        :param ar: Array - shared Array for all processes
        :param fr: Array - shared Array for all processes
        :param an: Annotators - shared Annotators object for all processes
        :param flt: filter functions - filter functions for pre-filtering reads
        :param f_on: bool - filter is on?
        :param report_num: int - how often to report progress
        :return: None
        """
        global lock
        global annotated_reads
        global filtered_reads
        global annotators
        global filters
        global filter_on
        global report_number
        lock = l
        annotated_reads = ar
        filtered_reads = fr
        annotators = an
        filters = flt
        filter_on = f_on
        report_number = report_num

        if PROFILE:
            global prof
            prof = cProfile.Profile()

            def fin():
                try:
                    with open('%s/profile-%s.out' % (config['general']['output_dir'], multiprocess.current_process().pid), "w") as f:
                        pstats.Stats(prof, stream=f).strip_dirs().sort_stats("time").print_stats()
                except TypeError:  # empty prof
                    pass

            multiprocess.util.Finalize(None, fin, exitpriority=1)

    def annotate_read(read):
        """
        Annotate a single read with global Annotator object.
        :param read: triple - read to be annotated
        :return: tuple - annotated read (output of Annotator), read
        """
        # start profiler
        if PROFILE:
            global prof
            prof.enable()

        global lock
        global annotated_reads
        global filtered_reads
        global annotators
        global filters
        global filter_on
        global report_number

        try:
            # filter and annotate a read
            annotation = [None for _ in range(len(annotators))]
            for i, (annotator, filter) in enumerate(zip(annotators, filters)):
                annot = None
                if filter.filter_read(read):
                    annot = annotator.annotate(read)
                    annotation[i] = annot
                    annotated_reads[i] += 1
                    if filter_on and not annotation[i].is_annotated_right():
                        annotation[i] = None

                    # was this read correctly annotated?
                    if config['general']['verbosity'] > 1:
                        lock.acquire()
                        report.log_str('Read %15s annotated with Annotator %s Filter %s - %s' % (
                            str(read), annotator, filter, 'OK' if annotation[i] is not None else ('FLT_FAIL' if annot is None else 'ANN_FAIL')))
                        lock.release()

            # write down progress
            filtered_reads.value += 1

            # print progress
            if report_number > 0 and filtered_reads.value % report_number == 0:
                lock.acquire()
                s = "\r    processed: %8d (passed filtering: %s) %s" % (filtered_reads.value, " ".join(map(str, annotated_reads)), '{duration}'.format(duration=datetime.now() - start_time))
                report.log_str(s, stdout_too=False, priority=logging.DEBUG, flush=True)
                sys.stdout.write(s)
                sys.stdout.flush()
                lock.release()

        except Exception as e:
            print("Exception in child process:")
            traceback.print_exc()
            raise e

        # stop profiler.
        if PROFILE:
            prof.disable()

        # return result
        return annotation

    # crate pool and initialize it with init_pool
    pool = multiprocess.Pool(threads, initializer=init_pool, initargs=(lock, annotated_reads, filtered_reads, annotators, filters, filter_on, report_every))

    # run all annotation and return
    res = [[] for _ in range(len(annotators))]
    for partial_res in pool.map(annotate_read, reader, chunksize=100):
        for i, pr in enumerate(partial_res):
            if pr is not None:
                res[i].append(pr)

    pool.close()
    pool.join()
    print("")
    return res, filtered_reads.value, annotated_reads


"""
Start of real code: whole dante annotation algorithm.
"""

# print time of the start:
start_time = datetime.now()
print(templates.START_TIME.format(start=start_time))

# load arguments
config = arguments.load_arguments()

# start profiling if needed
prog_prof = None
PROFILE = config['general']['profiler']
if PROFILE:
    prog_prof = cProfile.Profile()
    prog_prof.enable()

# initialize logging module:
report.configure_logger("%s/dante.log" % config['general']['output_dir'])
report.log_str(templates.START_TIME.format(start=start_time), stdout_too=False)

# print arguments (too clunky)
# report.log_str(arguments.save_arguments(config))

# deduplicated reads
dedup_ap = [[] for _ in range(len(config['motifs']))]

# go through all motifs and construct annotators and filters:
annotators = []
filters = []
for i, motif in enumerate(config['motifs']):
    sequence = ','.join([module_dict['seq'] for module_dict in motif['modules']])
    motif_tuple = report.seq_into_tuple(sequence)

    # initialize annotator
    report.log_str("Motif %12s: Constructing annotator for sequence %s" % (motif['full_name'], sequence))
    annotators.append(Annotator(motif_tuple, config['annotation']['delete_prob'], config['annotation']['insert_prob'],
                                config['annotation']['max_delete_skips'], config['annotation']['motif_frequency'],
                                config['annotation']['snp_chance']))

    # initialize filter
    read_filter = prefiltering.create_filter(motif['prefilter'], motif_tuple)
    report.log_str("Motif %12s: Filter constructed: %s" % (motif['full_name'], str(read_filter)))
    filters.append(read_filter)

# main annotation
all_reads = 0
annotations = [[] for _ in range(len(annotators))]
annotated_reads = [0] * len(annotations)
readers = []
for input_file in config['inputs']:

    # initialize reader
    read_filename = input_file['path']
    file_type = None if input_file['filetype'] == 'infer' else input_file['filetype']
    max_reads = None if input_file['max_reads'] == 'all' else input_file['max_reads']

    if read_filename == 'sys.stdin':
        read_file = ReadFile(None, config['general']['stranded'], max_reads, file_type, verbosity=config['general']['verbosity'])
    else:
        read_file = ReadFile(read_filename, config['general']['stranded'], max_reads, file_type, verbosity=config['general']['verbosity'])

    report.log_str("Parsing file %s (%s type) with parser %s %s " % (read_filename, read_file.file_type, read_file.__class__.__name__, read_file.reader.__name__))

    # run annotators
    report.log_str("Running %d annotator(s) on %2d process(es) for file %s" % (len(annotators), config['general']['cpu'], read_filename))
    readers.append(read_file)
    next_annotations, cur_reads, annotated_reads = run_annot_iter(read_file, annotators, filters, config['general']['cpu'], annotated_reads, all_reads,
                                                                  not config['general']['output_all'], config['general']['report_every'])
    new_reads = cur_reads - all_reads
    all_reads += new_reads
    for i, pr in enumerate(next_annotations):
        annotations[i].extend(pr)

    # write stats for current file
    report.log_str("Reads: %10d, annotated: %s" % (new_reads, ' '.join(map(str, map(len, next_annotations)))))

# write read distribution
report.write_read_distribution('%s/read_distr.npy' % config['general']['output_dir'], readers)

# write stats
for i, (motif, annot) in enumerate(zip(config['motifs'], annotations)):

    report.log_str('Motif %12s: Kept %8d/%8d reads' % (motif['full_name'], len(annot), all_reads))

    # setup motif sequences and dir
    motif_dir = '%s/%s' % (config['general']['output_dir'], motif['full_name'])
    sequence = ','.join([module_dict['seq'] for module_dict in motif['modules']])
    motif_tuple = report.seq_into_tuple(sequence)

    # convert to annotation pairs, deduplicate and convert back
    annotation_pairs = annotation.annotations_to_pairs(annot)
    dedup_ap[i], duplicates = annotation.remove_pcr_duplicates(annotation_pairs)

    # report them
    if not os.path.exists(motif_dir):
        os.makedirs(motif_dir)
    report.write_annotation_pairs('%s/annotation_pairs.txt' % motif_dir, dedup_ap[i])
    report.write_annotation_pairs('%s/annotation_pairs_duplicates.txt' % motif_dir, duplicates)

    # log it
    report.log_str('Motif %12s: %8d reads -- %8d pairs (%8d deduplicated + %8d PCR duplicates)' % (
        motif['full_name'], len(annot), len(annotation_pairs), len(dedup_ap[i]), len(duplicates)))

    # go through all motifs
    for j, pstf in enumerate(motif['postfilter']):
        # setup post filtering - no primers, insufficient quality, ...
        postfilter_class = postfilter.Postfilter(pstf, motif_tuple)

        # get indices
        index_rep, index_rep2 = postfilter_class.get_indexes()

        # write index into config back:
        config['motifs'][i]['postfilter'][j]['index_rep'] = index_rep

        # deduplicated annotations (for each pair we keep only one):
        dedup_annot = annotation.pairs_to_annotations_pick(dedup_ap[i], index_rep - 1)

        # log it
        report.log_str('Motif %12s: Extracted %8d reads from %8d pairs' % (motif['full_name'], len(dedup_annot), len(dedup_ap[i])))

        # get filtered stuff
        report.log_str("Motif %12s: Running post-filtering for %s required repetitions and %s required bases" % (motif['full_name'], pstf['repetitions'], pstf['bases']))
        qual_annot, primer_annot, filt_annot = postfilter_class.get_filtered(dedup_annot)
        report.log_str("Motif %12s: Post-filtered %5d (at least one primer %5d), remained %5d" % (motif['full_name'], len(filt_annot), len(primer_annot), len(qual_annot)))

        # write it to files
        report.log_str("Motif %12s: Generating output files into %s" % (motif['full_name'], motif_dir))
        report.write_all(qual_annot, primer_annot, filt_annot, dedup_ap[i], all_reads, motif_dir, motif['modules'], index_rep, index_rep2, j)

# -------- All_Call part of DANTE

# run all_call
for i, motif in enumerate(config['motifs']):

    sequence = ','.join([module_dict['seq'] for module_dict in motif['modules']])
    motif_tuple = report.seq_into_tuple(sequence)
    motif_dir = '%s/%s' % (config['general']['output_dir'], motif['full_name'])

    for j, pstf in enumerate(motif['postfilter']):

        # setup post filtering - no primers, insufficient quality, ...
        postfilter_class = postfilter.Postfilter(pstf, motif_tuple)

        # get indices
        index_rep, index_rep2 = postfilter_class.get_indexes()

        # write index into config back:
        config['motifs'][i]['postfilter'][j]['index_rep'] = index_rep

        # and strings
        rep_seq = motif['modules'][index_rep - 1]['seq']
        len_str = len(rep_seq.split('-')[1])

        # is it normal or 2-index?
        if index_rep2 is not None:
            rep_seq = motif['modules'][index_rep2 - 1]['seq']
            # len_str2 = len(rep_seq.split('-')[1])
            continue  # don't do a thing for now

        # deduplicated annotations (for each pair we keep only one):
        dedup_annot = annotation.pairs_to_annotations_pick(dedup_ap[i], index_rep - 1)

        # get filtered stuff
        qual_annot, primer_annot, filt_annot = postfilter_class.get_filtered(dedup_annot)
        postfilter_bases = postfilter_class.get_filtering_bases()

        # load read distribution
        read_distribution = report.load_read_distribution('%s/read_distr.npy' % config['general']['output_dir'])

        # run inference
        inference = all_call.Inference(read_distribution, config['allcall']['param_file'], str_rep=len_str,
                                       minl_primer1=postfilter_bases[index_rep - 2], minl_primer2=postfilter_bases[index_rep], minl_str=postfilter_bases[index_rep - 1])
        file_pcolor = '%s/pcolor_%d' % (motif_dir, j + 1)
        file_output = '%s/allcall_%d.txt' % (motif_dir, j + 1)
        inference.all_call(qual_annot, primer_annot, index_rep - 1, file_pcolor, file_output, motif['full_name'])

        # write the report
        confidence = report.read_all_call('%s/allcall_%d.txt' % (motif_dir, j + 1))
        if confidence is not None:
            conf, a1, a2, c1, c2 = confidence
            if isinstance(a1, (int, long)) and a1 > 0:
                report.write_alignment('%s/alignment_%d_a%d.fasta' % (motif_dir, j + 1, a1), qual_annot, index_rep - 1, allele=a1)
            if isinstance(a2, (int, long)) and a2 != a1 and a2 != 0:
                report.write_alignment('%s/alignment_%d_a%d.fasta' % (motif_dir, j + 1, a2), qual_annot, index_rep - 1, allele=a2)

# -------- generation of reports and finalizing

# generate report and output files for whole run
report.log_str('Generating final report')
report.write_report(config['general']['output_dir'], config['motifs'], config['general']['output_dir'])

# print the time of the end:
end_time = datetime.now()
report.log_str('DANTE Stopping    : {finish:%Y-%m-%d %H:%M:%S}'.format(finish=end_time))
report.log_str('Total time of run : {duration}'.format(duration=end_time - start_time))

# stop profiler:
if PROFILE:
    prog_prof.disable()
    with open('%s/profile-main.out' % config['general']['output_dir'], "w") as f:
        pstats.Stats(prog_prof, stream=f).strip_dirs().sort_stats("time").print_stats()
