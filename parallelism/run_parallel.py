import traceback
import logging
import multiprocess
import sys
import pstats
import typing
import cProfile
import datetime

import prefiltering
import report
import annotation

MOTIF_PRINT_LEN = 10


def shorten_str(string: str, max_length: int, ellipsis: str = '...') -> str:
    """
    Shorten string to max_length and include ellipsis.
    :param string: str - string to shorten
    :param max_length: int - maximum length to shorten to
    :param ellipsis: str - ellipsis to add to end of string
    :return: str - shortened string
    """
    if len(string) > max_length:
        return string[:max_length - len(ellipsis)] + ellipsis
    else:
        return string


def construct_annotator(motif: dict, config: dict) -> typing.Tuple[annotation.Annotator, prefiltering.DummyFilter]:
    """
    Construct Annotator and prefilter based on motif.
    :param motif: dict - motif specification
    :param config: dict - configuration dict
    :return: Annotator, *Filter - annotator and prefilter constructed
    """
    sequence = ','.join([module_dict['seq'] for module_dict in motif['modules']])
    motif_tuple = report.seq_into_tuple(sequence)

    # initialize annotator
    if not config['general']['quiet_mode']:
        report.log_str(f'Motif {shorten_str(motif["full_name"], MOTIF_PRINT_LEN):<{MOTIF_PRINT_LEN}s}: Constructing annotator for '
                       f'sequence {sequence}')
    annotator = annotation.Annotator(motif_tuple, config['annotation']['delete_prob'], config['annotation']['insert_prob'],
                          config['annotation']['max_delete_skips'], config['annotation']['motif_frequency'],
                          config['annotation']['snp_chance'])

    # initialize filter
    read_filter = prefiltering.create_filter(motif['prefilter'], motif_tuple)

    return annotator, read_filter


def create_annotators(config: dict) -> typing.Tuple[typing.List[annotation.Annotator], typing.List[prefiltering.DummyFilter]]:
    """
    Create all annotators. Parallel if more than 100.
    :param config: dict - configuration dict
    :return: list(Annotator), list(*Filter) - annotators and corresponding filters
    """
    # go through all motifs and construct annotators and filters:
    annotators = []
    filters = []

    # create annotator generator
    if len(config['motifs']) > 100 and config['general']['cpu'] > 1:
        report.log_str('Annotators processing in parallel on {cpu} cores.'.format(cpu=config['general']['cpu']))
        pool = multiprocess.Pool(config['general']['cpu'])
        results = pool.istarmap(construct_annotator, [(motif, config) for motif in config['motifs']])
    else:
        report.log_str('Annotators processing sequentially.')
        results = (construct_annotator(motif, config) for motif in config['motifs'])
    # fill annotators
    for i, (annotator, read_filter) in enumerate(results):
        if config['general']['quiet_mode'] and i % 100 == 0:
            report.log_str(
                'Constructing annotators: {i:6d}/{cnt:6d}   {date:%Y-%m-%d %H:%M:%S}'.format(i=i, cnt=len(config['motifs']),
                                                                                             date=datetime.datetime.now()))
        annotators.append(annotator)
        filters.append(read_filter)

    return annotators, filters


def run_annot_iter(reader, annotators, filters, threads, annotated_read_prev, filtered_reads_prev, config, start_time):
    """
    Runs the annotator iteratively and in parallel on all reads in reader.
    :param reader: iterator - iterative reader
    :param annotators: list(Annotator) - annotator classes for annotating reads
    :param filters: list(filters) - filters for pre-filtering reads
    :param threads: int - number of threads to use
    :param annotated_read_prev: list(int) - previously annotated reads
    :param filtered_reads_prev: int - previously filtered reads
    :param config: dict - configuration
    :param start_time: DateTime - start of Dante run
    :return: list(2D), int, int - list of annotations for filtered reads, number of annotated reads, number of filtered reads
    """
    # set up aliases
    filter_on = not config['general']['output_all'],
    report_every = 0 if config['general']['quiet_mode'] else config['general']['report_every']

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
        global report_every
        lock = l
        annotated_reads = ar
        filtered_reads = fr
        annotators = an
        filters = flt
        filter_on = f_on
        report_every = report_num

        if config['general']['profiler']:
            global prof
            prof = cProfile.Profile()

            def fin():
                try:
                    with open('%s/profile-%s.out' % (config['general']['output_dir'], multiprocess.current_process().pid), 'w') as f:
                        pstats.Stats(prof, stream=f).strip_dirs().sort_stats('time').print_stats()
                except TypeError:  # empty prof
                    pass

            multiprocess.util.Finalize(None, fin, exitpriority=1)

    def annotate_read_sequentially(read, annotated_reads, filtered_reads, annotators, filters, filter_on, report_every, lock=None):
        """
        Annotate a single read with global Annotator object.
        :param read: triple - read to be annotated
        :param l: Lock - shared lock object for all processes
        :param annotated_reads: Array - shared Array for all processes
        :param filtered_reads: Array - shared Array for all processes
        :param annotators: Annotators - shared Annotators object for all processes
        :param filters: filter functions - filter functions for pre-filtering reads
        :param filter_on: bool - filter is on?
        :param report_every: int - how often to report progress
        :param lock: lock object - if specified then use lock for printing
        :return: tuple - annotated read (output of Annotator), read
        """
        # filter and annotate a read
        annotations = [None for _ in range(len(annotators))]
        for i, (annotator, filter) in enumerate(zip(annotators, filters)):
            if filter.filter_read(read):
                annot = annotator.annotate(read)
                annotations[i] = annot
                annotated_reads[i] += 1
                if filter_on and not annotations[i].is_annotated_right():
                    annotations[i] = None

                # was this read correctly annotated?
                if config['general']['verbosity'] > 1:
                    if lock is not None:
                        lock.acquire()
                    report.log_str('Read %15s annotated with Annotator %s Filter %s - %s' % (
                        str(read), annotator, filter, 'OK' if annotations[i] is not None else ('FLT_FAIL' if annot is None else 'ANN_FAIL')))
                    if lock is not None:
                        lock.release()

        # write down progress
        filtered_reads.value += 1

        # print progress
        if report_every > 0 and filtered_reads.value % report_every == 0:
            if lock is not None:
                lock.acquire()
            s = f'\r    processed: {filtered_reads.value:8d} (passed filtering: {" ".join(map(str, annotated_reads))}) {datetime.datetime.now() - start_time}'
            report.log_str(s, stdout_too=False, priority=logging.DEBUG, flush=True)
            sys.stdout.write(s)
            sys.stdout.flush()
            if lock is not None:
                lock.release()

        # return result
        return annotations

    def annotate_read(read):
        """
        Annotate a single read with global Annotator object.
        :param read: triple - read to be annotated
        :return: tuple - annotated read (output of Annotator), read
        """
        # start profiler
        if config['general']['profiler']:
            global prof
            prof.enable()

        global lock
        global annotated_reads
        global filtered_reads
        global annotators
        global filters
        global filter_on
        global report_every

        try:
            annotations = annotate_read_sequentially(read, annotated_reads, filtered_reads, annotators, filters, filter_on, report_every, lock)
        except Exception as e:
            print('Exception in child process:')
            traceback.print_exc()
            raise e

        # stop profiler.
        if config['general']['profiler']:
            prof.disable()

        # return result
        return annotations

    def gather_res(results: typing.Iterable[typing.List[annotation.Annotation]], length: int) -> typing.List[typing.List[annotation.Annotation]]:
        """
        Filter results for None, transpose and return list.
        :param results: iterator(list) - Annotations as iterator of lists
        :param length: int - number of the annotators = number of the lists returned
        :return: list(list) - return Annotations in list of lists
        """
        res = [[] for _ in range(length)]
        for partial_res in results:
            for i, pr in enumerate(partial_res):
                if pr is not None:
                    res[i].append(pr)

        return res

    # crate pool and initialize it with init_pool
    if threads > 1 or config['general'].get('force_parallel', False):
        # print('Running in parallel ({threads} cores)'.format(threads=threads))
        pool = multiprocess.Pool(threads, initializer=init_pool,
                                 initargs=(lock, annotated_reads, filtered_reads, annotators, filters, filter_on, report_every))
        results = pool.imap(annotate_read, reader, chunksize=100)
        res = gather_res(results, len(annotators))
        pool.close()
        pool.join()
    else:
        # print('Running sequentially')
        results = (annotate_read_sequentially(read, annotated_reads, filtered_reads, annotators, filters, filter_on, report_every) for read in reader)
        res = gather_res(results, len(annotators))

    # adjust for report every:
    if report_every > 0:
        print('')

    # return
    return res, filtered_reads.value, annotated_reads
