import os
import cProfile
import pstats
from datetime import datetime

import arguments
from parser import ReadFile
import prefiltering
import templates
import report
import annotation
import postfilter
from annotation import Annotator
import all_call
from parallelism.run_parallel import run_annot_iter, MOTIF_PRINT_LEN, shorten_str, create_annotators


if __name__ == '__main__':
    # print time of the start:
    start_time = datetime.now()
    print(templates.START_TIME.format(start=start_time))

    # load arguments
    config = arguments.load_arguments()

    # start profiling if needed
    prog_prof = None
    if config['general']['profiler']:
        prog_prof = cProfile.Profile()
        prog_prof.enable()

    # initialize logging module:
    report.configure_logger('%s/dante.log' % config['general']['output_dir'])
    report.log_str(templates.START_TIME.format(start=start_time), stdout_too=False)

    # create all annotators and filters
    annotators, filters = create_annotators(config)

    # main annotation
    dedup_ap = [[] for _ in range(len(config['motifs']))]
    all_reads = 0
    annotations = [[] for _ in range(len(annotators))]
    annotated_reads = [0] * len(annotations)
    readers = []

    # read reads from specific region in bam?, read unmapped reads?
    bam = True
    include_unmapped = False

    for motif in config['motifs']:
        if 'chromosome' not in motif:
            bam = False
            include_unmapped = False
            break
        elif 'include_unmapped' in motif and motif['include_unmapped']:
            include_unmapped = True

    for input_file in config['inputs']:
        # initialize reader
        read_filename = input_file['path']
        file_type = None if input_file['filetype'] == 'infer' else input_file['filetype']
        max_reads = None if input_file['max_reads'] == 'all' else input_file['max_reads']

        if bam:
            # annotate specified reads for every motif
            report.log_str(f'Parsing file {read_filename} (bam type)')
            for i, motif in enumerate(config['motifs']):
                if config['general']['quiet_mode'] and i % 100 == 0:
                    report.log_str(
                        'Running annotators: {i:6d}/{cnt:6d}   {date:%Y-%m-%d %H:%M:%S}'.format(i=i, cnt=len(config['motifs']), date=datetime.now()))
                min_mapq = None if 'min_mapq' not in motif else motif['min_mapq']
                read_file = ReadFile(read_filename, config['general']['stranded'], None if max_reads is None else max_reads - all_reads, file_type,
                                     config['general']['verbosity'], motif['chromosome'], motif['ref_start'], motif['ref_end'], min_mapq=min_mapq)

                # initialize annotator and filter for current motif
                annot = [annotators[i]]
                filt = [filters[i]]
                areads = [sum(annotated_reads)]
                readers.append(read_file)

                # run annotator for current motif
                report.log_str(
                    f'Running annotator on {config["general"]["cpu"]} process(es) for file {read_filename} and motif {motif["description"]}',
                    stdout_too=not config['general']['quiet_mode'])
                next_annotations, cur_reads, areads_new = run_annot_iter(read_file, annot, filt, config['general']['cpu'], areads, all_reads,
                                                                         config, start_time)
                new_reads = cur_reads - all_reads
                all_reads += new_reads
                annotated_reads[i] = areads_new[0] - areads[0]
                annotations[i].extend(next_annotations[0])

                # write stats for current file and motif
                report.log_str('Reads: %10d, annotated: %s' % (new_reads, ' '.join(map(str, map(len, next_annotations)))),
                               stdout_too=not config['general']['quiet_mode'])

        if not bam or include_unmapped:
            read_file = ReadFile(None if read_filename == 'sys.stdin' else read_filename, config['general']['stranded'], max_reads, file_type,
                                 verbosity=config['general']['verbosity'], unmapped=include_unmapped)

            report.log_str('Parsing file %s (%s type) with parser %s %s ' % (read_filename, read_file.file_type, read_file.__class__.__name__,
                                                                             read_file.reader.__name__))

            # run annotators
            report.log_str('Running %d annotator(s) on %2d process(es) for file %s' % (len(annotators), config['general']['cpu'], read_filename))
            readers.append(read_file)
            next_annotations, cur_reads, annotated_reads = run_annot_iter(read_file, annotators, filters, config['general']['cpu'], annotated_reads,
                                                                          all_reads, config, start_time)
            new_reads = cur_reads - all_reads
            all_reads += new_reads
            for i, pr in enumerate(next_annotations):
                annotations[i].extend(pr)

            # write stats for current file
            report.log_str('Reads: %10d, annotated: %s' % (new_reads, ' '.join(map(str, map(len, next_annotations)))),
                           stdout_too=not config['general']['quiet_mode'])

    # write read distribution
    report.write_read_distribution('%s/read_distr.npy' % config['general']['output_dir'], readers)

    # write stats
    for i, (motif, annot) in enumerate(zip(config['motifs'], annotations)):

        report.log_str(f'Motif {shorten_str(motif["full_name"], MOTIF_PRINT_LEN):<{MOTIF_PRINT_LEN}s}: Kept {len(annot):8d}/{all_reads:8d} reads',
                       stdout_too=not config['general']['quiet_mode'])

        # setup motif sequences and dir
        motif_dir = '%s/%s' % (config['general']['output_dir'], motif['full_name'].replace('/', '_'))
        sequence = ','.join([module_dict['seq'] for module_dict in motif['modules']])
        motif_tuple = report.seq_into_tuple(sequence)

        # convert to annotation pairs, deduplicate and convert back
        annotation_pairs = annotation.annotations_to_pairs(annot)
        dedup_ap[i], duplicates = annotation.remove_pcr_duplicates(annotation_pairs)

        # report them
        if not config['general']['quiet_mode']:
            if not os.path.exists(motif_dir):
                os.makedirs(motif_dir)

            report.write_annotation_pairs('%s/annotation_pairs.txt' % motif_dir, dedup_ap[i])
            report.write_annotation_pairs('%s/annotation_pairs_duplicates.txt' % motif_dir, duplicates)

        # log it
        report.log_str(f'Motif {shorten_str(motif["full_name"], MOTIF_PRINT_LEN):<{MOTIF_PRINT_LEN}s}: {len(annot):8d} reads -- '
                       f'{len(annotation_pairs):8d} pairs ({len(dedup_ap[i]):8d} deduplicated + {len(duplicates):8d} PCR duplicates)',
                       stdout_too=not config['general']['quiet_mode'])

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
            if not config['general']['quiet_mode']:
                report.log_str(f'Motif {shorten_str(motif["full_name"], MOTIF_PRINT_LEN):<{MOTIF_PRINT_LEN}s}: Extracted {len(dedup_annot):8d} '
                               f'reads from {len(dedup_ap[i]):8d} pairs', stdout_too=not config['general']['quiet_mode'])

            # get filtered stuff
            report.log_str(f'Motif {shorten_str(motif["full_name"], MOTIF_PRINT_LEN):<{MOTIF_PRINT_LEN}s}: Running post-filtering for '
                           f'{pstf["repetitions"]} required repetitions and {pstf["bases"]} required bases',
                           stdout_too=not config['general']['quiet_mode'])
            qual_annot, primer_annot, filt_annot = postfilter_class.get_filtered(dedup_annot)
            report.log_str(f'Motif {shorten_str(motif["full_name"], MOTIF_PRINT_LEN):<{MOTIF_PRINT_LEN}s}: Post-filtered '
                           f'{len(filt_annot):5d} (at least one primer {len(primer_annot):5d}), remained {len(qual_annot):5d}',
                           stdout_too=not config['general']['quiet_mode'])

            # write it to files
            report.log_str(f'Motif {shorten_str(motif["full_name"], MOTIF_PRINT_LEN):<{MOTIF_PRINT_LEN}s}: Generating output files into {motif_dir}',
                           stdout_too=not config['general']['quiet_mode'])
            report.write_all(qual_annot, primer_annot, filt_annot, dedup_ap[i], all_reads, motif_dir, motif['modules'], index_rep, index_rep2, j,
                             config['general']['quiet_mode'])

    # -------- All_Call part of DANTE

    # run all_call
    for i, motif in enumerate(config['motifs']):

        sequence = ','.join([module_dict['seq'] for module_dict in motif['modules']])
        motif_tuple = report.seq_into_tuple(sequence)
        motif_dir = '%s/%s' % (config['general']['output_dir'], motif['full_name'].replace('/', '_'))

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
            qual_annot, primer_annot, _ = postfilter_class.get_filtered(dedup_annot)
            postfilter_bases = postfilter_class.get_filtering_bases()

            # load read distribution
            read_distribution = report.load_read_distribution('%s/read_distr.npy' % config['general']['output_dir'])

            # run inference
            inference = all_call.Inference(read_distribution, config['allcall']['param_file'], str_rep=len_str,
                                           minl_primer1=postfilter_bases[index_rep - 2], minl_primer2=postfilter_bases[index_rep],
                                           minl_str=postfilter_bases[index_rep - 1])
            file_pcolor = '%s/pcolor_%d' % (motif_dir, j + 1)
            if config['general']['quiet_mode']:
                file_pcolor = None
            file_output = '%s/allcall_%d.txt' % (motif_dir, j + 1)
            inference.all_call(qual_annot, primer_annot, index_rep - 1, file_pcolor, file_output, motif['full_name'])

            # write the report
            confidence = report.read_all_call('%s/allcall_%d.txt' % (motif_dir, j + 1))
            if confidence is not None:
                conf, a1, a2, c1, c2, _, _, _, _ = confidence
                if not config['general']['quiet_mode']:
                    if isinstance(a1, int) and a1 > 0:
                        report.write_alignment('%s/alignment_%d_a%d.fasta' % (motif_dir, j + 1, a1), qual_annot, index_rep - 1, allele=a1)
                    if isinstance(a2, int) and a2 != a1 and a2 != 0:
                        report.write_alignment('%s/alignment_%d_a%d.fasta' % (motif_dir, j + 1, a2), qual_annot, index_rep - 1, allele=a2)

    # -------- generation of reports and finalizing

    # generate report and output files for whole run
    report.log_str('Generating final report')
    report.write_report(config['general']['output_dir'], config['motifs'], config['general']['output_dir'], config['general']['quiet_mode'],
                        config['general']['skip_annotations'])

    # print the time of the end:
    end_time = datetime.now()
    report.log_str('DANTE Stopping    : {finish:%Y-%m-%d %H:%M:%S}'.format(finish=end_time))
    report.log_str('Total time of run : {duration}'.format(duration=end_time - start_time))

    # stop profiler:
    if config['general']['profiler']:
        prog_prof.disable()
        with open('%s/profile-main.out' % config['general']['output_dir'], 'w') as f:
            pstats.Stats(prog_prof, stream=f).strip_dirs().sort_stats('time').print_stats()
