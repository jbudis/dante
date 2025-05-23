import os
import re
import shutil
import typing
import gzip
import matplotlib
from collections import Counter
from matplotlib.colors import ListedColormap
import plotly.graph_objects as go
import numpy as np
import pandas as pd

import parser

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import annotation
import templates
import report.html_templates

# max repetitions on the graph
MAX_REPETITIONS = 40


def gen_annot_string(annot: annotation.Annotation) -> str:
    """
    Get the annotation description string
    :param annot: Annotation - to be described
    ::return: str - annotation description
    """
    likelihood = float(annot.probability)
    pairness = {None: 'Unknown', True: 'R1', False: 'R2'}[annot.read.left_pair]
    ann = templates.ANNOTATION.format(annotation=annot, likelihood=likelihood,
                                      errors=annot.n_insertions + annot.n_deletions + annot.n_mismatches, pairness=pairness)
    return ann


def write_annotations(out_file: str, annotations: typing.List[annotation.Annotation], zip_it: bool = True) -> None:
    """
    Stores annotations in alignment format into output text file
    :param out_file: Alignment file
    :param annotations: Annotated reads
    :param zip_it: bool - whether to gzip the resulting file
    """
    if zip_it and not out_file.endswith('.gz'):
        out_file += '.gz'
    with gzip.open(out_file, 'wt') if zip_it else open(out_file, 'w') as fw:
        for annot in annotations:
            fw.write(gen_annot_string(annot))


def write_annotation_pairs(out_file: str, annotation_pairs: typing.List[annotation.AnnotationPair], zip_it: bool = True) -> None:
    """
    Stores annotations in alignment format into output text file
    :param out_file: Alignment file
    :param annotation_pairs: Annotated pairs of reads
    :param zip_it: bool - whether to gzip the resulting file
    """
    if zip_it and not out_file.endswith('.gz'):
        out_file += '.gz'
    with gzip.open(out_file, 'wt') if zip_it else open(out_file, 'w') as fw:
        for ap in annotation_pairs:
            write_left = gen_annot_string(ap.ann1) if ap.ann1 is not None else 'Left None\n'
            fw.write(write_left)
            write_right = gen_annot_string(ap.ann2) if ap.ann2 is not None else 'Right None\n'
            fw.write(write_right + '\n')


def write_alignment(out_file: str, annotations: typing.List[annotation.Annotation], index_rep: int, index_rep2: int = None, allele: int = None,
                    allele2: int = None, zip_it: bool = True) -> None:
    """
    Creates a multi-alignment of all annotations into output text file
    :param out_file: str - alignment filename
    :param annotations: list(Annotation) - annotated reads
    :param index_rep: int - index of repetition module of a motif
    :param index_rep2: int - index of second repetition module of a motif
    :param allele: int/None - which allele to print only, if None print all of them
    :param allele2: int/None - which allele to print only from the second position, if None print all of them
    :param zip_it: bool - whether to gzip the resulting file
    """
    if allele is not None:

        if allele2 is not None and index_rep2 is not None:
            annotations = [a for a in annotations if a.module_repetitions[index_rep] == allele and a.module_repetitions[index_rep2] == allele2]
        else:
            annotations = [a for a in annotations if a.module_repetitions[index_rep] == allele]

    alignments = [''] * len(annotations)
    align_inds = np.zeros(len(annotations), dtype=int)
    states = []
    inserts = []

    while True:
        # get minimal state:
        min_comp = (True, 1e9, -1)
        total_done = 0
        for i, (annot, ai) in enumerate(zip(annotations, align_inds)):
            if ai >= len(annot.states):
                total_done += 1
                continue
            state = annot.states[ai]
            comparator = (not state.is_insert(), state.state_id, i)
            min_comp = min(comparator, min_comp)

        # if we have done every state, end:
        if total_done >= len(alignments):
            break

        states.append(min_comp[1])
        inserts.append('_' if min_comp[0] else 'I')

        # now print all states, that are minimal:
        for i, (annot, ai) in enumerate(zip(annotations, align_inds)):
            if ai >= len(annot.states):
                alignments[i] += '_'
                continue
            state = annot.states[ai]
            if state.state_id == min_comp[1]:
                alignments[i] += annot.read.sequence[ai]
                align_inds[i] += 1
            else:
                alignments[i] += '_'

    # sort according to motif count:
    if index_rep2 is not None:
        # trick to have sorting first with 1st allele then with second
        reps = np.array([ann.module_repetitions[index_rep] * 100000 + ann.module_repetitions[index_rep2] for ann in annotations])
    else:
        reps = np.array([ann.module_repetitions[index_rep] for ann in annotations])
    sort_inds = np.argsort(-reps)
    annotations = np.array(annotations)[sort_inds]
    alignments = np.array(alignments)[sort_inds]

    # move 0 states (first background) to the right:
    def move_right(alignment: str, bckg_length: int) -> str:
        """
        Shift first part of the alignment to the right.
        :param alignment: str - alignment of the read
        :param bckg_length: int - length of the alignment, that is generated by the first background state
        :return: str - alignment, where first part is shifted to right
        """
        background = alignment[:bckg_length]
        # find last empty:
        idx = 0
        for idx in reversed(range(bckg_length)):
            if background[idx] != '_':
                break
        idx += 1

        # return shifted alignment
        return ('_' * (bckg_length - idx)) + background[:idx] + alignment[bckg_length:]

    states = np.array(states)
    background_length = np.sum(states == 0)
    if background_length > 0:
        for i in range(len(alignments)):
            alignments[i] = move_right(alignments[i], background_length)

    def get_extended_num(num: int) -> str:
        """
        Convert to extended number (A = 10, B = 11, ...)
        :param num: int - number to be converted
        :return: str - the converted string representation of the number
        """
        if num < 10:
            return str(num)
        else:
            return chr(ord('A') + num - 10)

    # print to file
    if zip_it and not out_file.endswith('.gz'):
        out_file += '.gz'
    with gzip.open(out_file, 'wt') if zip_it else open(out_file, 'w') as fw:
        # print alignments
        for annot, align in zip(annotations, alignments):
            print('>%s' % annot.read.name, file=fw)
            print(align, file=fw)
        # print debug info
        print('# debug', file=fw)
        print(''.join(inserts), file=fw)
        print(''.join(map(lambda x: get_extended_num(x // 10), states)), file=fw)
        print(''.join(map(lambda x: str(x % 10), states)), file=fw)


def write_summary_statistics(out_file: str, annotations: typing.List[annotation.Annotation], n_parsed_reads: int) -> None:
    """
    Stores annotation statistics into output text file.
    :param out_file: File for output statistics
    :param annotations: Annotated reads
    :param n_parsed_reads: Number of reads in the input read file
    """
    remained = n_parsed_reads - len(annotations)
    inserts = sum(annot.n_insertions for annot in annotations)
    deletes = sum(annot.n_deletions for annot in annotations)
    bases = sum(annot.n_bases for annot in annotations)
    with open(out_file, 'w') as fw:
        fw.write(templates.OUT_STATS.format(inserts=inserts, deletes=deletes, bases=bases, remained=remained))


def sorted_repetitions(annotations: typing.List[annotation.Annotation]) -> typing.List[typing.Tuple[typing.Tuple[int, ...], int]]:
    """
    Aggregate same repetition counts for annotations and sort them according to quantity of repetitions of each module
    :param annotations: Annotated reads
    :return: list of (repetitions, count), sorted by repetitions
    """
    count_dict = Counter(tuple(annot.module_repetitions) for annot in annotations)
    return sorted(count_dict.items(), key=lambda k: k[0])


def write_histogram(out_file: str, annotations: typing.List[annotation.Annotation], profile_file: str = None,
                    index_rep: int = None, quiet: bool = False) -> None:
    """
    Stores quantity of different combinations of module repetitions into text file
    :param out_file: str - output file for repetitions
    :param annotations: Annotated reads
    :param profile_file: str - output file for profile
    :param index_rep: int - index repetition
    :param quiet: boolean - write profile?
    """
    # setup
    sorted_reps = sorted_repetitions(annotations)

    # write repetitions.txt
    with open(out_file, 'w') as fw:
        for repetitions, counts in sorted_reps:
            rep_code = '\t'.join(map(str, repetitions))
            fw.write('%s\t%s\n' % (counts, rep_code))

    # write profile if possible and needed
    if not quiet and profile_file is not None and index_rep is not None:
        length = max([0] + [x[0][index_rep] for x in sorted_reps])
        profile = np.zeros(length + 1, dtype=int)

        for repetitions, counts in sorted_reps:
            profile[repetitions[index_rep]] += counts

        with open(profile_file, 'w') as f:
            f.write('\t'.join(map(str, profile)))


def write_histogram_nomenclature(out_file: str, annotations: typing.List[annotation.Annotation],
                                 index_rep: int = None, index_rep2: int = None) -> None:
    """
    Stores quantity of different nomenclature strings into text file
    :param out_file: str - output file for repetitions
    :param annotations: Annotated reads
    :param index_rep: int - index of the first repetition (None if include all)
    :param index_rep2: int - index of the second repetition (None if include all)
    """
    # count nomenclature strings:
    count_dict = Counter(annot.get_nomenclature(index_rep, index_rep2, False) for annot in annotations)
    count_dict = sorted(count_dict.items(), key=lambda k: (-k[1], k[0]))

    # write nomenclatures to file
    with open(out_file, 'w') as fw:
        for nomenclature, count in count_dict:
            fw.write('%s\t%s\n' % (count, nomenclature))


def write_histogram_image2d(out_prefix: str, deduplicated: typing.List[annotation.Annotation],
                            index_rep: int, index_rep2: int, seq: str, seq2: str) -> None:
    """
    Stores quantity of different combinations of module repetitions, generates separate graph image for each module
    :param out_prefix: Output file prefix
    :param deduplicated: list[Annotation] - read pairs
    :param index_rep: int - index of repetition module of a motif
    :param index_rep2: int - index of the second repetition module of a motif
    :param seq: str - module of the repetition
    :param seq2: str - 2nd module of the repetition
    """
    if deduplicated is None or len(deduplicated) == 0:
        return

    dedup_reps = [(x.get_str_repetitions(index_rep), x.get_str_repetitions(index_rep2)) for x in deduplicated
                  if x.get_str_repetitions(index_rep) is not None and x.get_str_repetitions(index_rep2) is not None]

    if len(dedup_reps) == 0:
        return

    # assign maximals
    xm = max([r for (_, r), _ in dedup_reps])
    ym = max([r for _, (_, r) in dedup_reps])
    max_ticks = max(ym, xm) + 2
    xm = max(MAX_REPETITIONS, xm)
    ym = max(MAX_REPETITIONS, ym)

    # create data containers
    data = np.zeros((xm + 1, ym + 1), dtype=int)
    data_primer = np.zeros((xm + 1, ym + 1), dtype=int)
    for (c1, r1), (c2, r2) in dedup_reps:
        if c1 and c2:
            data[r1, r2] += 1
        if c1 and not c2:
            data_primer[r1, r2:] += 1
        if not c1 and c2:
            data_primer[r1:, r2] += 1

    # create colormap:
    cmap = matplotlib.cm.Reds
    my_cmap = cmap(np.arange(int(cmap.N * 0.15), int(cmap.N * 0.9)))  # start from orange to deep red (but not the almost-black red)
    my_cmap[0, -1] = 0.0  # Set alpha on the lowest element only
    my_cmap = ListedColormap(my_cmap)

    # plot pcolor
    plt.figure(figsize=(12, 8))
    img2 = plt.pcolor(data_primer[:max_ticks, :max_ticks], cmap='Blues', alpha=0.4, edgecolor=(1.0, 1.0, 1.0, 0.0), lw=0, vmin=np.min(data_primer),
                      vmax=np.max(data_primer) + 0.01)
    img1 = plt.pcolor(data[:max_ticks, :max_ticks], cmap=my_cmap, vmin=np.min(data), vmax=np.max(data) + 0.01)
    plt.xticks()
    plt.ylabel('STR %d [%s]' % (index_rep + 1, seq.split('-')[-1]))
    plt.xlabel('STR %d [%s]' % (index_rep2 + 1, seq2.split('-')[-1]))
    plt.colorbar(img1)
    plt.colorbar(img2)

    # setup ticks
    start_ticks = 5
    step_ticks = 5
    plt.xticks(np.array(range(start_ticks, max_ticks + 1, step_ticks)) + 0.5, range(start_ticks, max_ticks + 1, step_ticks))
    plt.yticks(np.array(range(start_ticks, max_ticks + 1, step_ticks)) + 0.5, range(start_ticks, max_ticks + 1, step_ticks))

    # output it
    plt.savefig(out_prefix + '.pdf')
    plt.savefig(out_prefix + '.png')
    plt.close()

    # ----- PLOTLY HISTOGRAM -----
    def parse_labels(num, num_primer):
        if num == 0 and num_primer == 0:
            return ''
        elif num == 0 and num_primer != 0:
            return '0 / %s' % str(num_primer)
        elif num != 0 and num_primer == 0:
            return '%s / 0' % str(num)
        else:
            return '%s/%s' % (str(num), str(num_primer))

    str1 = 'STR %d [%s]' % (index_rep + 1, seq.split('-')[-1])
    str2 = 'STR %d [%s]' % (index_rep2 + 1, seq2.split('-')[-1])

    text = [[parse_labels(data[i, j], data_primer[i, j]) for j in range(data.shape[1])] for i in range(data.shape[0])]

    fig = go.Figure()
    fig.add_trace(go.Heatmap(z=data_primer[:max_ticks, :max_ticks], name='Repetitions heatmap',
                             showscale=True, colorbar_x=1.3, colorbar_title='Partial reads', colorscale='Blues'))
    fig.add_trace(go.Heatmap(z=data[:max_ticks, :max_ticks], text=text, name='Repetitions heatmap',
                             showscale=True, colorbar_title='Full reads',
                             colorscale=[[0.0, 'rgba(255, 255, 255, 0.0)'],
                                         [0.01, 'rgba(249, 217, 201, 1.0)'],
                                         [0.5, 'rgba(229, 103, 76, 1.0)'],
                                         [1.0, 'rgba(143, 33, 29, 1.0)']]))

    fig.update_traces(texttemplate='%{text}', textfont_size=7,
                      hovertemplate='<b>{name1}:\t%{y}<br />{name2}:\t%{x}</b><br />Full / Partial:\t%{text}'.
                      format(name1=str1, y='{y}', name2=str2, x='{x}', text='{text}'))
    fig.update_layout(width=800, height=700, template='simple_white')
    fig.update_yaxes(title_text=str1)
    fig.update_xaxes(title_text=str2)

    with open(out_prefix + '.json', 'w') as f:
        f.write(fig.to_json())

    # fig.write_image(out_prefix + '_plotly.pdf')


def write_histogram_image(out_prefix: str, annotations: typing.List[annotation.Annotation],
                          filt_annot: typing.List[annotation.Annotation], index_rep: int) -> None:
    """
    Stores quantity of different combinations of module repetitions, generates separate graph image for each module
    :param out_prefix: Output file prefix
    :param annotations: Annotated reads.
    :param filt_annot: Annotated reads (filtered)
    :param index_rep: int - index of repetition module of a motif
    """
    if len(annotations) == 0 and len(filt_annot) == 0:
        return

    repetitions = sorted_repetitions(annotations)
    repetitions_filt = sorted_repetitions(filt_annot)

    # adjust variables
    width = 0.9
    plt.figure(figsize=(20, 8))
    xm = max([r[index_rep] for r, _ in repetitions] + [MAX_REPETITIONS])
    if repetitions_filt:
        xm = max([r[index_rep] for r, _ in repetitions_filt] + [xm])
    dist = [0] * (xm + 1)
    dist_filt = [0] * (xm + 1)

    # set data
    for r, c in repetitions:
        dist[r[index_rep]] += c
        dist_filt[r[index_rep]] += c
    for r, c in repetitions_filt:
        dist_filt[r[index_rep]] += c

    # create barplots
    rects_filt = plt.bar(np.arange(xm + 1), dist_filt, width, color='grey', alpha=0.4)
    rects = plt.bar(np.arange(xm + 1), dist, width)
    plt.xticks(np.arange(1, xm + 1))
    plt.ylabel('Counts')
    plt.xlabel('STR repetitions')
    _, max_y = plt.ylim()
    plt.xlim((0, xm + 1))

    # label numbers
    for rect, rect_filt in zip(rects, rects_filt):
        height = rect.get_height()
        height_filt = rect_filt.get_height()

        if height > 0:
            plt.text(rect.get_x() + rect.get_width() / 2., height + max_y / 100.0, '%d' % int(height), ha='center', va='bottom')
        if height_filt != height:
            plt.text(rect_filt.get_x() + rect_filt.get_width() / 2., height_filt + max_y / 100.0, '%d' % int(height_filt - height), ha='center',
                     va='bottom', color='grey')

    # output it
    plt.savefig(out_prefix + '.pdf')
    plt.savefig(out_prefix + '.png')
    plt.close()

    # ----- PLOTLY HISTOGRAM -----
    def parse_labels(_dist_filt, _dist):
        if _dist_filt == 0:
            return ''
        elif _dist_filt != 0 and _dist == 0:
            return str(_dist_filt)
        elif _dist_filt != 0 and _dist != 0:
            if _dist_filt > _dist:
                return str(_dist_filt - _dist)
            else:
                return ''
        else:
            return str(_dist_filt)

    dist_text = ["" if d == 0 else str(d) for d in dist]
    dist_filt_text = [parse_labels(df, d) for df, d in zip(dist_filt, dist)]

    fig = go.Figure()
    fig.add_bar(y=dist_filt, text=dist_filt_text, name='Partial reads',
                marker_color='rgb(204, 204, 204)', textfont_color='rgb(204, 204, 204)')
    fig.add_bar(y=dist, text=dist_text, marker_color='#636EFA', name='Full reads ')

    fig.update_traces(textposition='outside', texttemplate='%{text}', hovertemplate='%{text}', textfont_size=7)
    fig.update_layout(width=800, height=450,
                      title='Histogram of repetitions',
                      hovermode='x',
                      yaxis_fixedrange=True,
                      template='simple_white',
                      barmode='overlay')
    fig.update_yaxes(title_text='Read counts')
    fig.update_xaxes(title_text='STR repetitions', tickmode='array',
                     tickvals=list(range(5, len(dist), 5)),
                     ticktext=list(range(5, len(dist), 5)))

    with open(out_prefix + '.json', 'w') as f:
        f.write(fig.to_json())

    # fig.write_image(out_prefix + '_plotly.pdf')


def write_all(quality_annotations: typing.List[annotation.Annotation], filt_primer: typing.List[annotation.Annotation],
              filtered_annotations: typing.List[annotation.Annotation], all_reads: int, motif_dir: str, motif_modules: dict, index_rep: int,
              index_rep2: int, j: int, quiet: bool = False, zip_it: bool = True) -> None:
    """
    Write all output files: quality annotations, one-primer annotations, filtered annotations, statistics, repetitions + images.
    :param quality_annotations: list(Annotation) - list of blue annotations
    :param filt_primer: list(Annotation) - list of grey annotations
    :param filtered_annotations: list(Annotation) - list of filtered out annotation
    :param all_reads: int - number of all reads
    :param motif_dir: str - path to motif directory
    :param motif_modules: dict - motif modules dictionary from config
    :param index_rep: int - index of first repetition in modules
    :param index_rep2: int - index of second repetition in modules
    :param j: int - index of post-filter
    :param quiet: boolean - fewer files on the output?
    :param zip_it: bool - whether to gzip the resulting file
    :return: None
    """
    # create dir if not exists:
    if (not quiet or len(quality_annotations) > 0 or len(filt_primer) > 0) and not os.path.exists(motif_dir):
        os.makedirs(motif_dir)

    # write output files
    if not quiet:
        write_annotations('%s/annotations_%d.txt' % (motif_dir, j + 1), quality_annotations, zip_it=zip_it)
        write_annotations('%s/filtered_%d.txt' % (motif_dir, j + 1), filtered_annotations, zip_it=zip_it)
        write_annotations('%s/filtered_primer_%d.txt' % (motif_dir, j + 1), filt_primer, zip_it=zip_it)
        write_summary_statistics('%s/stats_%d.txt' % (motif_dir, j + 1), quality_annotations, all_reads)
        write_alignment('%s/alignment_%d.fasta' % (motif_dir, j + 1), quality_annotations, index_rep, index_rep2, zip_it=zip_it)
        write_alignment('%s/alignment_filtered_%d.fasta' % (motif_dir, j + 1), filt_primer, index_rep, index_rep2, zip_it=zip_it)

        if index_rep2 is not None:
            write_histogram_image2d('%s/repetitions_%d' % (motif_dir, j + 1), quality_annotations + filt_primer, index_rep, index_rep2,
                                    motif_modules[index_rep]['seq'], motif_modules[index_rep2]['seq'])
        else:
            write_histogram_image('%s/repetitions_%d' % (motif_dir, j + 1), quality_annotations, filt_primer, index_rep)

    if not quiet or len(quality_annotations) > 0:
        write_histogram('%s/repetitions_%d.txt' % (motif_dir, j + 1), quality_annotations, profile_file='%s/profile_%d.txt' % (motif_dir, j + 1),
                        index_rep=index_rep, quiet=quiet)
        write_histogram_nomenclature('%s/nomenclatures_%d.txt' % (motif_dir, j + 1), quality_annotations, index_rep=index_rep, index_rep2=index_rep2)
    if not quiet or len(filt_primer) > 0:
        write_histogram('%s/repetitions_grey_%d.txt' % (motif_dir, j + 1), filt_primer, quiet=quiet)
        write_histogram_nomenclature('%s/nomenclatures_grey_%d.txt' % (motif_dir, j + 1), filt_primer, index_rep=index_rep, index_rep2=index_rep2)


def get_seq_from_module(module_dict: dict) -> str:
    """
    Get one line sequence of a motif from its parameters.
    :param module_dict: dict - motif parameters
    :return: str - sequence of a motif
    """
    return ','.join([m['seq'] for m in module_dict])


def read_all_call(allcall_file: str) -> (float, int, int, float, float, float, float, float, float):
    """
    Read AllCall output and returns allele predictions and confidences.
    :param allcall_file: str - filename of AllCall output
    :return: tuple(5) - overall_confidence, allele numbers, and confidences for them, 0 for allele number if BG is the best
    """
    if not os.path.exists(allcall_file):
        return None

    with open(allcall_file) as f:
        lines = f.readlines()

    overall_conf = float(lines[0].strip().split()[-1].split('%')[0]) / 100

    def get_allele(line):
        """
        Get allele number and its confidence from a line.
        :param line: str - a single line of AllCall output
        :return: (int/str, float) - allele number and its confidence
        """
        split = line.strip().split()
        num = split[0]
        conf = float(split[-1].split('%')[0]) / 100
        try:
            num = int(num)
        except ValueError:
            pass
        return num, conf

    def get_probability(line):
        perc = line.strip().split()[-1][:-1]
        return float(perc)

    a1, c1 = get_allele(lines[1])
    a2, c2 = get_allele(lines[2])
    c3 = get_probability(lines[3])
    c4 = get_probability(lines[4])
    c5 = get_probability(lines[5])
    c6 = get_probability(lines[6])

    return overall_conf, a1, a2, c1, c2, c3, c4, c5, c6


def custom_format(template: str, **kwargs) -> str:
    """
    Custom format of strings for only those that we provide
    :param template: str - string to format
    :param kwargs: dict - dictionary of strings to replace
    :return: str - formatted string
    """
    for k, v in kwargs.items():
        template = template.replace('{%s}' % k, v)

    return template


def seq_into_tuple(sequence: str) -> typing.List[typing.Tuple[str, int]]:
    """
    Extract all repeating motifs.
    :param sequence: str - motif sequence in "ACATCAG,3-TGT,CATCGACT" format
    :return: list(tuple) - motif sequence in list((seq, rep)) format
    """
    modules = sequence.split(',')
    return [(m.split('-')[-1], 1 if len(m.split('-')) == 1 else int(m.split('-')[0])) for m in modules]


def tuple_into_seq(list_tuple: typing.List[typing.Tuple[str, int]]) -> str:
    """
    Extract all repeating motifs.
    :param list_tuple: list(tuple) - motif sequence in list((seq, rep)) format
    :return: str - motif sequence in "ACATCAG,3-TGT,CATCGACT" format
    """
    return ','.join([s if r == 1 else '%d-%s' % (r, s) for (s, r) in list_tuple])


def get_read_count(filename: str) -> int:
    """
    Get read count from repetitions filename.
    :param filename: str - filename
    :return: int - number of reads
    """
    # get reads count:
    if not os.path.exists(filename):
        return 0

    reads = 0
    with open(filename) as f:
        for line in f:
            reads += int(line.split()[0])

    return reads


# add to report table
def add_to_result_table(result_table: pd.DataFrame, motif_name: str, seq: str, postfilter: dict, reads_blue: int,
                        reads_grey: int, confidence: (float, float)) -> pd.DataFrame:
    """
    Create report table in pandas.
    :param result_table: pandas.DataFrame - table with the results
    :param motif_name: str - motif name
    :param seq: str - sequence
    :param postfilter: dict - post-filter options
    :param reads_blue: int - number of full reads
    :param reads_grey: int - number of partial reads
    :param confidence: tuple - allele predictions and their confidence
    :return: pandas.DataFrame - table with all the results
    """
    # write the results into a table in TSV format
    result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Motif'] = motif_name
    result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Sequence'] = seq
    result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Repetition index'] = postfilter['index_rep']
    result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Postfilter bases'] = postfilter['bases']
    result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Postfilter repetitions'] = postfilter['repetitions']
    result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Reads (full)'] = reads_blue
    result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Reads (partial)'] = reads_grey

    if confidence is not None:
        result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Overall confidence'] = confidence[0]
        result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Allele 1 prediction'] = confidence[1]
        result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Allele 2 prediction'] = confidence[2]
        result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Allele 1 confidence'] = confidence[3]
        result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Allele 2 confidence'] = confidence[4]
        result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Both Background prob.'] = confidence[5]
        result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'One Background prob.'] = confidence[6]
        result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'Background Expanded prob.'] = confidence[7]
        result_table.at['%s_%s' % (motif_name, postfilter['index_rep']), 'One Expanded prob.'] = confidence[8]

    return result_table


def find_file(filename: str, include_gzip: bool = False) -> typing.Optional[str]:
    """
    Find if we have a file with the provided filename.
    :param filename: str - file name
    :param include_gzip: bool - try the gzipped suffix
    :return: str/None - filename if exists, None if not
    """
    if os.path.exists(filename):
        return filename
    gzipped = filename + '.gz'
    if include_gzip and os.path.exists(gzipped):
        return gzipped
    return None


def write_report(report_dir: str, motifs: dict, output_dir: str, nomenclature=5, quiet: bool = False, skip_alignments: bool = False, include_alignments: bool = False) -> None:
    """
    Generate and write a report.
    :param report_dir: str - dir name for reports
    :param motifs: dict - parameters of motifs
    :param output_dir: str - output directory for searching of all_call outputs
    :param nomenclature: int - number of lines from nomenclature.txt to print
    :param quiet: boolean - fewer files on the output?
    :param skip_alignments: boolean - skip annotation generation
    :return: None
    """
    # tsv file with table:
    result_table = pd.DataFrame([], columns=['Motif', 'Sequence', 'Repetition index', 'Postfilter bases', 'Postfilter repetitions',
                                             'Overall confidence', 'Allele 1 prediction', 'Allele 1 confidence', 'Allele 2 prediction',
                                             'Allele 2 confidence', 'Reads (full)', 'Reads (partial)', 'Both Background prob.',
                                             'One Background prob.', 'Background Expanded prob.', 'One Expanded prob.'])

    # merge all_profiles:
    all_profiles = '%s/all_profiles.txt' % output_dir
    all_true = '%s/all_profiles.true' % output_dir

    mcs = {}
    ms = {}
    rows = {}
    alignments = {}
    mcs_static = []
    ms_static = []
    rows_static = []

    with open(all_profiles, 'w') as pf, open(all_true, 'w') as tf:
        for m in motifs:
            seq = get_seq_from_module(m['modules'])
            motif_name = m['full_name'].replace('/', '_')
            description = m['description']
            for i, postfilter in enumerate(m['postfilter']):
                # read files
                rep_file = find_file('%s/%s/repetitions_%d.json' % (output_dir, motif_name, i + 1))
                pcol_file = find_file('%s/%s/pcolor_%d.json' % (output_dir, motif_name, i + 1))
                align_file = None if skip_alignments else find_file('%s/%s/alignment_%d.fasta' % (output_dir, motif_name, i + 1), include_gzip=True)
                filt_align_file = None if skip_alignments else find_file('%s/%s/alignment_filtered_%d.fasta' % (output_dir, motif_name, i + 1),
                                                                         include_gzip=True)
                confidence = read_all_call('%s/%s/allcall_%d.txt' % (output_dir, motif_name, i + 1))

                # get number of reads:
                reads_blue = get_read_count('%s/%s/repetitions_%d.txt' % (output_dir, motif_name, i + 1))
                reads_grey = get_read_count('%s/%s/repetitions_grey_%d.txt' % (output_dir, motif_name, i + 1))

                # generate rows of table and images
                highlight = [postfilter['index_rep'] - 1]
                if postfilter['index_rep2'] != 'no':
                    highlight.append(postfilter['index_rep2'] - 1)
                row = report.html_templates.generate_row(motif_name, seq, confidence, postfilter, reads_blue, reads_grey, highlight=highlight)
                if motif_name in rows:
                    rows[motif_name].append(row)
                else:
                    rows[motif_name] = [row]

                rows_static.append(row)

                # add to csv table:
                result_table = add_to_result_table(result_table, motif_name, seq, postfilter, reads_blue, reads_grey, confidence)

                if not quiet:
                    mc, m, a = report.html_templates.generate_motifb64(motif_name, description, seq, rep_file, pcol_file, align_file, filt_align_file, confidence, postfilter, highlight=highlight, include_alignments=include_alignments)
                    if motif_name in mcs:
                        ms[motif_name].append(m)
                        alignments[motif_name][1].append(a[1])
                    else:
                        mcs[motif_name] = mc
                        ms[motif_name] = [m]
                        alignments[motif_name] = (a[0], [a[1]])

                    mc, m, a = report.html_templates.generate_motifb64(motif_name, description, seq, rep_file, pcol_file, align_file, filt_align_file, confidence, postfilter, highlight=highlight, include_alignments=include_alignments, static=True)
                    if mc not in mcs_static:
                        mcs_static.append(mc)
                    ms_static.append(m)

                    # add to profiles
                    if postfilter['index_rep2'] == 'no':
                        with open('%s/%s/profile_%d.txt' % (output_dir, motif_name.replace('/', '_'), i + 1)) as po:
                            line = po.readline()
                            pf.write('%s_%d\t%s\n' % (motif_name, i + 1, line))

                        # add to true
                        if confidence is None:
                            confidence = [0.0, 0, 0]
                        tf.write('%s_%d\t%s\t%s\n' % (motif_name, i + 1, str(confidence[1]), str(confidence[2])))

    # save the report file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    template = open('%s/report.html' % script_dir, 'r').read()
    template_static = open('%s/report.html' % script_dir, 'r').read()

    sample = os.path.basename(output_dir)

    template = custom_format(template, sample=sample)
    template_static = custom_format(template_static, sample=sample)
    tabs = []

    with open('%s/report.html' % report_dir, 'w') as f:
        contents_table = report.html_templates.contents.format(table='\n'.join(sorted(mcs.values())))
        template = custom_format(template, motifs_content=contents_table + '\n' + report.html_templates.make_datatable_string)

        for motif in sorted(motifs, key=lambda x: x['full_name']):
            motif_name = motif['full_name'].replace('/', '_')
            motif_clean = re.sub(r'[^\w_]', '', motif_name)

            with open('%s/%s/nomenclature.txt' % (report_dir, motif_name), 'r') as noms:
                lines = []
                for line in noms:
                    if line == '' or line is None:
                        break

                    line_split = line.split('\t')
                    motif_parts = [f'<td>{s}</td>' for s in line_split[1:]]

                    # replace chromosome_version if not available
                    if 'chromosome_version' not in motif and 'chromosome' in motif:
                        motif['chromosome_version'] = motif['chromosome']
                    # build reference string
                    ref = ''
                    if 'chromosome_version' in motif and 'ref_start' in motif and 'ref_end' in motif:
                        ref = f'{motif["chromosome_version"]}:g.{motif["ref_start"]}_{motif["ref_end"]}'

                    nom_row = report.html_templates.nomenclature_string.format(count=line_split[0] + 'x', ref=ref, parts='\n    '.join(motif_parts))
                    lines.append(nom_row)

                    # end?
                    if len(lines) >= nomenclature:
                        break

            tabs.append(report.html_templates.motif_summary.format(motif_id=motif_clean,
                                                                   nomenclatures='\n'.join(lines), table='\n'.join(rows[motif_name]),
                                                                   motifs='\n'.join(ms[motif_name])))

        f.write(custom_format(template, table='', motifs='\n'.join(tabs)))

    with open('%s/report_static.html' % report_dir, 'w') as f:
        contents_table = report.html_templates.contents.format(table='\n'.join(mcs_static))
        table = report.html_templates.motif_summary_static.format(table='\n'.join(rows_static))

        f.write(custom_format(template_static, motifs_content=contents_table,
                              table=table, motifs='\n'.join(ms_static)))

    if not include_alignments:
        for motif in alignments.keys():
            template_alignments = open('%s/alignments.html' % script_dir, 'r').read()
            template_alignments = custom_format(template_alignments, sample=motif,
                                                motif_desc=alignments[motif][0])

            with open('%s/%s/alignments.html' % (report_dir, motif), 'w') as f:
                f.write(custom_format(template_alignments, alignments='\n'.join(alignments[motif][1])))

    # copy javascript libraries
    if not quiet:
        shutil.copy2('%s/msa.min.gz.js' % script_dir, '%s/msa.min.gz.js' % report_dir)
        shutil.copy2('%s/plotly-2.14.0.min.js' % script_dir, '%s/plotly-2.14.0.min.js' % report_dir)
        shutil.copy2('%s/jquery-3.6.1.min.js' % script_dir, '%s/jquery-3.6.1.min.js' % report_dir)
        shutil.copy2('%s/datatables.min.js' % script_dir, '%s/datatables.min.js' % report_dir)

    # save the table(s)
    result_table.to_csv('%s/table.tsv' % report_dir, sep='\t')


def write_read_distribution(report_file: str, readers: typing.List[parser.ReadFile]) -> float:
    """
    Save the read distribution generated from all annotators. Return mean read length.
    :param report_file: str - file to write the distribution to
    :param readers: list(ReadFile) - ReadFiles with stored sub-distributions
    :return: float - mean read length
    """
    # get maximal read and create the empty distribution
    max_length = 0
    all_distribs = [reader.distribution for reader in readers]
    for distribution in all_distribs:
        max_length = max(list(distribution.keys()) + [max_length])

    read_distribution = np.zeros(max_length + 1, dtype=int)

    # fill the distribution
    total_count = 0
    total_length = 0
    for reader in readers:
        for k, v in reader.distribution.items():
            read_distribution[k] += v
            total_count += v
            total_length += k * v

    # save the distribution
    np.save(report_file, read_distribution)

    # as figure, too
    plt.plot(read_distribution)
    plt.savefig(report_file + '.png')
    plt.savefig(report_file + '.pdf')
    plt.close()

    # return mean length
    return total_length / total_count


def load_read_distribution(report_file: str) -> np.ndarray:
    """
    Load the read_distribution.
    :param report_file: str - distribution filename
    :return: ndarray - read length distribution
    """
    return np.load(report_file)
