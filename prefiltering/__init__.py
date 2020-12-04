from prefiltering.DummyFilter import DummyFilter
from prefiltering.RegexFilter import RegexFilter
from prefiltering.LevenshteinFilter import LevenshteinFilter
from prefiltering.SimpleFilter import SimpleFilter
from prefiltering.BamFilter import BamFilter
from prefiltering.MultiFilter import MultiFilter

import report


def get_prefilter_seq(prefilter, motif_tuple):
    """
    Get prefilter sequence from the prefilter config structure and motif tuple.
    :param prefilter: dict - prefilter config
    :param motif_tuple: tuple(str, int) - motif tuple (sequence, times)
    :return: tuple(str, int) - prefilter tuple (sequence, times)
    """
    prefilter_seq = motif_tuple if prefilter['seq'] == 'infer' else report.seq_into_tuple(prefilter['seq'])
    if prefilter['seq'] == 'infer' and prefilter['type'] == "SimpleFilter":
        # take only those that repeats and lower them by one
        repeating_seq = list(filter(lambda seq, rep: rep > 1, motif_tuple))
        prefilter_seq = list(map(lambda seq, rep: (seq, rep - 1), repeating_seq))
    return prefilter_seq


def create_filter(prefilter, motif_tuple):
    """
    Create the correct filter according to the info in the config file.
    :param prefilter: dict - prefilter config
    :param motif_tuple: tuple(str, int) - motif tuple (sequence, times)
    :return: *Filter - filter for reads
    """

    prefilter_seq = get_prefilter_seq(prefilter, motif_tuple)

    if prefilter['type'] == "RegexFilter":
        return RegexFilter(prefilter_seq, prefilter['mismatches'])
    elif prefilter['type'] == "LevenshteinFilter":
        return LevenshteinFilter(prefilter_seq, prefilter['mismatches'])
    elif prefilter['type'] == "SimpleFilter":
        return SimpleFilter(prefilter_seq)
    elif prefilter['type'] == "BamFilter":
        for item in ['chromosome', 'ref_start', 'ref_end']:
            assert item in prefilter, 'Error: %s is missing in BamFilter specification, add "%s: smth" line to config!' % (item, item)
        overlap = prefilter['overlap'] if 'overlap' in prefilter else 1
        min_mapq = prefilter['min_mapq'] if 'min_mapq' in prefilter else None
        return BamFilter(prefilter['chromosome'], prefilter['ref_start'], prefilter['ref_end'], overlap=overlap, min_mapq=min_mapq)
    elif prefilter['type'] == "MultiFilter":
        filters = list(map(lambda x: create_filter(x, motif_tuple), prefilter['subfilters']))
        return MultiFilter(filters)
    else:
        return DummyFilter()
