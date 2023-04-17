import arguments
import report


class Postfilter:
    """
    Class that encapsulated postfiltering.
    """

    def __init__(self, postfilter, motif_tuple):
        """
        Initialize postfilter class.
        :param postfilter: dict - postfilter dictionary - part of config file that tells about current postfiltering
        :param motif_tuple: tuple - motif sequence and its repetitions
        """
        self.postfilter = postfilter
        self.postfilter_bases = None if postfilter['bases'] == 'no' else list(map(int, postfilter['bases'].split(',')))
        self.postfilter_repetitions = None if postfilter['repetitions'] == 'no' else list(map(int, postfilter['repetitions'].split(',')))
        self.postfilter_errors = None if postfilter['max_errors'] == 'no' else float(postfilter['max_errors'])
        self.pf = self.postfilter_bases if self.postfilter_bases is not None else self.postfilter_repetitions
        self.index_rep = arguments.infer_index_rep(motif_tuple, self.pf) - 1 \
            if postfilter['index_rep'] == 'infer' else int(postfilter['index_rep']) - 1
        # is it normal or 2-index?
        self.index_rep2 = None
        if postfilter['index_rep2'] != 'no':
            self.index_rep2 = int(postfilter['index_rep2']) - 1

    def get_indexes(self):
        """
        Get indices of repetitions.
        :return: int, int - first and second index, second may be None
        """
        return self.index_rep, self.index_rep2

    def get_filtering_bases(self):
        """
        Get filtering bases.
        :return: list(int) - list of minimal numbers of bases to be present
        """
        return self.postfilter_bases

    def get_filtered(self, dedup_annot):
        """
        Get filtered annotations.
        :param dedup_annot: list(Annotation) - deduplicated annotations
        :return: list(Annotation), list(Annotation), list(Annotation) - blue, grey, and filtered out annotations
        """
        quality_annotations = [an for an in dedup_annot if an.has_required_modules(self.postfilter_repetitions) and
                                                     an.has_required_bases(self.postfilter_bases) and
                                                     an.has_less_errors(self.postfilter_errors, self.postfilter['max_errors_relative'])]
        filtered_annotations = [an for an in dedup_annot if an not in quality_annotations]
        filt_primer = [an for an in filtered_annotations if
                       an.has_one_primer(self.postfilter_bases, self.postfilter_repetitions, self.index_rep, self.index_rep2) and
                       an.has_less_errors(self.postfilter_errors, self.postfilter['max_errors_relative'])]
        filtered_annotations = [an for an in dedup_annot if an not in quality_annotations and an not in filt_primer]

        return quality_annotations, filt_primer, filtered_annotations
