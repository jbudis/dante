from annotation.Annotation import Annotation
from operator import attrgetter


def annotations_to_pairs(annotations):
    """
    Convert an array of annotations to annotation pairs array.
    :param annotations: list(Annotation) - annotations
    :return: list(AnnotationPair)
    """

    # sort:
    sorted_list = sorted(annotations, key=lambda ann: (ann.read.sort_name, ann.read.is_left_current()))

    # remove duplicates:
    seen = set()
    deduplicated = []
    for ann in sorted_list:
        if (ann.read.sort_name, ann.read.is_left_current()) not in seen:
            seen.add((ann.read.sort_name, ann.read.is_left_current()))
            deduplicated.append(ann)

    result = []

    i = 0
    while i < len(deduplicated):
        if i + 1 < len(deduplicated) and deduplicated[i].read.sort_name == deduplicated[i + 1].read.sort_name:
            assert deduplicated[i].read.is_left_current() is not None and deduplicated[i + 1].read.is_left_current() is not None and \
                   deduplicated[i].read.is_left_current() != deduplicated[i + 1].read.is_left_current()
            result.append(AnnotationPair(deduplicated[i], deduplicated[i + 1]))
            i += 2
        else:
            result.append(AnnotationPair(deduplicated[i], None))
            i += 1

    return result


def pairs_to_annotations(annotation_pairs):
    """
    Convert an array of annotations pairs to annotation array.
    :param annotation_pairs: list(AnnotationPair) - annotations
    :return: list(Annotation)
    """
    annotations = []
    for ap in annotation_pairs:
        if ap.ann1 is not None:
            annotations.append(ap.ann1)
        if ap.ann2 is not None:
            annotations.append(ap.ann2)

    return annotations


def pairs_to_annotations_pick(annotation_pairs, index_str):
    """
    Convert an array of annotations pairs to annotation array. Leave only the more informative one.
    :param annotation_pairs: list(AnnotationPair) - annotations
    :param index_str: int: index of the STR to look at
    :return: list(Annotation)
    """
    annotations = []

    for ap in annotation_pairs:
        if ap.ann1 is None:
            annotations.append(ap.ann2)
            continue
        if ap.ann2 is None:
            annotations.append(ap.ann1)
            continue
        ann1_comparison = (ap.ann1.primers(index_str), ap.ann1.module_repetitions[index_str], ap.ann1.module_bases[index_str])
        ann2_comparison = (ap.ann2.primers(index_str), ap.ann2.module_repetitions[index_str], ap.ann2.module_bases[index_str])
        if ann1_comparison > ann2_comparison:
            annotations.append(ap.ann1)
        else:
            annotations.append(ap.ann2)

    return annotations


def remove_pcr_duplicates(annot_pairs):  # TODO make this faster/parallelize - takes too long when number of found reads is more than 10.000-100.000 (1min at 7.500, 4h at 400.000)
    """
    Remove PCR duplicates -- deduplicate the annotation pair list.
    :param annot_pairs: list(AnnotationPair) - list of Annotation Pairs
    :return: list(AnnotationPair), list(AnnotationPair) - deduplicated list and duplications
    """
    def remove_none(ann_pairs: list, first: bool = True) -> list:
        """
        Remove annotation pairs with None first (second) annotation from the pair.
        :param ann_pairs: list(AnnotationPair) - list of Annotation Pairs
        :param first: bool - look at the first annotation from the pair?
        :return: list(AnnotationPair) - list of Annotation Pairs without None pairs
        """
        arr = []
        for ap in ann_pairs:
            if first and ap.ann1 is not None:
                arr.append(ap)
            if not first and ap.ann2 is not None:
                arr.append(ap)
        return arr

    def restore_none(pairs_with_none: list, pairs_without_none: list, first: bool = True) -> list:
        """
        Restore annotation pairs with None first (second) annotation from the pair.
        :param pairs_with_none: list(AnnotationPair) - list of Annotation Pairs with None annotations
        :param pairs_without_none: list(AnnotationPair) - list of Annotation Pairs without None annotations
        :param first: bool - look at the first annotation from the pair?
        :return: list(AnnotationPair) - list of Annotation Pairs with restored Annotation Pairs with None annotation
        """
        ann_pairs = pairs_without_none.copy()

        for ap in pairs_with_none:
            if first and ap.ann1 is None:
                ann_pairs.append(ap)
            if not first and ap.ann2 is None:
                ann_pairs.append(ap)
        return ann_pairs

    def deduplicate(ann_pairs: list):
        """
        Remove PCR duplicates -- deduplicate the annotation pair list.
        :param ann_pairs: list(AnnotationPair) - list of Annotation Pairs sorted by annotation_1 or annotation_2
        :return: list(AnnotationPair), list(AnnotationPair) - deduplicated list and duplications
        """
        dedup = []
        duplic = []

        if not ann_pairs:
            return [], []

        # Find duplicates by comparing neighbours in sorted list
        prev_ap = ann_pairs[0]
        for curr_ap in ann_pairs[1:]:
            if prev_ap == curr_ap:
                if prev_ap.more_info_than(curr_ap):
                    duplic.append(curr_ap)
                else:
                    duplic.append(prev_ap)
                    prev_ap = curr_ap
            else:
                dedup.append(prev_ap)
                prev_ap = curr_ap

        dedup.append(prev_ap)

        return dedup, duplic

    if not annot_pairs:
        return [], []

    # Deduplication according to first annotation in pair
    curr_pairs = remove_none(annot_pairs, True)
    curr_pairs = sorted(curr_pairs, key=attrgetter('ann1.read.sequence'))
    deduplicated_1, duplicates_1 = deduplicate(curr_pairs)

    curr_pairs = restore_none(annot_pairs, deduplicated_1, True)

    # Deduplication according to second annotation in pair
    curr_pairs = remove_none(curr_pairs, False)
    curr_pairs = sorted(curr_pairs, key=lambda ann: ann.ann2.read.sequence[::-1])
    deduplicated_2, duplicates_2 = deduplicate(curr_pairs)

    deduplicated = restore_none(deduplicated_1, deduplicated_2, False)
    duplicates = duplicates_1 + duplicates_2

    return deduplicated, duplicates


class AnnotationPair:
    """
    Encapsulate annotation pairs.
    """

    def __init__(self, ann1, ann2):
        """
        Initialize the AnnotationPair object
        :param ann1: Annotation - first annotation of a pair
        :type ann1: Annotation
        :param ann2: Annotation - second annotation of a pair
        :type ann2: Annotation | None
        """
        # assert isinstance(ann1, Annotation), "Ann1 is no Annotation, but %s" % type(ann1).__name__ # TODO why this throws?

        self.ann1 = ann1
        self.ann2 = ann2
        self.contains_none = ann2 is None

        if not self.contains_none:
            assert self.ann1.read.is_left_current() != self.ann2.read.is_left_current()

        if not self.ann1.read.is_left_current():
            self.ann1, self.ann2 = self.ann2, self.ann1

    def __eq__(self, second_pair):
        """
        Annotation pairs equal when they are produced by the same fragment.
        :param second_pair: AnnotationPair - second annotation pair
        :type second_pair: AnnotationPair
        :return: bool - whether the annotation pairs are produced by the same fragment
        """
        # first check if we deal with simple annotations
        if self.ann1 is None:
            return second_pair.ann2 is not None and self.ann2.same_end_fragment(second_pair.ann2)

        if self.ann2 is None:
            return second_pair.ann1 is not None and self.ann1.same_start_fragment(second_pair.ann1)

        # return the full comparison
        return (second_pair.ann1 is None or self.ann1.same_start_fragment(second_pair.ann1)) and (second_pair.ann2 is None or self.ann2.same_end_fragment(second_pair.ann2))

    def has_required_modules(self, required_repetitions):
        """
        Validate, if read modules are sufficiently annotated (in at least one of the reads)
        :param required_repetitions: list of number of required annotated modules
        :return: True, if one of the annotations has required number of annotated modules
        """
        left = self.ann1 is not None and self.ann1.has_required_modules(required_repetitions)
        rght = self.ann2 is not None and self.ann2.has_required_modules(required_repetitions)
        return left or rght

    def has_required_bases(self, required_bases):
        """
        Validate, if read bases are sufficiently annotated (in at least one of the reads)
        :param required_bases: list of number of required annotated bases, one for each module
        :return: True, if one of the annotations has required number of annotated bases for each module
        """
        left = self.ann1 is not None and self.ann1.has_required_bases(required_bases)
        rght = self.ann2 is not None and self.ann2.has_required_bases(required_bases)
        return left or rght

    def has_one_primer(self, required_bases, required_repetitions, index_rep, index_rep2=None):
        """
        Validate, if at least one primer is sufficiently annotated (in at least one of the reads)
        :param required_bases: list of number of required annotated bases, one for each module
        :param required_repetitions: list of number of required annotated modules
        :param index_rep: int - index of the repetition, that we are looking at
        :param index_rep2: int - index of the second repetition, that we are looking at
        :return: True, if one of the annotations has at least one primer is sufficiently annotated
        """
        left = self.ann1 is not None and self.ann1.has_one_primer(required_bases, required_repetitions, index_rep, index_rep2)
        rght = self.ann2 is not None and self.ann2.has_one_primer(required_bases, required_repetitions, index_rep, index_rep2)
        return left or rght

    def more_info_than(self, second_pair):
        """
        Check if this AnnotationPair has more info than second AnnotationPair
        :param AnnotationPair - second annotation pair
        :return: bool - True if it has
        """

        # first compare if both has annotations:
        if self.contains_none and not second_pair.contains_none:
            return False

        if not self.contains_none and second_pair.contains_none:
            return True

        # then return those that have more annotated modules, or bases:

        def eval_info_value(ann: Annotation):
            """
            Evaluate the info value of an Annotation.
            :param ann: Annotation - to evaluate
            :return: int, int
            """
            return 0 if ann is None else sum(ann.module_repetitions), 0 if ann is None else sum(ann.module_bases)

        m1f, b1f = eval_info_value(self.ann1)
        m2f, b2f = eval_info_value(self.ann2)

        m1s, b1s = eval_info_value(second_pair.ann1)
        m2s, b2s = eval_info_value(second_pair.ann2)

        return (m1f + m2f, b1f + b2f) > (m1s + m2s, b1s + b2s)

    def get_str_repetitions(self, index_str):
        """
        Get the number of str repetitions for a particular index.
        :param index_str: int - index of a str
        :return: (bool, int) - closed?, number of str repetitions
        """

        def add_primers(annotation, index_str):
            """
            Add str repetitions from one annotation into an array of results.
            :param annotation: Annotation - input annotation
            :param index_str: int - index of a str
            :return: list(bool, int) - closed?, number of str repetitions
            """
            primer1 = index_str > 0 and annotation.module_repetitions[index_str - 1] > 0
            primer2 = index_str + 1 < len(annotation.module_repetitions) and annotation.module_repetitions[index_str + 1] > 0
            if primer1 or primer2:
                return [(primer1 and primer2, annotation.module_repetitions[index_str])]
            return []

        results = []
        if self.ann1 is not None and self.ann1.is_annotated_right():
            results.extend(add_primers(self.ann1, index_str))
        if self.ann2 is not None and self.ann2.is_annotated_right():
            results.extend(add_primers(self.ann2, index_str))

        # return the highest (first see if it is closed and then pick the highest number):
        if len(results) > 0:
            return sorted(results)[-1]
        return None
