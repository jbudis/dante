import numpy as np
from numba import jit
from prefiltering.DummyFilter import DummyFilter


@jit
def levenshtein_partial(target_str, search_str, maximal_dist=10):
    """
    Search for Levenshtein partial distance between target and search strings. (Searches search string in target string)
    :param target_str: str - target string
    :param search_str: str - string to search in target
    :param maximal_dist: int - maximal distance to search
    :return: int, 2D-ndarray - Levenshtein partial distance, dynamic programming matrix
    """

    # setup variables
    n = len(target_str)
    m = len(search_str)

    # setup array
    max_mn = max(m, n)
    dist = np.full((n + 1, m + 1), max_mn, dtype=int)
    dist[:, 0] = np.zeros(n + 1)
    dist[0, :] = np.arange(m + 1)

    # dynamic programming search
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if dist[i - 1, j - 1] > maximal_dist:
                continue
            add = 1 if search_str[j - 1] != target_str[i - 1] else 0
            dist[i, j] = min(dist[i - 1, j] + 1, dist[i, j - 1] + 1, dist[i - 1, j - 1] + add)

    # return distance and matrix
    return np.min(dist[:, m]), dist


@jit
def levenshtein_backward(target_str, search_str, dist):
    """
    Backward pass to obtain best match from levenshtein_partial search. Slow - use only for debug.
    :param target_str: str - target string
    :param search_str: str - string to search in target
    :param dist: 2D-ndarray - dynamic programming matrix from levenshtein_partial
    :return: str - best match string ("_" means deletions from search string, lowercase means insertions to search string)
    """

    m = len(search_str)
    pos = np.argmin(dist[:, m])
    if dist[pos, m] == 0: return search_str, pos - len(search_str)

    res = ""
    i = m
    while i > 0:
        print(res[::-1])
        minimal = np.argmin([dist[pos, i], dist[pos, i - 1], dist[pos - 1, i]])
        if minimal == 0:
            pos -= 1
            i -= 1
            res += target_str[pos]
        elif minimal == 1:
            i -= 1
            res += "_"
        elif minimal == 2:
            pos -= 1
            res += target_str[pos].lower()

    return res[::-1]


class LevenshteinFilter(DummyFilter):
    """
    Filter class for Levenshtein partial filtering. All partial strings must be found with some insertions/deletions/substitutions allowed.
    """

    def __init__(self, filters, fuzziness):
        """
        Initialize the filter with target sequence and allowed deviation.
        :param filters: str - target sequence
        :param fuzziness: int/list(int) in string - allowed deviation
        """
        super(self.__class__, self).__init__()
        self.filters = filters
        self.fuzziness = fuzziness
        if len(self.fuzziness) == 1:
            self.fuzziness *= len(self.filters)
        assert len(self.filters) == len(self.fuzziness), "Number of filters should correspond to number of deviations"

    def filter_read(self, read):
        """
        Test whether the read has passed the filter.
        :param read: (3tuple) - the specified read.
        :return: bool - True if read has passed the filter
        """
        for (seq, repetitions), max_dev in zip(self.filters, self.fuzziness):
            deviations, _ = levenshtein_partial(read.sequence, seq * repetitions, maximal_dist=max_dev)
            if deviations > max_dev:
                return False
        return True

    def encapsulate_reader(self, reader):
        """
        Encapsulates reader object with fuzzy filter.
        :param reader: iterator over sequences
        :return: reader with fuzzy filter
        """
        for read in reader:
            if self.filter_read(read):
                yield read

    def __str__(self):
        return 'LevenshteinFilter: fuzziness=(%s) ' % (','.join(map(str, self.fuziness))) + ','.join(map(lambda x: x[0] * x[1], self.filters))
