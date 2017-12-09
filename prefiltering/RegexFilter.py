import regex
from DummyFilter import DummyFilter


class RegexFilter(DummyFilter):
    """
    Filter class for regex filtering.
    """

    def __init__(self, filters, mismatches):
        """
        Initialize the filter with target sequences and allowed deviation.
        :param filters: list(tuple) - target sequences and their repetitions
        :param mismatches: int - allowed deviation
        """
        super(self.__class__, self).__init__()
        self.filters = filters
        self.mismatches = int(mismatches)

        # build reference
        reference = "".join(["%s.*" % (seq * rep) for (seq, rep) in self.filters])
        reference = "(.*%s){e<%d}" % (reference, self.mismatches)

        self.regex_compiled = regex.compile(reference, regex.BESTMATCH)

    def filter_read(self, read):
        """
        Test whether the read has passed the filter.
        :param read: (3tuple) - the specified read.
        :return: bool - True if read has passed the filter
        """
        match = self.regex_compiled.match(read.sequence)
        return match is not None

    def encapsulate_reader(self, reader):
        """
        Encapsulates reader object with regex filter.
        :param reader: iterator over sequences
        :return: reader with regex filter
        """
        for read in reader:
            if self.filter_read(read):
                yield read

    def __str__(self):
        return 'RegexFilter: mismatches=%d ' % self.mismatches + ','.join(map(lambda x: x[0] * x[1], self.filters))
