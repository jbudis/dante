from .DummyFilter import DummyFilter


class SimpleFilter(DummyFilter):
    """
    Filter class for simple filtering - if the all sequence are found in read, the read is accepted.
    """

    def __init__(self, filters):
        """
        Initialize the filter with target sequences and allowed deviation.
        :param filters: list(tuple) - target sequences and their repetitions
        """
        super(self.__class__, self).__init__()
        self.filters = filters

    def filter_read(self, read):
        """
        Test whether the read has passsed the filter.
        :param read: (3tuple) - the specified read.
        :return: bool - True if read has passed the filter
        """
        for seq, repetitions in self.filters:
            if seq * repetitions not in read.sequence:
                return False
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
        return 'SimpleFilter:' + ','.join(map(lambda x: x[0] * x[1], self.filters))
