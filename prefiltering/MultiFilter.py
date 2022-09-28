from prefiltering.DummyFilter import DummyFilter


class MultiFilter(DummyFilter):
    """
    Filter class for filtering with multiple filters
    """

    def __init__(self, filters):
        """
        Initialize BamFilter object
        :param filters - list(*Filter) - filters to accept
        """
        super(self.__class__, self).__init__()
        self.filters = filters
        assert len(self.filters) > 0

    def filter_read(self, read):
        """
        Test whether the read has passed the filter.
        :param read: Read - the specified read.
        :return: bool - True if read has passed the filter
        """
        passed = True
        for filter in self.filters:
            if not filter.filter_read(read):
                passed = False
                break

        return passed

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
        return 'MultiFilter: ' + ' + '.join(map(str, self.filters))
