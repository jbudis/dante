class DummyFilter(object):
    """
    Filter class for dummy filtering (no sequences are to be filtered out)
    """

    def __init__(self):
        pass

    def filter_read(self, read):
        """
        Test whether the read has passed the filter.
        :param read: (3tuple) - the specified read.
        :return: bool - True if read has passed the filter
        """
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
        return 'DummyFilter'
