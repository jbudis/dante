import re

from prefiltering.DummyFilter import DummyFilter

# define base mapping to regex for nucleotide symbols
base_mapping = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': '[GA]',
    'Y': '[CT]',
    'K': '[GT]',
    'M': '[AC]',
    'S': '[GC]',
    'W': '[AT]',
    'D': '[GAT]',
    'H': '[ACT]',
    'V': '[GCA]',
    'N': '[ACTG]'
}


class SimpleFilter(DummyFilter):
    """
    Filter class for simple filtering - if the all sequence are found in read, the read is accepted
    """

    def __init__(self, filters):
        """
        Initialize the filter with target sequences and allowed deviation
        :param filters: list(tuple) - target sequences and their repetitions
        """
        super(self.__class__, self).__init__()
        self.filters = [seq * repetitions for seq, repetitions in filters]

        # mode of operation - use simpler, quicker solution if only ACTG are present
        self.simple = all(all(ch in 'ACTG' for ch in test_string) for test_string in self.filters)

        # if mode is not simple, build regex filters
        self.regex_filters = None
        if not self.simple:
            self.regex_filters = [re.compile(r''.join(base_mapping.get(ch, '[ACTG]') for ch in filter)) for filter in self.filters]

    def filter_read(self, read):
        """
        Test whether the read has passed the filter
        :param read: (3tuple) - the specified read
        :return: bool - True if read has passed the filter
        """
        if self.simple:
            for test_string in self.filters:
                if test_string in read.sequence:
                    return True
        else:
            for regex_compiled in self.regex_filters:
                if regex_compiled.search(read.sequence):
                    return True

        return False

    def encapsulate_reader(self, reader):
        """
        Encapsulates reader object with dummy filter
        :param reader: iterator over sequences
        :return: reader with dummy filter
        """
        for read in reader:
            if self.filter_read(read):
                yield read

    def __str__(self):
        return 'SimpleFilter:' + ','.join(self.filters)
