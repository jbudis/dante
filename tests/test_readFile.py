from unittest import TestCase
from parser import ReadFile
import os
from glob import glob

N_READS = 10
READ_FILE_PREFIX = 'n10'
READ_1_SEQ = 'CCGAGTAGTTGGGATTACAGGCATGTGCCACCATGCCCCAGCTAATTTTGTATTTTTAGTTAGAGATGGGGTTTCTCTGTGCTGGTCAGGCTGGTCTCGAACTCCTGACCTCAGGTGAGCTCGAGGATCACCGACTGCCCATAGAGAGGCTGAGACTGCCAAGGC'
READ_2_SEQ = 'TCGCCCCTCACTCAC'
RC_1_SEQ = 'GCCTTGGCAGTCTCAGCCTCTCTATGGGCAGTCGGTGATCCTCGAGCTCACCTGAGGTCAGGAGTTCGAGACCAGCCTGACCAGCACAGAGAAACCCCATCTCTAACTAAAAATACAAAATTAGCTGGGGCATGGTGGCACATGCCTGTAATCCCAACTACTCGG'
RC_2_SEQ = 'GTGAGTGAGGGGCGA'
READ_1_QUAL = 'D@FDBA:::4::16C3:<DDAFDEEB?>2::2::::::/:999@<@@@43833399-93928C>>>8///*/91777///--88788877:00*/7::78-98=:>>?B8@=22-----------(-5---2666>>?E7>?BBB999;A?==AB8822-22825'
READ_2_QUAL = '+.....(556666;;'
RC_1_QUAL = '52822-2288BA==?A;999BBB?>7E?>>6662---5-(-----------22=@8B?>>:=89-87::7/*00:77888788--///77719/*///8>>>C82939-99333834@@@<@999:/::::::2::2>?BEEDFADD<:3C61::4:::ABDF@D'
RC_2_QUAL = ';;666655(.....+'
MISSING_1_QUAL = 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
MISSING_2_QUAL = 'IIIIIIIIIIIIIII'
MISSING_RC_1_QUAL = MISSING_1_QUAL
MISSING_RC_2_QUAL = MISSING_2_QUAL


class TestReadFile(TestCase):
    """
    Unit tests for reading of files.
    """

    read_file_dir = '%s/read_files' % os.path.dirname(__file__)

    def iter_test_files(self):
        """
        Iterator for test case files.
        :return: iterator(str) over filenames
        """
        for rf in glob('%s/%s.*' % (self.read_file_dir, READ_FILE_PREFIX)):
            yield rf

    def test_n_reads(self):
        """
        test if number of reads in test files is N_READS.
        :return: None
        """

        for read_file in self.iter_test_files():
            print 'Processing file %s' % read_file
            reader = ReadFile(read_file, 'yes')
            read_count = 0
            for _ in reader.iter_reads():
                read_count += 1

            self.assertEqual(read_count, N_READS,
                             'Number of parsed reads in %s does not match' % read_file)

    def test_parsed_sequences(self):
        """
        Test if the first sequences in test files are exactly READ_X_SEQ.
        :return: None
        """
        for read_file in self.iter_test_files():
            print 'Processing file %s' % read_file
            reader = ReadFile(read_file, 'yes')
            iterator = reader.iter_reads()
            _, s1, _ = iterator.next()
            self.assertEqual(s1, READ_1_SEQ,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, s1, READ_1_SEQ))
            _, s2, _ = iterator.next()
            self.assertEqual(s2, READ_2_SEQ,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, s2, READ_2_SEQ))

    def test_parsed_quality(self):
        """
        Test if the first qualities in test files are exactly READ_X_QUAL/MISSING_X_QUAL.
        :return: None
        """
        for read_file in self.iter_test_files():
            print 'Processing file %s' % read_file
            has_quality = not (read_file.endswith('.fa') or read_file.endswith('.fa.gz') or
                               read_file.endswith('.fasta') or read_file.endswith('.fasta.gz'))
            reader = ReadFile(read_file, 'yes')
            r1_qual = READ_1_QUAL if has_quality else MISSING_1_QUAL
            r2_qual = READ_2_QUAL if has_quality else MISSING_2_QUAL
            iterator = reader.iter_reads()
            _, _, q1 = iterator.next()
            self.assertEqual(q1, r1_qual,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, q1, r1_qual))
            _, _, q2 = iterator.next()
            self.assertEqual(q2, r2_qual,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, q2, r2_qual))

    def test_stranded_reverse(self):
        """
        Test if the first qualities and sequences in reverse read test files are RC_X_QUAL/MISSING_RC_X_QUAL.
        :return: None
        """
        for read_file in self.iter_test_files():
            print 'Processing file %s' % read_file
            has_quality = not (read_file.endswith('.fa') or read_file.endswith('.fa.gz') or
                               read_file.endswith('.fasta') or read_file.endswith('.fasta.gz'))
            reader = ReadFile(read_file, 'reverse')
            r1_qual = RC_1_QUAL if has_quality else MISSING_RC_1_QUAL
            r2_qual = RC_2_QUAL if has_quality else MISSING_RC_2_QUAL
            iterator = reader.iter_reads()
            _, s1, q1 = iterator.next()
            self.assertEqual(q1, r1_qual,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, q1, r1_qual))
            self.assertEqual(s1, RC_1_SEQ,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, s1, RC_1_SEQ))
            _, s2, q2 = iterator.next()
            self.assertEqual(q2, r2_qual,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, q2, r2_qual))
            self.assertEqual(s2, RC_2_SEQ,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, s2, RC_2_SEQ))

    def test_stranded_both(self):
        """
        Test if the first qualities and sequences in both reverse and forward read test files are RC_X_QUAL/MISSING_RC_X_QUAL and READ_X_QUAL/MISSING_X_QUAL.
        :return: None
        """
        for read_file in self.iter_test_files():
            print 'Processing file %s' % read_file
            has_quality = not (read_file.endswith('.fa') or read_file.endswith('.fa.gz') or
                               read_file.endswith('.fasta') or read_file.endswith('.fasta.gz'))
            reader = ReadFile(read_file, 'both')
            r1_rc_qual = RC_1_QUAL if has_quality else MISSING_RC_1_QUAL
            r2_rc_qual = RC_2_QUAL if has_quality else MISSING_RC_2_QUAL
            r1_qual = READ_1_QUAL if has_quality else MISSING_1_QUAL
            r2_qual = READ_2_QUAL if has_quality else MISSING_2_QUAL
            iterator = reader.iter_reads()
            _, s1, q1 = iterator.next()
            self.assertEqual(q1, r1_qual,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, q1, r1_qual))
            self.assertEqual(s1, READ_1_SEQ,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, s1, READ_1_SEQ))

            _, s1_rc, q1_rc = iterator.next()
            self.assertEqual(q1_rc, r1_rc_qual,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, q1_rc, r1_rc_qual))
            self.assertEqual(s1_rc, RC_1_SEQ,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, s1_rc, RC_1_SEQ))

            _, s2, q2 = iterator.next()
            self.assertEqual(q2, r2_qual,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, q2, r2_qual))
            self.assertEqual(s2, READ_2_SEQ,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, s2, READ_2_SEQ))

            _, s2_rc, q2_rc = iterator.next()
            self.assertEqual(q2_rc, r2_rc_qual,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, q2_rc, r2_rc_qual))
            self.assertEqual(s2_rc, RC_2_SEQ,
                             'Sequence of parsed reads in %s does not match %s, %s' % (read_file, s2_rc, RC_2_SEQ))
