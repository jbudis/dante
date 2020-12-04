import annotation.Decoder
import numpy as np
import annotation.state
from annotation.Annotation import Annotation
import report

QUAL_BEG = 33
QUAL_END = 127
QUALITY_CODES = list(map(chr, range(QUAL_BEG, QUAL_END)))
QUALITY_NUMS = range(QUAL_END - QUAL_BEG)
NUCLEOTIDES = ['A', 'C', 'G', 'T', 'N']
MOTIF_NUCLEOTIDES = ['A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'X', 'N']
N_QUALITY_CODES = len(QUALITY_CODES)
N_NUCLEOTIDES = len(NUCLEOTIDES)

BASE_N_PROB = 0.001


def quality_index(asc):
    """
    Returns integer equivalent of quality (0-based)
    :param asc: char - quality of a single base
    :return: int - quality
    """
    return ord(asc) - QUAL_BEG


def code_base(base):
    """
    Changes nucleotide to integer code
    :param base: nucleotide , e.g C
    :return: integer code, e.g. 1
    """
    return NUCLEOTIDES.index(base)


class Annotator:
    """
    Class for annotation of reads.
    """

    def __init__(self, motif, delete_prob=0.01, insert_prob=0.01, max_delete_skips=2, motif_frequency=0.001, snp_chance=0.02):
        """
        Initializes the annotator with all of the constants and variables.
        :param motif: list(tuple) - motif to annotate tuple has (sequence, number of repeats)
        :param delete_prob: float - deletion probability
        :param insert_prob: float - insertion probability
        :param max_delete_skips: int - maximal skips in deletions
        :param motif_frequency: float - frequency of motif
        :param snp_chance: float - probability of SNP (single nucleotide polymorphism)
        """

        # Store argument parameters into instance variables
        self.motif = motif
        self.delete_prob = delete_prob
        self.insert_prob = insert_prob
        self.max_delete_skips = max_delete_skips + 1
        self.motif_frequency = motif_frequency
        self.n_modules = len(self.motif)

        self.initial = dict()
        self.trans = dict()

        # Sequence states that deletion edge skips
        self.deleted_states = dict()
        self.states = []
        self.repeats = []
        self.parts = []

        # construct background frequency
        background_prob = (1 - BASE_N_PROB) / 4.0
        self.background_frequency = {'A': background_prob, 'C': background_prob,
                                     'G': background_prob, 'T': background_prob,
                                     'N': BASE_N_PROB}

        def key_position_prob(quality, n_key_positions):
            """
            Calculates the position probability based on the quality.
            :param quality: int - quality number
            :param n_key_positions: int - number of preferred nucleotides at position
            :return: emission probability of base
            """

            def decode_quality_num(quality):
                """
                Converts quality letter from fastq file into quality of calling associated base
                See https://en.wikipedia.org/wiki/FASTQ_format#Variations for details
                :param quality: Quality code in one number 0-based
                :return: Probability of base to be correctly called
                """
                return 10.0 ** (-quality / 10.0)

            position_prob = (1 - snp_chance - decode_quality_num(quality) - BASE_N_PROB) / n_key_positions
            return max(background_prob, position_prob)

        # construct the emission probabilities: (sums to 1 for every quality)
        kp1 = {qual: key_position_prob(qual, 1) for qual in QUALITY_NUMS}
        kn3 = {qual: (1 - kp1[qual] - BASE_N_PROB) / 3.0 for qual in QUALITY_NUMS}

        kp2 = {qual: key_position_prob(qual, 2) for qual in QUALITY_NUMS}
        kn2 = {qual: (1 - kp2[qual] - BASE_N_PROB) / 2.0 for qual in QUALITY_NUMS}

        kp3 = {qual: key_position_prob(qual, 3) for qual in QUALITY_NUMS}
        kn1 = {qual: (1 - kp3[qual] - BASE_N_PROB) / 1.0 for qual in QUALITY_NUMS}

        nd = {qual: BASE_N_PROB for qual in QUALITY_NUMS}

        self.key_frequency = {'A': {'A': kp1, 'C': kn3, 'G': kn3, 'T': kn3, 'N': nd},
                              'C': {'A': kn3, 'C': kp1, 'G': kn3, 'T': kn3, 'N': nd},
                              'G': {'A': kn3, 'C': kn3, 'G': kp1, 'T': kn3, 'N': nd},
                              'T': {'A': kn3, 'C': kn3, 'G': kn3, 'T': kp1, 'N': nd},
                              'M': {'A': kp2, 'C': kp2, 'G': kn2, 'T': kn2, 'N': nd},
                              'R': {'A': kp2, 'C': kn2, 'G': kp2, 'T': kp2, 'N': nd},
                              'W': {'A': kp2, 'C': kn2, 'G': kn2, 'T': kp2, 'N': nd},
                              'S': {'A': kn2, 'C': kp2, 'G': kp2, 'T': kn2, 'N': nd},
                              'Y': {'A': kn2, 'C': kp2, 'G': kn2, 'T': kp2, 'N': nd},
                              'K': {'A': kn2, 'C': kn2, 'G': kp2, 'T': kp2, 'N': nd},
                              'V': {'A': kp3, 'C': kp3, 'G': kp3, 'T': kn1, 'N': nd},
                              'H': {'A': kp3, 'C': kp3, 'G': kn1, 'T': kp3, 'N': nd},
                              'D': {'A': kp3, 'C': kn1, 'G': kp3, 'T': kp3, 'N': nd},
                              'B': {'A': kn1, 'C': kp3, 'G': kp3, 'T': kp3, 'N': nd},
                              'X': self.background_frequency,
                              'N': self.background_frequency}

        assert set(self.key_frequency.keys()) == set(MOTIF_NUCLEOTIDES)

        self.__prepare_hmm()

    def __str__(self):
        """
        Prints his own motif.
        :return: str - motif that it annotates
        """
        return report.tuple_into_seq(self.motif)

    def __add_state(self, state):
        """
        Add new state into HMM
        :param state: name of the state
        """
        self.trans[state] = dict()
        self.states.append(state)

    def __add_sequence(self, sequence, repeats, module_id):
        """
        Prepare states and their emission probabilities for either primer or STR
        :param sequence: sequence of nucleotides
        :param repeats: number of sequence repetitions, 1 for primers
        """
        is_motif = repeats > 1
        state_constructor = annotation.state.MotifState if is_motif else annotation.state.SequenceState
        code = (len(self.states), len(self.states) + len(sequence) - 1, repeats)
        self.parts.append(code)
        if is_motif:
            self.repeats.append(code)
        for i, nucleotide in enumerate(sequence):
            # last = if it is not motif then True, else last - True if it is last.
            self.__add_state(state_constructor(nucleotide, self.key_frequency[nucleotide], module_id, len(self.states), i == 0 or not is_motif, i == len(sequence) - 1 or not is_motif))

    def __trans_probability(self, from_pos, to_pos):
        """
        Probability of transaction between two states
        :param from_pos: Index of state on the start of the transaction
        :param to_pos: Index of ending state
        :return: Probability of transaction, 0 if they are not connected
        """
        not_connected = 0
        from_state, to_state = self.states[from_pos], self.states[to_pos]
        if from_state not in self.trans:
            return not_connected
        if to_state not in self.trans[from_state]:
            return not_connected
        return self.trans[from_state][to_state]

    def __connect(self, from_pos, to_pos, prob, join=False):
        """
        Connect two states in HMM with an edge
        :param from_pos: name of the first state to connect
        :param to_pos: name of the second state to connect
        :param prob: transition probability of going from the first state to the second state
        :param join: if True, add probability to already existing edge probability, else use higher value
        """
        existing_prob = self.__trans_probability(from_pos, to_pos)
        self.trans[self.states[from_pos]][self.states[to_pos]] = prob + existing_prob if join else max(prob, existing_prob)

    """
    def __crange(self, motif_start, i, j, n):
        k = ((i - motif_start) + 1) % n
        ret = []
        while motif_start + k != j:
            ret.append(motif_start + k)
            k = (k + 1) % n
        return ret
    """

    # TODO verify this method, mainly skips inside cycles
    def __connect_deletions(self):
        """
        Prepare edges that can skip states in the linear state sequence or motifs, part of profile HMM extension
        """
        # linear skips
        for k in range(2, self.max_delete_skips + 1):
            for i in range(len(self.states) - k):  # skip to the k-th next state
                if i == 0:  # if we are connecting the module, we have to apply the motif frequency also
                    self.__connect(i, i + k, self.delete_prob ** (k - 1) * self.motif_frequency)
                else:
                    self.__connect(i, i + k, self.delete_prob ** (k - 1))
                self.deleted_states[(i, i + k)] = self.states[i + 1: i + k]

        # skips inside cycles
        for (motif_start, motif_end, motif_repetitions) in self.repeats:  # for all cycles
            motif_len = motif_end - motif_start + 1
            max_skip = min(self.max_delete_skips, motif_len - 1)
            for k in range(2, max_skip + 1):
                for i in range(0, motif_len):
                    del_start, del_end = motif_start + i, motif_start + (i + k) % motif_len
                    already_connected = (del_start, del_end) in self.deleted_states
                    if already_connected:
                        continue
                    self.__connect(del_start, del_end, self.delete_prob ** (k - 1))
                    self.deleted_states[(del_start, del_end)] \
                        = self.states[del_start + 1:motif_end + 1] + self.states[motif_start:del_end]

    # TODO validate this method
    def __connect_insertions(self):
        """
        Prepare states that can insert nucleotides between key positions, part of profile HMM extension
        """
        for key_state in range(1, len(self.states) - 2):
            insert_state = len(self.states)
            module_id = self.states[key_state].module_id
            self.__add_state(annotation.state.InsertState(self.background_frequency, module_id, len(self.states)))
            self.__connect(key_state, insert_state, self.insert_prob)  # go to insert state
            ends = [self.repeats[j][1] for j in range(len(self.repeats))]
            if key_state in ends:
                current_repeat_length = self.repeats[ends.index(key_state)][1] - self.repeats[ends.index(key_state)][0] + 1
                current_repeat_count = self.repeats[ends.index(key_state)][2]
                self.__connect(insert_state, insert_state, self.insert_prob)  # cycle in insert state
                self.__connect(insert_state, key_state - current_repeat_length + 1,
                               1 - self.insert_prob - 1 / current_repeat_count)  # go to the repeat beginning
                self.__connect(insert_state, key_state + 1, 1 / current_repeat_count)  # go to next repeat/primer
            else:
                self.__connect(insert_state, insert_state, self.insert_prob)  # cycle in insert state
                self.__connect(insert_state, key_state + 1, 1 - self.insert_prob)  # go to next state

    def __connect_next(self):
        """
        Create transactions between key states of the HMM
        """
        for state_idx, state in enumerate(self.states):

            # Insert states are already fully connected
            if state.is_insert():
                continue

            # Last motif state, do not further connect insert states
            next_state = self.states[state_idx + 1]
            if next_state.is_insert():
                break

            # remaining shift probability for a state (count the added insertions and deletions and other connections)
            remaining = 1 - sum(self.trans[state].values())
            assert remaining >= 0, 'Cannot connect %s-th state %s, already connected with probability %s' % \
                                   (state_idx, state, 1 - remaining)

            # already fully connected
            if remaining == 0:
                continue

            self.__connect(state_idx, state_idx + 1, remaining)

    def __connect_cycles(self):
        """
        Create cycles at the str motifs, ends with starts
        """
        for motif_start, motif_end, motif_repetitions in self.repeats:
            if motif_repetitions > 1:
                self.__connect(motif_end, motif_start, 1 - (1.0 / motif_repetitions))

    def __connect_modules(self):
        """
        Add skips of whole modules
        """
        ends = [0] + [end for _, end, _ in self.parts[:-1]]
        starts = [start for start, _, _ in self.parts] + [self.parts[-1][1] + 1]
        for i, end in enumerate(ends):
            remaining = 1 - sum(self.trans[self.states[end]].values())
            next_starts = starts[i + 1:]
            # dirty trick: we model first state differently: does not go into last background
            if i == 0:
                next_starts = starts[i + 1:-1]
            trans_prob = remaining / (len(next_starts) + 1)
            for start in next_starts:
                self.__connect(end, start, trans_prob)

    def __connect_first_state(self):
        """
        Creates loop at first background node
        """
        first_state = self.states[0]

        # Loop to itself - generate background random sequence of nucleotides
        remaining = 1 - sum(self.trans[first_state].values())
        self.__connect(0, 0, remaining - self.n_modules * self.motif_frequency)

    def __connect_states(self):
        """
        Add edges between states in the HMM
        """
        last_state = len(self.states) - 1

        # profile hmm specific connections
        self.__connect_deletions()
        self.__connect_insertions()

        # cycle at initial, background state
        self.__connect_first_state()

        # cycle at final, background state
        self.__connect(last_state, last_state, 1)

        # connect starts and ends of modules
        self.__connect_cycles()
        self.__connect_modules()

        # join subsequent states with remaining transition probabilities
        self.__connect_next()

    def __calculate_initial(self):
        """
        Assign to states probabilities of starting annotation in them
        """
        # key_states = [state for state in self.states if (state.is_key() and state.module_id != self.n_modules - 1) or state.is_background()]
        motif_states = [state for state in self.states if state.is_key()]
        backg_state = self.states[0]
        # All motif states should have same initial probability, background should have the rest
        for state in motif_states:
            self.initial[state] = self.motif_frequency
        # background takes the rest:
        self.initial[backg_state] = 1 - self.motif_frequency * len(motif_states)

    def __create_initial(self):
        """
        Transform initial probabilities into format suitable for the Decoder class
        :return: numpy array of initial probabilities (Nx1 float)
        """
        initial = np.array([self.initial.get(state, 0) for state in self.states])
        assert np.isclose(np.sum(initial), np.array(1.0)), "Initial probabilities does not sum to 1"
        return initial

    def __create_trans(self):
        """
        Transform transaction probabilities into format suitable for the Decoder class
        :return: numpy array of initial probabilities (NxN float)
        """
        trans = []
        for from_state in self.states:
            trans.append([self.trans[from_state].get(to_state, 0) for to_state in self.states])
        trans_np = np.array(trans)
        trans_sums = np.sum(trans_np, axis=1)
        assert np.isclose(trans_sums, np.ones(len(trans_sums))).all(), \
            "Transaction probabilities does not sum to 1.\n%s" % trans_sums
        return trans_np

    def __create_emit(self):
        """
        Transform emission probabilities into format suitable for the Decoder class
        :return: numpy array of initial probabilities (Nx5 float)
        """
        emit = []
        for state in self.states:
            emission = [state.emission[nucleotide] for nucleotide in NUCLEOTIDES]
            emit.append(emission)

        emit_np = np.array(emit)
        for s, _ in enumerate(self.states):
            emit_state = emit_np[s, :]
            if type(emit_state[0]) == dict:
                for q, _ in enumerate(QUALITY_CODES):
                    assert np.isclose(np.sum([emit_state[n][q] for n, _ in enumerate(NUCLEOTIDES)]), np.array(1.0)), \
                        "Emissions does not sum to 1, %s, %s" % (state, emission)
            else:
                assert np.isclose(np.sum(emission), np.array(1.0)), "Emissions does not sum to 1, %s, %s" % (state, emission)
        return emit_np

    def __transform_into_decoder_format(self):
        """
        Transform HMM into format suitable for the Viterbi Decoder
        """
        self.initial_np = self.__create_initial()
        self.trans_np = self.__create_trans()
        self.emit_np = self.__create_emit()

    def __prepare_hmm(self):
        """
        Create states and edges of the HMM based on motif code defined by user
        """
        self.__add_state(annotation.state.BackgroundState(self.background_frequency, len(self.states)))
        for module_id, (sequence, repeats) in enumerate(self.motif):
            self.__add_sequence(sequence, repeats, module_id)
        self.__add_state(annotation.state.BackgroundState(self.background_frequency, len(self.states)))
        self.__connect_states()
        self.__calculate_initial()
        self.__transform_into_decoder_format()
        self.decoder = annotation.Decoder.Decoder(self.initial_np, self.trans_np, self.emit_np)

    def annotate(self, read):
        """
        Annotate nucleotide sequence with states of the HMM model
        :param read: Read - read to be annotated
        :type read: Read
        :return: Highest probability of generating sequence with HMM states,
                 annotation for the input sequence in a string representation
        """
        seq = list(map(code_base, read.sequence))
        qual = list(map(quality_index, read.quality))
        probability, predicted_states = self.decoder.decode_log(seq, qual)
        state_seq = [self.states[i] for i in predicted_states]
        skips = [self.deleted_states.get((i, j), []) for i, j in
                 zip(predicted_states, predicted_states[1:])]
        annotation = Annotation(read, state_seq, skips, probability, self.motif)
        return annotation
