CODE_GAP = '-'
CODE_BACKGROUND = '-'


class State:
    """
    Node in the HMM architecture, with emission probability
    """

    def __init__(self, emission, module_id, state_id):
        """
        :param emission: Multinomial distribution of nucleotide generation, one float value for each nucleotide
        :param module_id: Id of the module which contains the state
        :param state_id: Id of the state in the ordered state sequence
        """
        self.emission = emission
        self.module_id = module_id
        self.state_id = state_id

    def is_background(self):
        """
        Specify, if state is outside primers and motifs and emits random nucleotides without preference
        :return: True, if state is background type
        """
        return False

    def is_insert(self):
        """
        Specify, if state is inside primers or motifs and emits random nucleotides without preference
        to simulate insertions between key states
        :return: True, if state is insert type
        """
        return False

    def is_motif(self):
        """
        Specify, if state is in motif sequence with preferred nucleotide in emission probabilities
        :return: True, if state is motif type
        """
        return False

    def is_sequence(self):
        """
        Specify, if state is in primer sequence with preferred nucleotide in emission probabilities
        :return: True, if state is sequence type
        """
        return False

    def is_key(self):
        """
        Specify, if state has specified nucleotide that is preferred in emission probabilities
        :return: True, if state is key type
        """
        return False

    def is_last_in_module(self):
        """
        Is this State last in the module?
        :return: True, if it is last - only for MotifStates relevant
        """
        return True

    def is_first_in_module(self):
        """
        Is this State last in the module?
        :return: True, if it is last - only for MotifStates relevant
        """
        return True


class KeyState(State):
    """
    State with specified nucleotide that is preferred in emission probabilities
    """

    def __init__(self, nucleotide, emission, module_id, state_id, first, last):
        State.__init__(self, emission, module_id, state_id)
        self.nucleotide = nucleotide
        self.first = first
        self.last = last

    def __str__(self):
        return self.nucleotide.lower()

    def is_sequence(self):
        return True

    def is_key(self):
        return True

    def is_first_in_module(self):
        """
        Is this State last in the module?
        :return: True, if it is last - only for MotifStates relevant
        """
        return self.first

    def is_last_in_module(self):
        """
        Is this State last in the module?
        :return: True, if it is last - only for MotifStates relevant
        """
        return self.last


class MotifState(KeyState):
    """
    State in motif sequence with preferred nucleotide in emission probabilities
    """

    def __init__(self, nucleotide, emission, module_id, state_id, first, last):
        KeyState.__init__(self, nucleotide, emission, module_id, state_id, first, last)
        self.module_id = module_id

    def __str__(self):
        return self.nucleotide.upper()

    def is_motif(self):
        return True


class SequenceState(KeyState):
    """
    State in primer sequence with preferred nucleotide in emission probabilities
    """

    def __init__(self, nucleotide, emission, module_id, state_id, first, last):
        KeyState.__init__(self, nucleotide, emission, module_id, state_id, True, True)
        self.module_id = module_id

    def __str__(self):
        return self.nucleotide.lower()

    def is_sequence(self):
        return True


class BackgroundState(State):
    """
    State outside primers and motifs, that emits random nucleotides without preference
    """

    def __init__(self, emission, state_id):
        State.__init__(self, emission, CODE_BACKGROUND, state_id)

    def __str__(self):
        return CODE_GAP

    def is_background(self):
        return True


class InsertState(State):
    """
    State inside primers or motifs, that emits random nucleotides without preference
    to simulate insertions between key states
    """

    def __init__(self, emission, module_id, state_id):
        State.__init__(self, emission, module_id, state_id)

    def __str__(self):
        return CODE_GAP

    def is_insert(self):
        return True
