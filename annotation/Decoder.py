import numpy as np
import Annotator


class Decoder(object):
    """
    Decoder object for Viterbi algorithm for annotation.
    """

    def __init__(self, initial_prob, trans_prob, obs_prob_dict):
        """
        Initializes the Decoder object with all matrices:
        :param initial_prob: 1D ndarray - initial probability of all states (length=N)
        :param trans_prob: 2D ndarray - transition probabilities for all state pairs (length=NxN)
        :param obs_prob_dict: 2D ndarray - observation probability dicts for all combinations of quality and nucleotide (length=QUALxNUCL)
        """
        self.num_states = initial_prob.shape[0]  # == N
        # we expect some probabilities to be 0 (ignore warnings)
        old_errhandling = np.seterr(divide="ignore")
        self.initial_prob = np.log(initial_prob)
        self.trans_prob = np.log(trans_prob).T
        np.seterr(**old_errhandling)

        assert self.initial_prob.shape == (self.num_states,)
        assert self.trans_prob.shape == (self.num_states, self.num_states)
        assert obs_prob_dict.shape[0] == self.num_states

        def obs_fnct(obs, asc, obs_prob):
            """
            Return the correct dictionary for observed nucleotide number and quality number.
            :param obs: int - nucleotide number
            :param asc: int - quality number
            :param obs_prob: 2D ndarray - observation probability dicts for all combinations of quality and nucleotide (length=QUALxNUCL)
            :return: ndarray - emission probabilities vector (length=N)
            """
            return np.array([item[0][asc] if isinstance(item[0], dict) else item[0] for item in obs_prob[:, obs, None]])

        # create the emission probability matrix (length=NUCLxQUALxN)
        self.obs_prob = np.zeros((obs_prob_dict.shape[1], Annotator.N_QUALITY_CODES, self.num_states))  # nucleotides * quality * num_states
        for i in range(obs_prob_dict.shape[1]):
            for j in Annotator.QUALITY_NUMS:
                self.obs_prob[i, j, :] = np.log(obs_fnct(i, j, obs_prob_dict))

    def decode_log(self, obs, quality):
        """
        Decode the optimal path using Viterbi algorithm.
        :param obs: list(int) - observed nucleotide numbers
        :param quality: list(int) - observed quality numbers
        :return: float, list(int) - likelihood of the best path, states on the best path
        """
        # initialization
        trellis = np.zeros((len(obs), self.num_states))
        backpt = np.ones((len(obs), self.num_states), 'int32') * -1
        trellis[0, :] = self.initial_prob + self.obs_prob[obs[0], quality[0]]

        # dynamic programming:
        for t in range(1, len(obs)):
            temp_mat = np.tile(trellis[t - 1, :], (self.num_states, 1)) + self.trans_prob
            backpt[t, :] = np.argmax(temp_mat, axis=1)
            trellis[t, :] = temp_mat[np.arange(self.num_states), backpt[t, :]] + self.obs_prob[obs[t], quality[t]]  # first element is the same as np.max(temp_mat, axis=1)

        # termination
        best_states = [np.argmax(trellis[-1, :])]
        for i in range(len(obs) - 1, 0, -1):
            best_states.append(backpt[i, best_states[-1]])

        # return likelihood and best states
        return np.exp(np.max(trellis[-1, :])), best_states[::-1]
