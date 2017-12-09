from __future__ import print_function

import math
import functools
from scipy.stats import binom
import numpy as np
import itertools
# from annotation import Annotation
import sys

import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from copy import copy


def combine_distribs(deletes, inserts):
    """
    Combine insert and delete models/distributions
    :param deletes: ndarray - delete distribution
    :param inserts: ndarray - insert distribution
    :return: ndarray - combined array of the same length
    """

    # how much to fill?
    to_fill = sum(deletes == 0.0) + 1
    while to_fill < len(inserts) and inserts[to_fill] > 0.0001:
        to_fill += 1

    # create the end array
    len_del = len(deletes)
    end_distr = np.zeros_like(deletes, dtype=float)

    # fill it!
    for i, a in enumerate(inserts[:to_fill]):
        # print i,a,(deletes*a)[:len_del-i]
        end_distr[i:] += (deletes * a)[:len_del - i]

    # print("end_distr", end_distr[:3], deletes[:3], inserts[:3])
    return end_distr


def const_rate(n, p1=0.0, p2=1.0, p3=1.0):
    """
    Constant rate function.
    :param n: int - allele number (unused)
    :param p1: float - constant parameter
    :param p2: float - linear parameter (unused)
    :param p3: float - additional parameter (unused)
    :return: float - p1
    """
    return p1


def linear_rate(n, p1=0.0, p2=1.0, p3=1.0):
    """
    Linear rate function.
    :param n: int - allele number
    :param p1: float - constant parameter
    :param p2: float - linear parameter
    :param p3: float - additional parameter (unused)
    :return: float - p1 + p2 * n
    """
    return p1 + p2 * n


def n2_rate(n, p1=0.0, p2=1.0, p3=1.0):
    """
    Quadratic rate function.
    :param n: int - allele number
    :param p1: float - constant parameter
    :param p2: float - linear parameter
    :param p3: float - quadratic parameter
    :return: float - p1 + p2 * n + p3 * n * n
    """
    return p1 + p2 * n + p3 * n * n


def exp_rate(n, p1=0.0, p2=1.0, p3=1.0):
    """
    Exponential rate function.
    :param n: int - allele number
    :param p1: float - constant parameter
    :param p2: float - linear parameter
    :param p3: float - exponential parameter
    :return: float - p1 + p2 * e^(p3 * n)
    """
    return p1 + p2 * math.exp(p3 * n)


def clip(value, minimal, maximal):
    """
    Clips value to range <minimal, maximal>
    :param value: ? - value
    :param minimal: ? - minimal value
    :param maximal: ? - maximal value
    :return: ? - clipped value
    """
    return min(max(minimal, value), maximal)


def model_full(rng, model_params, n, rate_func=linear_rate):
    """
    Create binomial model for both deletes and inserts of STRs
    :param rng: int - max_range of distribution
    :param model_params: 4-tuple - parameters for inserts and deletes
    :param n: int - target allele number
    :param rate_func: function - rate function for deletes
    :return: ndarray - combined distribution
    """
    p1, p2, p3, q = model_params
    deletes = binom.pmf(np.arange(rng), n, clip(1 - rate_func(n, p1, p2, p3), 0.0, 1.0))
    inserts = binom.pmf(np.arange(rng), n, q)
    return combine_distribs(deletes, inserts)


def model_template(rng, model_params, rate_func=linear_rate):
    """
    Partial function for model creation.
    :param rng: int - max_range of distribution
    :param model_params: 4-tuple - parameters for inserts and deletes
    :param rate_func: function - rate function for deletes
    :return: partial function with only 1 parameter - n - target allele number
    """
    return functools.partial(model_full, rng, model_params, rate_func=rate_func)


class Inference:
    """ Class for inference of alleles. """

    MIN_REPETITIONS = 1

    # default parameters for inference
    DEFAULT_MODEL_PARAMS = (-0.0107736, 0.00244419, 0.0, 0.00440608)
    DEFAULT_FIT_FUNCTION = "linear"

    def __init__(self, read_distribution, params_file, str_rep=3, minl_primer1=5, minl_primer2=5, minl_str=5, p_bckg_closed=None, p_bckg_open=None, p_expanded=None):
        """
        Initialization of the Inference class + setup of all models and their probabilities.
        :param read_distribution: ndarray(int) - read distribution
        :param params_file: str - filename of parameters
        :param str_rep: int - length of the STR
        :param minl_primer1: int - minimal length of the left primer
        :param minl_primer2: int - minimal length of the right primer
        :param minl_str: int - minimal length of the STR
        :param p_bckg_closed: float - probability of the background model for closed observation
        :param p_bckg_open: float - probability of the background model for open observation
        :param p_expanded: float - probability of the expanded model (if None it is equal to other models)
        """
        # assign variables
        self.str_rep = str_rep
        self.minl_primer1 = minl_primer1
        self.minl_primer2 = minl_primer2
        self.minl_str = minl_str
        self.read_distribution = read_distribution
        self.sum_reads_log = np.log(np.sum(read_distribution))
        self.sum_reads = np.sum(read_distribution)
        self.params_file = params_file
        self.p_expanded = p_expanded
        self.p_bckg_closed = p_bckg_closed
        self.p_bckg_open = p_bckg_open

    def construct_models(self, min_rep, max_rep, e_model):
        """
        Construct all models needed for current inference.
        :param min_rep: int - minimal allele to model
        :param max_rep: int - maximal allele to model
        :param e_model: int - model for expanded alleles
        :return: None
        """

        # extract params
        model_params, rate_func_str = self.read_params(self.params_file)
        str_to_func = {"linear": linear_rate, "const": const_rate, "exponential": exp_rate, "square": n2_rate}
        rate_func = const_rate
        if rate_func_str in str_to_func.keys():
            rate_func = str_to_func[rate_func_str]

        # save min_rep and max_rep
        self.min_rep = min_rep
        self.max_rep = max_rep  # non-inclusive
        self.max_with_e = e_model + 1  # non-inclusive

        # get models
        mt = model_template(self.max_with_e, model_params, rate_func)
        self.background_model = np.concatenate([np.zeros(self.min_rep, dtype=float), np.ones(self.max_with_e - self.min_rep, dtype=float) / float(self.max_with_e - self.min_rep)])
        self.expanded_model = mt(self.max_with_e - 1)
        self.allele_models = {i: mt(i) for i in range(min_rep, max_rep)}
        self.models = {'E': self.expanded_model, 'B': self.background_model}
        self.models.update(self.allele_models)

        # get model likelihoods
        open_to_closed = 10.0

        l_others = 1.0
        l_bckg_open = 0.01
        l_exp = 1.01

        l_bckg_model_open = 1.0

        if self.p_expanded is None:
            self.p_expanded = l_exp
        if self.p_bckg_open is None and self.p_bckg_closed is None:
            self.p_bckg_open = l_bckg_open
            self.p_bckg_closed = self.p_bckg_open / open_to_closed
        if self.p_bckg_closed is None:
            self.p_bckg_closed = self.p_bckg_open / open_to_closed
        if self.p_bckg_open is None:
            self.p_bckg_open = self.p_bckg_closed * open_to_closed

        self.model_probabilities = {'E': self.p_expanded, 'B': l_bckg_model_open}
        self.model_probabilities.update({i: l_others for i in self.allele_models.keys()})

    def read_params(self, params_file):
        """
        Reads all parameters written with write_params(print_all=True)
        :param params_file: str - filename to read parameters from, if None, load default params
        :return: 4-tuple, 2-tuple, function - parameters for model, read count drop, and error function for model distributions
        """
        if params_file is None:
            return self.DEFAULT_MODEL_PARAMS, self.DEFAULT_FIT_FUNCTION

        # read 2nd and last line of the file
        with open(params_file) as f:
            lines = f.readlines()
            fit_function = lines[1].strip().split()[1]
            split = map(float, lines[-1].strip().split())

        if len(split) < 4:
            print("ERROR: parameters were not read successfully, using defaults!", file=sys.stderr)
            return self.DEFAULT_MODEL_PARAMS, self.DEFAULT_FIT_FUNCTION

        # extract parameters from last line of file
        model_params = tuple(split[0:4])

        return model_params, fit_function

    def likelihood_rl(self, rl):
        """
        Likelihood of a read with this length.
        :param rl: int - read length
        :return: float - likelihood of a read this long
        """
        # print('rl', self.read_distribution[rl] / float(self.sum_reads))
        return self.read_distribution[rl] / float(self.sum_reads)

    @staticmethod
    def likelihood_model(model, g):
        """
        Likelihood of a generated allele al from a model of
        :param model: ndarray - model that we evaluate
        :param g: int - observed read count
        :return: float - likelihood of a read coming from this model
        """
        return model[g]

    def likelihood_intersection(self, model_i, model_j, g):
        return min(model_i[g], model_j[g])

    def likelihood_coverage(self, true_length, rl, closed=True):
        """
        Likelihood of generating a read with this length and this allele.
        :param true_length: int - true number of repetitions of an STR
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return: float - likelihood of a read being generated with this attributes
        """
        whole_inside_str = max(0, true_length * self.str_rep + self.minl_primer1 + self.minl_primer2 - rl + 1)
        # closed_overlapping = max(0, rl - self.minl_primer1 - self.minl_primer2 - true_length * self.str_rep + 1)
        open_overlapping = max(0, rl + true_length * self.str_rep - 2 * self.minl_str + 1)

        assert open_overlapping > whole_inside_str, '%d open %d whole inside %d %d %d' % (open_overlapping, whole_inside_str, true_length, rl, self.minl_str)

        return 1.0 / float(open_overlapping - whole_inside_str)

    def likelihood_read_allele(self, model, observed, rl, closed=True):
        """
        Likelihood of generation of read with observed allele count and rl.
        :param model: ndarray - model for the allele
        :param observed: int - observed allele count
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return:
        """
        if closed:
            return self.likelihood_rl(rl) * self.likelihood_model(model, observed) * self.likelihood_coverage(observed, rl, True)
        else:
            number_of_options = 0
            partial_likelihood = 0
            for true_length in itertools.chain(range(observed, self.max_rep), [self.max_with_e - 1]):
                partial_likelihood += self.likelihood_model(model, true_length) * self.likelihood_coverage(true_length, rl, False)
                number_of_options += 1

            return self.likelihood_rl(rl) * partial_likelihood / float(number_of_options)

    def likelihood_read_intersection(self, model_i, model_j, observed, rl, closed=True):
        """
        Likelihood of generation of read with observed allele count and rl.
        :param model: ndarray - model for the allele
        :param observed: int - observed allele count
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return:
        """
        if closed:
            return self.likelihood_rl(rl) * self.likelihood_intersection(model_i, model_j, observed) * self.likelihood_coverage(observed, rl, True)
        else:
            number_of_options = 0
            partial_likelihood = 0
            for true_length in itertools.chain(range(observed, self.max_rep), [self.max_with_e - 1]):
                partial_likelihood += self.likelihood_intersection(model_i, model_j, true_length) * self.likelihood_coverage(true_length, rl, False)
                number_of_options += 1

            return self.likelihood_rl(rl) * partial_likelihood / float(number_of_options)

    def likelihood_read(self, observed, rl, model_index1, model_index2, closed=True):
        """
        Compute likelihood of generation of a read from either of those models.
        :param observed: int - observed allele count
        :param rl: int - read length
        :param model_index1: char/int - model index for left allele
        :param model_index2: char/int - model index for right allele
        :param closed: bool - if the read is closed - i.e. both primers are therse
        :return: float - likelihood of this read generation
        """

        # print('testing', model_index1, model_index2)

        model_i = self.models[model_index1]
        model_j = self.models[model_index2]

        model_prob_i = self.model_probabilities[model_index1]
        model_prob_j = self.model_probabilities[model_index2]

        # TODO: tuto podla mna nemoze byt len tak +, chyba tam korelacia modelov, ale v ramci zjednodusenia asi ok
        allele1_likelihood = model_prob_i * self.likelihood_read_allele(model_i, observed, rl, closed)
        allele2_likelihood = model_prob_j * self.likelihood_read_allele(model_j, observed, rl, closed)
        p_bckg = self.p_bckg_closed if closed else self.p_bckg_open
        bckgrnd_likelihood = p_bckg * self.likelihood_read_allele(self.models['B'], observed, rl, closed)

        # alleles_intersection = min(model_prob_j, model_prob_i) * self.likelihood_read_intersection(model_i, model_j, observed, rl, closed)
        # if alleles_intersection > 0.0:
        #    print('%g %g %g %s %s %d' % (alleles_intersection, allele2_likelihood, allele1_likelihood, str(model_index1), str(model_index2), observed))

        assert not np.isnan(allele2_likelihood)
        assert not np.isnan(allele1_likelihood)
        assert not np.isnan(bckgrnd_likelihood)
        # assert alleles_intersection <= max(allele1_likelihood, allele2_likelihood), '%g %g %g %s %s %d' % (
        #    alleles_intersection, allele2_likelihood, allele1_likelihood, str(model_index1), str(model_index2), observed)

        # print('read_%s' % (str(closed)), observed, 'all1_lh', allele1_likelihood, 'all2_lh', allele2_likelihood)

        return allele1_likelihood + allele2_likelihood + bckgrnd_likelihood  # - alleles_intersection

    def infer(self, annotations, filt_annotations, index_rep, verbose=True):
        """
        Does all of the inference, computes for which 2 combination of alleles are these annotations and parameters the best.
        argmax_{G1, G2} P(G1, G2 | AL, COV, RL) ~ P(AL, COV, RL | G1, G2) * P(G1, G2) = prod_{read_i} P(al_i, cov_i, rl_i | G1, G2) * P(G1, G2) =independent G1 G2=
         = prod_{read_i} P(al_i, cov_i, rl_i | G1) * P(al_i, cov_i, rl_i | G2) * P(G1) * P(G2) {here G1, G2 is from possible alleles, background, and expanded, priors are from params}
         P(al_i, cov_i, rl_i | G1) - 2 options: 1. closed evidence (al_i = X), we know X; 2. open evidence (al_i >= X), cl_i == True if i is closed
         1.: P(al_i, cov_i, rl_i, cl_i | G1) = P(rl_i is from read distribution) * p(allele is al_i | G1) * P(read generated closed evidence | rl_i, al_i)
         2.: P(rl_i is from r.distr.) * P(allele is >= al_i | G1) * P(read generated open evidence | rl_i, al_i)
        :param annotations: iterator(reads) - closed reads (both primers set)
        :param filt_annotations: iterator(reads) - open reads (only one primer set)
        :param index_rep: int - index of a repetition
        :param verbose: bool - print more stuff?
        :return: dict(tuple(int, int):float) - directory of model indices to their likelihood
        """
        # generate closed observed and read_length arrays
        observed_annots = map(lambda x: x.module_repetitions[index_rep], annotations)
        rl_annots = map(lambda x: len(x.read.sequence), annotations)
        closed_annots = np.ones_like(observed_annots, dtype=bool)

        # generate open observed and read_length arrays
        observed_fa = map(lambda x: x.module_repetitions[index_rep], filt_annotations)
        rl_fa = map(lambda x: len(x.read.sequence), filt_annotations)
        closed_fa = np.zeros_like(observed_fa, dtype=bool)

        # join them and keep the information if they are open or closed
        observed_arr = np.concatenate([observed_annots, observed_fa]).astype(int)
        rl_arr = np.concatenate([rl_annots, rl_fa]).astype(int)
        closed_arr = np.concatenate([closed_annots, closed_fa]).astype(bool)

        # generate the boundaries:
        overhead = 3
        if len(observed_annots) == 0:
            max_rep = max(observed_fa) + overhead  # non-inclusive
            min_rep = max(self.MIN_REPETITIONS, max(observed_fa) - overhead)  # inclusive
        else:
            max_rep = max(observed_annots) + overhead + 1  # non-inclusive
            min_rep = max(self.MIN_REPETITIONS, min(observed_annots) - overhead)  # inclusive

        # expanded allele
        e_allele = max_rep
        if len(observed_fa) > 0:
            e_allele = max(max_rep, max(observed_fa) + 1)

        # generate all the models
        self.construct_models(min_rep, max_rep, e_allele)

        tested_models = []

        for model_index1 in range(min_rep, max_rep):
            for model_index2 in range(model_index1, max_rep):
                tested_models.append((model_index1, model_index2))
            tested_models.append((model_index1, 'E'))
            # tested_models.append(('B', model_index1))

        tested_models.append(('B', 'B'))
        tested_models.append(('E', 'E'))

        # go through every model and evaluate:
        evaluated_models = {}
        for m1, m2 in tested_models:

            evaluated_models[(m1, m2)] = 0

            if verbose:
                print('model', m1, m2)

            # go through every reads
            for obs, rl, closed in zip(observed_arr, rl_arr, closed_arr):
                lh = self.likelihood_read(obs, rl, m1, m2, closed=closed)
                # TODO weighted sum according to the closeness/openness of reads?
                evaluated_models[(m1, m2)] += np.log(lh)

            if verbose:
                print('model', m1, m2, 'log-likelihood', evaluated_models[(m1, m2)])

        return evaluated_models

    def print_pcolor(self, lh_dict, display_file, name, lognorm=True):
        """
        Get maximum likelihood option and alternatively print it to image file.
        :param lh_dict: dict(tuple(int, int):float) - directory of model indices to their likelihood
        :param display_file: str - filename for pcolor image output
        :param name: str - name to use in title
        :param lognorm: bool - use loglog scale in displaying likelihood array
        :return: tuple(int, int) - option with highest likelihood
        """
        # convert to a numpy array:
        lh_array = np.zeros((self.max_rep, self.max_rep + 1))
        for (k1, k2), v in lh_dict.items():
            if k1 == 'B':
                k1 = 0
            if k2 == 'B':
                k2 = 0
            if k1 == 'E':
                k1 = 0
            if k2 == 'E':
                k2 = self.max_rep
            lh_array[k1, k2] = v

        # print(lh_dict, lh_array)

        # get minimal and maximal likelihood
        ind_good = (lh_array < 0.0) & (lh_array > -1e10) & (lh_array != np.nan)
        if len(lh_array[ind_good]) == 0:
            return lh_array, (0, 0)
        lh_array[~ind_good] = np.NINF
        z_min, z_max = min(lh_array[ind_good]), max(lh_array[ind_good])

        max_str = len(lh_array)

        # generate image file if specified:
        if display_file is not None:
            plt.figure()

            if lognorm:
                lh_view = -np.log(-lh_array)
                z_min = -np.log(-z_min)
                z_max = -np.log(-z_max)
            else:
                lh_view = lh_array

            # background:
            bg_size = max(2, (len(lh_view) - self.min_rep) // 6)
            if len(lh_view) - self.min_rep <= 6:
                bg_size = 1
            lh_view[-bg_size:, self.min_rep:self.min_rep + bg_size] = lh_view[0, 0]
            # expanded
            lh_view[-bg_size:, self.min_rep + bg_size:self.min_rep + 2 * bg_size] = lh_view[0, self.max_rep]

            # plotting
            plt.title("%s likelihood of each option for %s" % ("Loglog" if lognorm else "Log", name))
            plt.xlabel('2nd allele')
            plt.ylabel('1st allele')
            start_ticks = 5
            step_ticks = 5
            plt.xticks(np.concatenate([np.array(range(start_ticks - self.min_rep, max_str - self.min_rep, step_ticks)), [max_str - self.min_rep]]) + 0.5,
                       list(range(start_ticks, max_str, step_ticks)) + ['E(>%d)' % (self.max_with_e - 2)])
            plt.yticks(np.array(range(start_ticks - self.min_rep, max_str - self.min_rep, step_ticks)) + 0.5, range(start_ticks, max_str, step_ticks))
            palette = copy(plt.cm.jet)
            palette.set_under('gray', 1.0)
            plt.pcolor(lh_view[self.min_rep:, self.min_rep:], cmap=palette, vmin=z_min, vmax=z_max)
            plt.colorbar()

            # draw dividing line:
            plt.plot([max_str - self.min_rep, max_str - self.min_rep], [0, max_str - self.min_rep], 'k', linewidth=3)

            # background:
            plt.text(float(bg_size) / 2.0, max_str - self.min_rep - float(bg_size) / 2.0, 'BG', size=20, horizontalalignment='center',
                     verticalalignment='center', path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")])
            # expanded
            plt.text(bg_size + float(bg_size) / 2.0, max_str - self.min_rep - float(bg_size) / 2.0, 'Exp', size=20, horizontalalignment='center',
                     verticalalignment='center', path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")])

            # save
            plt.savefig(display_file + '.pdf')
            plt.savefig(display_file + '.png')
            plt.close()

        # output best option
        best = sorted(np.unravel_index(np.argmax(lh_array), lh_array.shape))

        # and convert it to symbols
        if best[0] == 0 and best[1] == 0:
            best_sym = ('B', 'B')
        else:
            best_sym = map(lambda x: 'E' if x == self.max_rep or x == 0 else x, best)

        return lh_array, best, best_sym

    @staticmethod
    def get_confidence(lh_array, predicted):
        """
        Get confidence of a prediction.
        :param lh_array: 2D-ndarray - log likelihoods of the prediction
        :param predicted: tuple(int, int) - predicted alleles
        :return: tuple(float, float, float) - prediction confidence of all, first, and second allele(s)
        """

        # get confidence
        lh_corr_array = lh_array - np.max(lh_array)
        lh_sum = np.sum(np.exp(lh_corr_array))
        confidence = np.exp(lh_corr_array[predicted[0], predicted[1]]) / lh_sum
        confidence1 = np.sum(np.exp(lh_corr_array[predicted[0], :])) / lh_sum
        confidence2 = np.sum(np.exp(lh_corr_array[:, predicted[1]])) / lh_sum

        return confidence, confidence1, confidence2

    @staticmethod
    def write_output(file_desc, predicted, conf, name):
        """
        Write result of one prediction.
        :param file_desc: file descriptor - where to write to
        :param predicted: tuple(int/char, int/char) - predicted alleles
        :param conf: tuple(float, float, float) - confidence of prediction (whole, 1st allele, 2nd allele)
        :param name: str/int - name/number of the sample
        :return: None
        """

        def write_output_fd(f, predicted, conf, name):
            print("Predicted alleles for %s: (confidence = %5.1f%%)" % (str(name), conf[0] * 100.0), file=f)
            print("\t%3s (confidence = %5.1f%%)" % (str(predicted[0]), conf[1] * 100.0), file=f)
            print("\t%3s (confidence = %5.1f%%)\n" % (str(predicted[1]), conf[2] * 100.0), file=f)

        if type(file_desc) is str:
            with open(file_desc, 'w') as f:
                write_output_fd(f, predicted, conf, name)
        else:
            write_output_fd(file_desc, predicted, conf, name)

    def all_call(self, annotations, filt_annotations, index_rep, file_pcolor, file_output, name):
        """
        Run All_call - inference of likelihoods, printing of pcolor and writing output.
        :param annotations: list(Annotation) - good (blue) annotations
        :param filt_annotations: list(Annotation) - (grey) annotations with one primer
        :param index_rep: int - index of a repetition
        :param file_pcolor: str - file prefix for a pcolor image
        :param file_output: str - file for all_call output
        :param name: str - name of the sample
        :return: None
        """

        # if we do not have any good annotations, then quit
        if len(annotations) == 0 and len(filt_annotations) == 0:
            # write output
            # self.write_output(file_output, ('B', 'B'), (0.0, 0.0, 0.0), name)

            return None

        # infer likelihoods
        lh_dict = self.infer(annotations, filt_annotations, index_rep, verbose=False)

        # print pcolor image
        lh_array, predicted, predicted_sym = self.print_pcolor(lh_dict, file_pcolor, name)

        # get confidence of our prediction
        conf = self.get_confidence(lh_array, predicted)

        # write output
        self.write_output(file_output, predicted_sym, conf, name)
