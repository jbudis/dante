from __future__ import print_function

import sys
import math
import copy
import functools
from scipy.stats import binom
import numpy as np
import matplotlib.pyplot as plt
import collections as col

import scipy.optimize

nan = float("nan")
minf = float("-inf")


def is_smth_near(profile, v1, rng=2):
    """
    Return index of something near if something exists near the target position. Else return None.
    :param profile: list(int) - profile
    :param v1: int - target place
    :param rng: int - how far to search from target
    :return: int/None - either a place where we have found something or None
    """
    # fill the search space
    to_search = [v1]
    for i in range(1, rng + 1):
        to_search.extend([v1 - i, v1 + i])
    # search the search space
    for i in to_search:
        if profile[i] != 0:
            return i

    return None


def good_for_sampling(sample, v1, v2, profiles, single=False, extract_all=False):
    """
    Tests if the sample is good for training purposes.
    :param sample: str - sample name
    :param v1: int - allele 1
    :param v2: int - allele 2
    :param profiles: pd.DataFrame - profile (counts of STRs)
    :param single: bool - extract single allele samples
    :param extract_all: bool - extract all samples (found in profiles)
    :return: bool - if the samples are "good"
    """
    if sample not in profiles.index:
        return False
    if v1 not in profiles.columns or v2 not in profiles.columns:
        return False
    if is_smth_near(profiles.loc[sample], v1) is None or is_smth_near(profiles.loc[sample], v2) is None:
        return False
    if v1 == 0 and v2 == 0:
        return False
    if extract_all:
        return True
    elif single:
        if v2 == 0:
            v2 = v1
        return v1 == v2
    else:
        ret = max(v2, v1) - min(v2, v1) > 1
        # print("Good for sampling:", sample, v1, v2, ret, profiles.loc[sample])
        return ret


def split_profile(sample, v1, v2, profiles, verbose=False):
    """
    Split profiles into 2 each for corresponding allele.
    :param sample: str - sample name
    :param v1: int - allele 1
    :param v2: int - allele 2
    :param profiles: pd.DataFrame - profiles
    :param verbose: bool - be verbose?
    :return: tuple(ndarray, ndarray) - arrays for both alleles
    """

    def fill_array(src, loc):
        """
        Copy non-empty sub-array from src starting at loc.
        :param src: ndarray - source array
        :param loc: int - starting position
        :return: ndarray - array with slice of src array around loc position
        """
        dst = np.zeros_like(src)

        dst[loc] = src[loc]
        i = loc + 1
        while i < len(src) and src[i] > 0:
            dst[i] = src[i]
            i += 1
        i = loc - 1
        while i >= 0 and src[i] > 0:
            dst[i] = src[i]
            i -= 1
        return dst

    # fill partial arrays:
    test_array = np.array(profiles.loc[sample])
    v1 = is_smth_near(profiles.loc[sample], v1)
    v2 = 0 if v2 == 0 else is_smth_near(profiles.loc[sample], v2)
    if v1 is None:
        print("ERROR: this should not happen")
        exit(-1)

    left = fill_array(test_array, v1)
    if v2 == v1 or v2 == 0 or v2 is None:
        right = []
    else:
        right = fill_array(test_array, v2)

    # check if ok:
    if sum(test_array) - sum(left) - sum(right) != 0:
        if verbose:
            print("Sample %20s: (%2d, %2d) - %5d outliers detected..." % (sample, v1, v2, sum(test_array) - sum(left) - sum(right)), file=sys.stderr,
                  end=" " if sum(test_array) - sum(left) - sum(right) < 0 else "\n")
    if sum(test_array) - sum(left) - sum(right) < 0:
        right[:(v1 + v2) // 2] = 0
        left[(v1 + v2) // 2:] = 0
        if verbose:
            print("Recounting...\nSample %20s: (%2d, %2d) - %5d outliers detected..." % (sample, v1, v2, sum(test_array) - sum(left) - sum(right)), file=sys.stderr)

    # return both arrays
    return left, right


def write_params(file_desc, model_params, read_drop_params, read_drop_params_rel, fit_function, print_all=False):
    """
    Write trained parameters to output.
    :param file_desc: file descriptor - where to write to
    :param model_params: 4-tuple(floats) - parameters of fit_function
    :param read_drop_params: 2-tuple(floats) - parameters for linear read count drop model absolute
    :param read_drop_params_rel: 2-tuple(floats) - parameters for linear read count drop model relative
    :param fit_function: function - error function for model distributions
    :param print_all: bool - whether to print all parameters in the end in computer-friendly output
    :return: None
    """
    strings = {"linear": "%f + %f * x" % (model_params[0], model_params[1]), "const": "%f" % (model_params[0]), "n2": "%f + %f * x + %f * x * x" % (model_params[0], model_params[1], model_params[2]),
               "exp": "%f + %f * exp(%f * x)" % (model_params[0], model_params[1], model_params[2])}
    print("Model parameters:", file=file_desc)
    print("\tfunction: %s" % fit_function, file=file_desc)
    print("\tdeletes: %s" % strings[fit_function], file=file_desc)
    print("\tinserts: linear with parameter %f" % model_params[3], file=file_desc)

    print("Read loss parameters (absolute):", file=file_desc)
    print("\tfunction: linear", file=file_desc)
    print("\tread count drop: %f - %f * x" % (read_drop_params[0], -read_drop_params[1]), file=file_desc)

    print("Read loss parameters (relative):", file=file_desc)
    print("\tfunction: linear", file=file_desc)
    print("\tread count drop: %f - %f * x" % (read_drop_params_rel[0], -read_drop_params_rel[1]), file=file_desc)

    if print_all:
        all_params = np.hstack((model_params, read_drop_params[:2], read_drop_params_rel[:2]))
        print("All parameters:", file=file_desc)
        print("\t%g %g %g %g %g %g %g %g" % tuple(all_params), file=file_desc)


def extract_alleles(sample, true_values, profiles, verbose=False):
    """
    Splits sample into 1 or 2 subsamples for training
    :param sample: int - sample number
    :param true_values: pd.DataFrame - true values of alleles
    :param profiles: ndarray - profile (counts of STRs)
    :param verbose: bool - be verbose?
    :return: dict{int: ndarray} - true_values to corresponding sub-profiles
    """
    v1, v2 = true_values.get_value(index=sample, col=0), true_values.get_value(index=sample, col=1)

    # get subprofiles
    left, right = split_profile(sample, v1, v2, profiles, verbose)

    # construct dict
    d = {v1: left}
    if len(right) != 0:
        d[v2] = right

    return d


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


def norm_profiles(all_profiles):
    """
    Normalize profiles dictionary
    :param all_profiles: list(dict{int:ndarray}) - list of all profiles
    :return: list(dict{int:ndarray}) - normalized all_profiles
    """
    res = copy.deepcopy(all_profiles)
    for td in res:
        for k, v in td.items():
            if sum(v) == 0:
                print("divided by 0.", sum(v), v, k, td[k])
            else:
                td[k] = v / float(sum(v))

    return res


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


def relative_read_drop_train_linear((a1, a2), p1=0.0, p2=1.0, p3=1.0):
    """
    Ratio of linear functions for relative read drop quantification.
    :param a: (list(int), list(int)) - allele numbers
    :param p1: float - constant parameter
    :param p2: float - linear parameter
    :return: float - ratio of the two linear functions
    """
    return np.array([linear_rate(aa1, p1, p2) / linear_rate(aa2, p1, p2) for aa1, aa2 in zip(a1, a2)])


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
    :param rng: ndarray - range of distribution
    :param model_params: 4-tuple - parameters for inserts and deletes
    :param n: int - target allele number
    :param rate_func: function - rate function for deletes
    :return: ndarray - combined distribution
    """
    p1, p2, p3, q = model_params
    deletes = binom.pmf(rng, n, clip(1 - rate_func(n, p1, p2, p3), 0.0, 1.0))
    inserts = binom.pmf(rng, n, q)
    return combine_distribs(deletes, inserts)


def model_template(rng, model_params, rate_func=linear_rate):
    """
    Partial function for model creation.
    :param rng: ndarray - range of distribution
    :param model_params: 4-tuple - parameters for inserts and deletes
    :param rate_func: function - rate function for deletes
    :return: partial function with only 1 parameter - n - target allele number
    """
    return functools.partial(model_full, rng, model_params, rate_func=rate_func)


def likelihood_multinomial(model_d, true_nums):
    """
    Calculates likelihood P(model_d | true_nums).  (sum(t(x))!/ (sum(t(x)))! )*prod(m(x)^t)
    :param model_d: ndarray - model distribution
    :param true_nums: ndarray - true distribution
    :return: float - log-likelihood of generation of the true distribution out of the model distribution according to multinomial
    """
    res = likelihood_simple(model_d, true_nums)
    res += math.log(math.factorial(np.sum(true_nums)))
    res -= np.sum(map(lambda x: math.log(math.factorial(x)), true_nums))
    return minf if np.isnan(res) else res


def likelihood_simple(model_d, true_d):
    """
    Calculates likelihood P(model_d | true_nums)
    :param model_d: ndarray - model distribution
    :param true_d: ndarray - true distribution
    :return: float - log-likelihood of generation of the true distribution out of the model distribution
    """
    res = sum(map(lambda (m, t): 0.0 if t == 0.0 else (minf if m == 0.0 else t * math.log(m)), zip(model_d, true_d)))
    return minf if np.isnan(res) else res


def likelihood(model_tmp, samples):
    """
    Calculated likelihood for all samples
    :param model_tmp: partial function - model template with current parameters
    :param samples: list(dict{int:ndarray}) - all profiles
    :return: float - log likelihood of current model summed on all samples
    """
    loglike = 0.0
    for true_dict in samples:
        for k, v in true_dict.items():
            md = model_tmp(k)
            loglike += likelihood_multinomial(md, v) / np.sum(v)
            # loglike += likelihood_one(md, v) / sum(v)

    return loglike


def comparison(params, samples, rate_func):
    """
    Get minus log likelihood of model specified with params on all samples
    :param params: 4-tuple - model parameters
    :param samples: list(dict{int:ndarray}) - all profiles
    :param rate_func: function - rate function for deletes
    :return: float - -log likelihood
    """
    x = np.arange(len(samples[0].values()[0]))
    model_tmp = model_template(x, params, rate_func)
    return -likelihood(model_tmp, samples)


def display(model_tmp, samples, filename, min_allele=4):
    """
    Generates a file with model vs true distribution comparison
    :param model_tmp: partial function - model template with current parameters
    :param samples: list(dict{int:ndarray}) - all profiles
    :param filename: str - filename prefix to write model comparison to
    :param min_allele: int - lowest allele to display
    :return: None
    """
    data = []
    for true_dict in samples:
        for k, v in true_dict.items():
            data.append((k, v))

    data = sorted(data, key=lambda xx: xx[0])
    data = filter(lambda xx: xx[0] >= min_allele, data)

    x = np.arange(len(data[0][1]))

    for i, (k, v) in enumerate(data):
        md = model_tmp(k)
        plt.figure()
        plt.plot(x, md * sum(v), "r-")
        plt.plot(x, v, "b-")
        plt.title(k)

        last_dot = filename.rfind('.')
        filename_curr = "%s%03d%s" % (filename[:last_dot], i, filename[last_dot:])

        plt.savefig(filename_curr)
        plt.close()


def count_occurrences(samples, len_repeating=1):
    """
    Count read counts and relative read counts.
    :param samples: list(dict{int:ndarray}) - all profiles
    :param len_repeating: int - length of the STR
    :return: defaultdict(list), defaultdict(list) - read counts and relative read counts
    """
    numbers = col.defaultdict(list)
    numbers_prop = col.defaultdict(list)

    for td in samples:

        num_v1 = 0
        allele_v1 = 0

        for k, v in td.items():
            numbers[k * len_repeating].append(sum(v))

            # we found second allele
            if allele_v1 != 0 and k != 0 and k != allele_v1:
                allele_v2 = k
                num_v2 = sum(v)

                if allele_v2 != 0:
                    if allele_v2 < allele_v1:
                        allele_v1, allele_v2 = allele_v2, allele_v1
                        num_v1, num_v2 = num_v2, num_v1

                    numbers_prop[(allele_v1 * len_repeating, allele_v2 * len_repeating)].append(num_v1 / float(num_v2))

            # first allele
            if allele_v1 == 0:
                allele_v1 = k
                num_v1 = sum(v)

    return numbers, numbers_prop


def plot_occurrences(x, y, filename, read_drop_params):
    """
    Plot read counts into file.
    :param x: list - x coordinate to plot
    :param y: list - y coordinate to plot
    :param filename: str - filename where to plot to
    :param read_drop_params: 2-tuple(floats) - parameters for linear read count drop model
    :return: None
    """
    plt.figure()
    plt.plot(x, y, "r.")
    v_lin_rate = np.vectorize(linear_rate)
    plt.plot(x, v_lin_rate(x, *read_drop_params), "g-")
    plt.savefig(filename)
    plt.close()


def flatten(numbers):
    """
    Flatten the read count dictionary into arrays for curve fitting.
    :param numbers: defaultdict(list) - read counts
    :return: list, list - list of allele numbers and list of corresponding read counts
    """
    x = []
    y = []
    for k, v in numbers.items():
        y.extend(v)
        x.extend([k] * len(v))

    return x, y


def train_read_drop_abs(nums):
    """
    Train the read drop parameters of a linear drop.
    :param nums: defaultdict(list) - read counts
    :return: list(2) - parameters of linear read drop function
    """
    params_read_drop = None
    x, y = flatten(nums)
    v_lin_err = np.vectorize(linear_rate)
    try:
        params_read_drop, _ = scipy.optimize.curve_fit(v_lin_err, x, y, (0., 1.), maxfev=20000)
    except scipy.optimize.OptimizeWarning:
        pass

    return params_read_drop


def train_read_drop_rel(nums):
    """
    Train the read drop parameters of a linear drop.
    :param nums: defaultdict(list) - read counts
    :return: list(2) - parameters of linear read drop function
    """
    params_read_drop = None
    x, y = flatten(nums)

    # p1 = 1.0
    # p2 = (1-c)/(c*a2-a1)

    p2 = [(1 - c) / (c * a2 - a1) for c, (a1, a2) in zip(y, x)]

    return 1.0, np.mean(p2)


def train_read_drop_rel2(nums):
    """
    Train the read drop parameters of a linear drop.
    :param nums: defaultdict(list) - read counts
    :return: list(2) - parameters of linear read drop function
    """
    params_read_drop = None
    x, y = flatten(nums)
    x = map(lambda a: a[0] / float(a[1]), x)
    print(x, y)
    v_lin_err = np.vectorize(linear_rate)
    try:
        params_read_drop, _ = scipy.optimize.curve_fit(v_lin_err, x, y, (0.0, 1.0), maxfev=20000)
    except scipy.optimize.OptimizeWarning:
        pass

    return params_read_drop


def train_read_drop_rel3(nums):
    """
    Train the read drop parameters of a linear drop.
    :param nums: defaultdict(list) - read counts
    :return: list(2) - parameters of linear read drop function
    """
    params_read_drop = None
    x, y = flatten(nums)
    x = map(list, zip(*x))

    # v_rdl_err = np.vectorize(relative_read_drop_train_linear)
    try:
        params_read_drop, _ = scipy.optimize.curve_fit(relative_read_drop_train_linear, x, y, (0.1, 0.9), maxfev=20000)
    except scipy.optimize.OptimizeWarning:
        pass

    return params_read_drop
