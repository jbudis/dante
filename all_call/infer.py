import all_call.train
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from copy import copy

# range of display for pcolor
min_str = 4

# NaN, -inf
nan = float("nan")
minf = float("-inf")


def print_pcolor(lh_array, display_file, name, lognorm=True):
    """
    Get maximum likelihood option and alternatively print it to image file.
    :param lh_array: 2D ndarray - likelihood array
    :param display_file: str - filename for pcolor image output
    :param name: str - name to use in title
    :param lognorm: bool - use loglog scale in displaying likelihood array
    :return: tuple(int, int) - option with highest likelihood
    """
    # get minimal and maximal likelihood
    ind_good = (lh_array < 0.0) & (lh_array > -1e10)
    lh_array[~ind_good] = minf
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
        bg_size = max(2, len(lh_view) // 8)
        if len(lh_view) - min_str <= 4:
            bg_size = 1
        lh_view[-bg_size:, min_str:min_str + bg_size] = lh_view[0, 0]

        # plotting
        plt.title("%s likelihood of each option for %s" % ("Loglog" if lognorm else "Log", name))
        plt.xlabel('2nd allele')
        plt.ylabel('1st allele')
        start_ticks = 5
        step_ticks = 5
        plt.xticks(np.array(range(start_ticks - min_str, max_str - min_str + 1, step_ticks)) + 0.5, range(start_ticks, max_str + 1, step_ticks))
        plt.yticks(np.array(range(start_ticks - min_str, max_str - min_str + 1, step_ticks)) + 0.5, range(start_ticks, max_str + 1, step_ticks))
        palette = copy(plt.cm.jet)
        palette.set_under('gray', 1.0)
        plt.pcolor(lh_view[min_str:, min_str:], cmap=palette, vmin=z_min, vmax=z_max)
        plt.colorbar()

        # background:
        plt.text(float(bg_size) / 2.0, max_str - min_str - float(bg_size) / 2.0, 'BG', size=20, horizontalalignment='center',
                 verticalalignment='center', path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")])

        # save
        plt.savefig(display_file + '.pdf')
        plt.savefig(display_file + '.png')
        plt.close()

    # output best option
    best = np.array(np.unravel_index(np.argmax(lh_array), lh_array.shape))

    return sorted([best[0], best[1]])


def combine_models(i, j, model_tmp, read_drop_params, minimal_prob, minimal_weight, len_repeating=1):
    """
    Combine 2 single-allele models/distributions into one according to read count drop specified in read_drop_params
    :param i: int, allele 1 number
    :param j: int, allele 2 number
    :param model_tmp: partial function - model template
    :param read_drop_params: 2-tuple(floats) - parameters for linear read count drop model
    :param minimal_prob: float - base chance to generate any STR
    :param minimal_weight: float - minimal weight of the model for second allele
    :param len_repeating: int - length of the STR
    :return: ndarray - combined model for both alleles
    """
    # get models
    if i == 0:
        i = j
    if j == 0:
        j = i
    model_i = model_tmp(i)
    if i == j:
        model_ij = model_i
    else:
        model_j = model_tmp(j)

        # get weights
        weight_i = all_call.train.linear_rate(i * len_repeating, *read_drop_params)
        weight_j = all_call.train.linear_rate(j * len_repeating, *read_drop_params)

        sum_w = weight_i + weight_j
        weight_i /= sum_w
        weight_j /= sum_w

        # apply some minimal weight, otherwise it will negate higher allele models
        if weight_j < minimal_weight:
            weight_j = minimal_weight
            weight_i = 1.0 - weight_j
        if weight_i < minimal_weight:
            weight_i = minimal_weight
            weight_j = 1.0 - weight_i

        # get combined model
        model_ij = weight_i * model_i + weight_j * model_j

    # assign minimal probability to all (first assign min_prob to all, such as minimal_prob = min_prob/(1-n*min_prob):
    min_prob = minimal_prob / (1.0 - len(model_ij) * minimal_prob)
    model_ij += min_prob

    return model_ij / np.sum(model_ij)


def generate_likelihoods(profile, model_params, fit_function, read_drop_params, minimal_prob, minimal_weight, background_prior=None, len_repeating=1):
    """
    Generate likelihoods of all options.
    :param profile: ndarray - profile (counts of STRs)
    :param model_params: 4-tuple(floats) - parameters of fit_function
    :param fit_function: function - error function for model distributions
    :param read_drop_params: 2-tuple(floats) - parameters for linear read count drop model
    :param minimal_prob: float - base chance to generate any STR
    :param minimal_weight: float - minimal weight of the model for second allele
    :param background_prior: float - how much more likely is to have background than one concrete allele
    :param len_repeating: int - length of the STR
    :return: 2D ndarray - likelihoods of all options
    """
    # create options for all up to +2 more than last
    max_str = len(profile) + 2
    lh_array = np.zeros((max_str, max_str))
    num_models = 0

    # crete model template
    model_tmp = all_call.train.model_template(np.arange(max_str), model_params, fit_function)

    # compute likelihoods for all 2 allele options
    for i in range(min_str, max_str):
        for j in range(i, max_str):
            model_ij = combine_models(i, j, model_tmp, read_drop_params, minimal_prob, minimal_weight, len_repeating=len_repeating)
            lh_array[i, j] = all_call.train.likelihood_multinomial(model_ij, profile)
            num_models += 1

    # compute likelihood with no alleles (background)
    if background_prior is None:
        background_prior = num_models
    lh_array[0, 0] = min(all_call.train.likelihood_multinomial(np.ones_like(model_ij, dtype=float) / len(model_ij), profile) + np.log(background_prior), -0.000001)

    return lh_array


def write_output(file_desc, predicted, conf, name):
    """
    Write result of one prediction.
    :param file_desc: file descriptor - where to write to
    :param predicted: tuple(int, int) - predicted alleles
    :param conf: tuple(float, float, float) - confidence of prediction (whole, 1st allele, 2nd allele)
    :param name: str/int - name/number of the sample
    :return: None
    """
    print("Predicted alleles for %s: (confidence = %5.1f%%)" % (str(name), conf[0] * 100.0), file=file_desc)
    print("\t%3d (confidence = %5.1f%%)" % (predicted[0], conf[1] * 100.0), file=file_desc)
    print("\t%3d (confidence = %5.1f%%)\n" % (predicted[1], conf[2] * 100.0), file=file_desc)
