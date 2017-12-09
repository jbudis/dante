from __future__ import print_function

import sys
import numpy as np
from datetime import datetime

from all_call import train, infer, input
import scipy.optimize

MIN_HITS = 10


def end_it(start_time):
    """
    Write end time.
    :param start_time: datetime object - start time o f the program
    :return: None
    """
    # print the time of the end:
    end_time = datetime.now()
    print('AllCall Stopping  : {finish:%Y-%m-%d %H:%M:%S}'.format(finish=end_time), file=sys.stderr)
    print('Total time of run : {duration}'.format(duration=end_time - start_time), file=sys.stderr)


def change_suffix(filename, suffix):
    """
    Change suffix into a new one.
    :param filename: str - filename
    :param suffix: str - string to append
    :return: str - new filename with new suffix
    """
    return filename[:filename.rfind('.')] + suffix


# print time of the start:
start_time = datetime.now()
print('AllCall = Allele Caller (for Dante)', file=sys.stderr)
print('AllCall Starting  : {start:%Y-%m-%d %H:%M:%S}'.format(start=start_time), file=sys.stderr)

args = input.load_arguments()
print("Arguments read successfully...", file=sys.stderr)

# training:
if args.train:
    print("Training mode...", file=sys.stderr)

    # read profile(s)
    if args.old_profiles:
        profiles = input.read_profiles(args.profiles)
    else:
        profiles = input.read_dante(args.profiles)

    # read_true values
    true_values = input.read_true(args.true_values)

    if profiles is not None and true_values is not None:
        print("Training data read successfully:\n"
              "  %3d samples read\n"
              "  %3d true values read" % (len(profiles.index), len(true_values)), file=sys.stderr)
    else:
        if profiles is None:
            print("ERROR: wrong type of file: %s", args.profiles)
        if true_values is None:
            print("ERROR: wrong type of file: %s", args.true_values)
        end_it(start_time)
        exit(-1)

    # filter profiles for good for training:
    test_samples = filter(lambda sample: train.good_for_sampling(sample, true_values[sample][0], true_values[sample][1], profiles, single=False), true_values)
    single_samples = filter(lambda sample: train.good_for_sampling(sample, true_values[sample][0], true_values[sample][1], profiles, single=True), true_values)
    print("Extracted:\n"
          "  %3d good training samples\n"
          "  %3d single values training samples" % (len(test_samples), len(single_samples)), file=sys.stderr)

    # samples_struct is list(dict(int:ndarray(int))) - list of all samples with dicts of alleles and their respective arrays
    all_samples = test_samples + single_samples
    samples_struct = map(lambda sample: train.extract_alleles(sample, true_values, profiles, args.verbosity_level > 1), all_samples)

    # normalize (for model training):
    samples_struct_norm = train.norm_profiles(samples_struct)

    # train the model:
    print("Training of the binomial model...", file=sys.stderr)
    start = [0.005, 0.0005, 0.01, 0.005]
    params = scipy.optimize.fmin(train.comparison, start, args=(samples_struct, input.fit_functions[args.fit_function]), xtol=0.000001, maxfun=1000000, maxiter=1000000)

    # display?
    if args.model_fig is not None:
        print("Generating model figures...", file=sys.stderr)
        rng = np.arange(len(samples_struct[0].values()[0]))
        model_tmp = train.model_template(rng, params, input.fit_functions[args.fit_function])
        train.display(model_tmp, samples_struct_norm, args.model_fig)

    # count the occurrences:
    nums, nums_prop = train.count_occurrences(samples_struct, len_repeating=args.len_repeating)

    # train read count drop:
    print("Training of linear read count drop (absolute)...", file=sys.stderr)
    params_read_drop = train.train_read_drop_abs(nums)
    if params_read_drop is None:
        params_read_drop = input.DEFAULT_READ_DROP

    print("Training of linear read count drop (relative)...", file=sys.stderr)
    params_read_drop_rel = train.train_read_drop_rel(nums_prop)
    if params_read_drop_rel is None:
        params_read_drop_rel = input.DEFAULT_READ_DROP

    # display?
    if args.model_fig is not None:
        print("Generating read count drop fit...", file=sys.stderr)
        last_dot = args.model_fig.rfind('.')
        filename_curr = "%s_reads%s" % (args.model_fig[:last_dot], args.model_fig[last_dot:])
        x, y = train.flatten(nums)
        train.plot_occurrences(x, y, filename_curr, params_read_drop)
        filename_curr = "%s_reads_relative%s" % (args.model_fig[:last_dot], args.model_fig[last_dot:])
        x5, y5 = train.flatten(nums_prop)
        train.plot_occurrences(x5, y5, filename_curr, params_read_drop_rel)

    # write parameters to the file:
    print("Writing the parameters to output (%s)..." % args.params_file if args.params_file is not None else "stdout", file=sys.stderr)
    if args.params_file is not None:
        with open(args.params_file, 'w') as f:
            train.write_params(f, params, params_read_drop, params_read_drop_rel, args.fit_function, print_all=True)
    else:
        train.write_params(sys.stdout, params, params_read_drop, params_read_drop_rel, args.fit_function, print_all=True)

# inference:
else:
    print("Inference mode...", file=sys.stderr)

    # read profile(s)
    if args.old_profiles:
        profiles = input.read_profiles(args.profiles)
    else:
        profiles = input.read_dante(args.profiles)

    # read true values
    if args.true_values is not None:
        true_values = input.read_true(args.true_values)
    else:
        true_values = None

    # print status
    if profiles is not None:
        print("Training data read successfully:\n"
              "  %3d samples read" % len(profiles.index), file=sys.stderr)
    else:
        if profiles is None:
            print("ERROR: wrong type of file or file is empty: %s", args.profiles)
        end_it(start_time)
        exit(-1)

    # read parameters:
    model_params, read_drop_params, read_drop_params_rel, args.fit_function = input.read_params(args.params_file)
    if args.verbosity_level > 1:
        train.write_params(sys.stderr, model_params, read_drop_params, read_drop_params_rel, args.fit_function)
    used_read_drop = read_drop_params_rel if args.relative else read_drop_params

    # infer the numbers with their corresponding likelihood
    correct = 0
    semi_correct = 0
    found_trues = 0

    correct_min = 0
    semi_correct_min = 0
    found_trues_min = 0

    # filter based on index_rep
    if args.index_rep is not None:
        profiles = profiles.iloc[args.index_rep - 1:args.index_rep]

    # iterate through profiles
    for i, row in enumerate(profiles.itertuples()):
        name = row[0]
        profile = np.array(row[1:])
        ind_last = np.arange(len(profile))[profile.nonzero()][-1]
        profile = profile[:ind_last + 1]

        # change name of the output
        if args.index_rep is None:
            if args.pcolor is not None:
                args.pcolor = '%s_%d' % (args.pcolor, i + 1)
            if args.output is not None:
                args.output = '%s_%d.txt' % (args.output, i + 1)

        # generate likelihood array
        lh_array = infer.generate_likelihoods(profile, model_params, input.fit_functions[args.fit_function], used_read_drop, args.minimal_prob, args.minimal_weight, len_repeating=args.len_repeating)

        # get best option
        predicted = infer.print_pcolor(lh_array, args.pcolor, name)

        # get confidence (background is at [0,0])
        lh_corr_array = lh_array - np.max(lh_array)
        lh_sum = np.sum(np.exp(lh_corr_array))
        confidence = np.exp(lh_corr_array[predicted[0], predicted[1]]) / lh_sum
        confidence1 = np.sum(np.exp(lh_corr_array[predicted[0], :])) / lh_sum
        confidence2 = np.sum(np.exp(lh_corr_array[:, predicted[1]])) / lh_sum
        conf = (confidence, confidence1, confidence2)

        # compare to true value:
        if true_values is not None:
            if name in true_values:
                tv = true_values[name]
                if tv[0] == 0:
                    tv[0] = tv[1]
                if tv[1] == 0:
                    tv[1] = tv[0]
                correct1 = tv[0] == predicted[0]
                correct2 = tv[1] == predicted[1]
                correctness = "CORRECT" if correct1 and correct2 else ("WRONG!!" if sum(profile) >= MIN_HITS else "%d hits!" % sum(profile))
                print("True values for  %15s: (%s)" % (name, correctness), file=sys.stderr)
                print("  Allele 1: %3d (%s)" % (tv[0], " OK" if correct1 else "BAD"), file=sys.stderr)
                print("  Allele 2: %3d (%s)" % (tv[1], " OK" if correct2 else "BAD"), file=sys.stderr)

                correct += correct1 and correct2
                semi_correct += correct1 or correct2

                found_trues += 1

                if sum(profile) >= MIN_HITS:
                    correct_min += correct1 and correct2
                    semi_correct_min += correct1 or correct2

                    found_trues_min += 1

        # output...
        if args.output is not None:
            with open(args.output + '.txt', 'w') as f:
                infer.write_output(f, predicted, conf, name)
        else:
            infer.write_output(sys.stdout, predicted, conf, name)

        # optionally output old profiles and true_values for training
        if args.output_profile is not None:
            with open(args.output_profile, 'a') as f:
                print('\t'.join(map(str, row)), file=f)
            with open(change_suffix(args.output_profile, '.true'), 'a') as f:
                print('%s\t%s\t%s' % (name, predicted[0], predicted[1]), file=f)

    if true_values is not None:
        print("Prediction accuracy:", file=sys.stderr)
        print("  correct samples: %3d/%3d (%5.1f%%)" % (correct, found_trues, float('Nan') if found_trues == 0 else correct / float(found_trues) * 100.0), file=sys.stderr)
        print("  correct alleles: %3d/%3d (%5.1f%%)" % (correct + semi_correct, found_trues * 2, float('Nan') if found_trues == 0 else (correct + semi_correct) / float(found_trues) / 2 * 100.0),
              file=sys.stderr)

    if true_values is not None:
        print("Prediction accuracy for profiles with more than %d hits:" % MIN_HITS, file=sys.stderr)
        print("  correct samples: %3d/%3d (%5.1f%%)" % (correct_min, found_trues_min, float('Nan') if found_trues_min == 0 else correct_min / float(found_trues_min) * 100.0), file=sys.stderr)
        print("  correct alleles: %3d/%3d (%5.1f%%)" % (
            correct_min + semi_correct_min, found_trues_min * 2, float('Nan') if found_trues_min == 0 else (correct_min + semi_correct_min) / float(found_trues_min) / 2 * 100.0),
              file=sys.stderr)

# write end time
end_it(start_time)

# training:
# python allcall.py -t --model-fig all_call/training_input/ion_torrent.png --profiles all_call/training_input/ion_torrent.tsv --old-profiles
#                   --true-values all_call/training_input/ion_torrent_truevals.json --params-file all_call/training_input/ion_torrent.params
