import sys
from datetime import datetime
from all_call import train, input
import scipy.optimize
import pandas as pd

MIN_HITS = 10


def start_it():
    """
    Write start time.
    :return: datetime - start time
    """
    # print the time of the end:
    # print time of the start:
    start_time = datetime.now()
    print('AlleG = Allele Gatherer (for Dante)', file=sys.stderr)
    print('AlleG Starting  : {start:%Y-%m-%d %H:%M:%S}'.format(start=start_time), file=sys.stderr)
    return start_time


def end_it(start_time):
    """
    Write end time.
    :param start_time: datetime object - start time o f the program
    :return: None
    """
    # print the time of the end:
    end_time = datetime.now()
    print('AlleG Stopping  : {finish:%Y-%m-%d %H:%M:%S}'.format(finish=end_time), file=sys.stderr)
    print('Total time of run : {duration}'.format(duration=end_time - start_time), file=sys.stderr)


def change_suffix(filename, suffix):
    """
    Change suffix into a new one.
    :param filename: str - filename
    :param suffix: str - string to append
    :return: str - new filename with new suffix
    """
    return filename[:filename.rfind('.')] + suffix


### main program:

if __name__ == "__main__":
    # start
    start_time = start_it()

    # load arguments
    args = input.load_arguments()
    print("Arguments read successfully...", file=sys.stderr)

    # get all data files:
    configs, profiles, true_vals = input.crawl_dante(args.dir_structure)

    # merge the profile files (and output them)
    merged_profile = input.merge_profiles(profiles, args.output_profile)
    merged_profile = merged_profile.astype(int)

    # load or merge the true_vals files (and output them)
    if args.input_true is None:
        merged_true = input.merge_profiles(true_vals, args.output_true)
    else:
        merged_true = pd.read_csv(args.input_true, sep='\t', index_col=0, parse_dates=True, engine='python')
        merged_true.columns = [0, 1]

    # fix 'B' and 'E'
    def convert_to_int(x):
        if x == 'B':
            return 0
        if x == 'E':
            return 0
        return int(x)

    merged_true = merged_true.applymap(convert_to_int)

    # write progress
    if merged_profile is not None and merged_true is not None:
        print("Training data read successfully:\n"
              "  %3d samples read\n"
              "  %3d true values read" % (len(merged_profile.index), len(merged_true.index)), file=sys.stderr)
    else:
        print("ERROR: directory %s does not contain Dante subdirectory structure: <dir_name>/config.yaml,all_profiles.txt,all_profiles.true" % args.dir_structure)
        end_it(start_time)
        exit(-1)

    if not args.prepare:

        # filter profiles good for training:
        test_samples = filter(lambda sample: train.good_for_sampling(sample, merged_true.get_value(index=sample, col=0), merged_true.get_value(index=sample, col=1),
                                                                     merged_profile, single=False), merged_true.index)
        single_samples = filter(lambda sample: train.good_for_sampling(sample, merged_true.get_value(index=sample, col=0), merged_true.get_value(index=sample, col=1),
                                                                       merged_profile, single=True), merged_true.index)
        print("Extracted:\n"
              "  %3d good training samples\n"
              "  %3d single values training samples" % (len(test_samples), len(single_samples)), file=sys.stderr)

        # samples_struct is list(dict(int:ndarray(int))) - list of all samples with dicts of alleles and their respective arrays
        all_samples = test_samples + single_samples
        samples_struct = map(lambda sample: train.extract_alleles(sample, merged_true, merged_profile, args.verbosity_level > 1), all_samples)

        if len(all_samples) >= 0:
            # normalize (for model training):
            if args.verbosity_level >= 2:
                print(samples_struct)
            print("Normalizing...", file=sys.stderr)
            samples_struct_norm = train.norm_profiles(samples_struct)

            # train the model:
            print("Training of the binomial model...", file=sys.stderr)
            start = [0.005, 0.0005, 0.01, 0.005]
            params = scipy.optimize.fmin(train.comparison, start, args=(samples_struct, input.fit_functions[args.fit_function]), xtol=0.000001, maxfun=1000000, maxiter=1000000)

            # count the occurrences:
            nums, nums_prop = train.count_occurrences(samples_struct, len_repeating=3)

            # train read count drop:
            print("Training of linear read count drop (absolute)...", file=sys.stderr)
            params_read_drop = train.train_read_drop_abs(nums)
            if params_read_drop is None:
                params_read_drop = input.DEFAULT_READ_DROP

            print("Training of linear read count drop (relative)...", file=sys.stderr)
            params_read_drop_rel = train.train_read_drop_rel(nums_prop)
            if params_read_drop_rel is None:
                params_read_drop_rel = input.DEFAULT_READ_DROP

            # save the params
            print("Writing the parameters to output (%s)..." % args.output_params if args.output_params is not None else "stdout", file=sys.stderr)
            if args.output_params is not None:
                with open(args.output_params, 'w') as f:
                    train.write_params(f, params, params_read_drop, params_read_drop_rel, args.fit_function, print_all=True)
            else:
                train.write_params(sys.stdout, params, params_read_drop, params_read_drop_rel, args.fit_function, print_all=True)

            # create new config files (if needed)
            if args.config_dir is not None:
                print("Creating new config files in %s..." % args.config_dir, file=sys.stderr)
                for config in configs:
                    # read the config and output new one
                    input.update_config(config, args.config_dir, args.output_params)

    # write end time
    end_it(start_time)

    # # display?
    # if args.model_fig is not None:
    #     print("Generating model figures...", file=sys.stderr)
    #     rng = np.arange(len(samples_struct[0].values()[0]))
    #     model_tmp = train.model_template(rng, params, input.fit_functions[args.fit_function])
    #     train.display(model_tmp, samples_struct_norm, args.model_fig)
    #
    # # display?
    # if args.model_fig is not None:
    #     print("Generating read count drop fit...", file=sys.stderr)
    #     last_dot = args.model_fig.rfind('.')
    #     filename_curr = "%s_reads%s" % (args.model_fig[:last_dot], args.model_fig[last_dot:])
    #     x, y = train.flatten(nums)
    #     train.plot_occurrences(x, y, filename_curr, params_read_drop)
    #     filename_curr = "%s_reads_relative%s" % (args.model_fig[:last_dot], args.model_fig[last_dot:])
    #     x5, y5 = train.flatten(nums_prop)
    #     train.plot_occurrences(x5, y5, filename_curr, params_read_drop_rel)
    #
    #
