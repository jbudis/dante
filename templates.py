OUT_STATS = '''ins={inserts}
del={deletes}
base={bases}
reads={remained}'''

ANNOTATION = '''>{annotation.read.name}
pair={pairness}
complement={annotation.read.complement}
likelihood={likelihood}
bases={annotation.module_bases}
repetitions={annotation.module_repetitions}
errors={annotation.n_insertions}I, {annotation.n_deletions}D, {annotation.n_mismatches}M = {errors}E
{annotation.ann_sequence}
{annotation.ann_motif}
{annotation.ann_module}
{annotation.error_line}
'''

DANTE_DESCRIPTION = '''
    DANTE = Da Amazing NucleoTide Exposer
    ---------------------------------------
Sequence annotator based on Hidden Markov Models (HMM). All of the 
configuration is done via a (required) configuration file in YAML format.
A default configuration file with inline explanations is a necessary part 
of Dante distribution (usually in input_data/default.yaml). If an option 
is not present in the configuration file, this default will be used. 
Reads from the input file(s) supplied in the configuration file are 
pre-filtered, annotated by HMM, and post-filtered. The base quality scores 
are used to determine the emit probability of each base. DANTE creates a 
directory for each of the annotated motifs and a summary report in HTML 
(report.html). For each post-filter options (X), the motif directory contains
the following files:
1. annotations_X.txt -- correctly annotated reads, that passed all filters
2. filtered_primer_X.txt -- correctly annotated reads, that have only one 
                flanking region
3. filtered_X.txt -- annotated reads, that didn't pass the post-filtering
4. repetitions_X.txt -- number of correctly annotated reads with both 
                flanking regions and their number of STR repetitions
5. repetitions_grey_X.txt -- number of correctly annotated reads with one 
                flanking region and their number of STR repetitions
6. alignment_X.txt -- alignment of number of correctly annotated reads with 
                both flanking regions
7. all_call_X.txt -- genotype and its confidence
8. pcolor_X.{pdf/png} -- plot of all tested genotype options and their 
                likelihoods
'''

START_TIME = '''DANTE = "Da Amazing NucleoTide Exposer"
DANTE Starting    : {start: %Y-%m-%d %H:%M:%S}'''
