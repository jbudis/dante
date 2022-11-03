# Dante
 Dante ("Da Amazing NucleoTide Exposer") is an alignment-free algorithm for genotyping STR alleles based on NGS reads originating from the STR locus of interest. The method accounts for natural deviations from the expected sequence, such as variation in the repeat count, sequencing errors, ambiguous bases, and complex loci containing several motifs.

Dante profiles sequenced DNA fragments against STR motifs defined in a user-friendly configuration file. For each of the target motifs, [the final report](http://158.195.68.48/dante/example/report.html) contains the following: 
- pair of the most probable genotypes with confidence scores
- observed number of repetitions, probabilities of possible genotypes in a color plot
- sequence logos corresponding to called alleles

Reported figures provide an evidence for an expanded allele, which is too long to be captured by a single NGS read as well as for allelic single point mutations, small insertions, and deletions that may be relevant for a diagnostic evaluation.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine. 

### Installation

Dante is developed and tested on Python 3.7 with libraries determined in the `conda_env.yaml`. Create the conda environment as following:
```bash
conda env create -n dante -f conda_env.yaml
```

Furthermore, if you need to use the config converter program `make_config.py`, you need to additionally install `ensembl_rest` python package (not available directly through conda):

```bash
pip install ensembl_rest==0.3.3
```

If you want to use `create_motif_report.py` script to aggregate reports by motif, you will need to install additional `beautifulsoup4` dependency:

```bash
conda install beautifulsoup4
```

(To do this on some very old systems, you need to install openssl to your conda environment with `conda install -c conda-forge openssl==3.0.5`.)

### Example dataset

Cloned repository contains set of example reads in Fastq format. You may try Dante with: 
```
cd <dante_dir>/example
python ../dante.py config.yaml
```

Dante will annotate and genotype 7 selected STR motifs using reads in files example_R1.fastq.gz and example_R2.fastq.gz. Results will be stored in the <dante_dir>/example/report directory. Summary web page with genotyping calls and supporting figures would be generated at <dante_dir>/example/report.html. 

## Config file

Dante requires (and accepts) only one input: configuration file in YAML format. All the configuration options are described directly in the default configuration file ("input_data/default.yaml") and in the example config file ("example/config.yaml").

## Methods

Dante accepts various file formats for input NGS reads:
- Fasta(.gz)
- Fastq(.gz)
- SAM/BAM
- .txt output from dante (debug purposes) 

Input sequence reads are pre-filtered to keep only those that have potential to originate from the STR locus of interest. Reads that passed the pre-filtering are annotated by the profile HMM specifically tailored for the STR motif. The annotation process outlines segments of read sequence that are most likely parts of the motif. Reads with consistent annotation are stored and the number of STR repetitions is supplied to the genotyping process.

The genotyping process takes annotated reads with full or partial overlap of the STR sequence and infers the most probable underlying genotype. It assumes single or two alleles of a diploid organism, but takes into account possibility of the shared, homozygous genotype. Prediction of expanded alleles that cannot be fully captured due to technical limitations of NGS sequencing is supported by reads with partial overlaps. 

## Training of stutter parameters (allcall.py)

To train the PCR stutter insertion and deletion probability, we provide a separate program called allcall.py. This program requires several samples (directories) with already completed Dante annotation. For example, assume we have 3 config files, each corresponding to one sample. After running Dante once for each config, we have 3 directories: "Dante_outputs/sample1", "Dante_outputs/sample2", "Dante_outputs/sample3" (defined in corresponding config files as general: output_dir: Dante_outputs/sampleX). Each of these directories contains results of Dante for one sample and all motifs defined in the corresponding config. 

Then you should run allcall.py on a directory where these results reside (Dante_outputs) in preparation mode as: 
```bash
python allcall.py Dante_outputs --prepare
```  

This call generates files: 
- Dante_outputs/all_profiles.txt -- all repetition profiles of all samples and all motifs
- Dante_outputs/all_profiles.true -- all predicted results (by Dante) 

At this point you should manually fix results, that are not correctly predicted in the file "Dante_outputs/all_profiles.true" (if everything is correct, you do not need training at all ;-) ). After that, you are ready to train the stutter parameters (and optionally output new config files for Dante based on the previous config files with parameter --config-dir):

```bash
python allcall.py Dante_outputs --input-true=Dante_outputs/all_profiles.true [--config-dir=Dante_outputs]
```

This call creates a file "Dante_outputs/params.txt" with new stutter parameters and (optionally) new config files for Dante. These new config files differ from the previous ones by the stutter parameter file (the newly generated is used), output directory ("_retrained" is appended), and input files (output of the previous run of Dante is used as input). The change of input guarantees that Dante calls with these config files are much faster, since only reads  correctly annotated in the previous phase are processed. 

Finally, in case you did not mess with generated config files, directories "Dante_outputs/sample1_retrained", "Dante_outputs/sample2_retrained", "Dante_outputs/sample3_retrained" contain Dante results with trained stutter parameters.

Dante uses the same stutter parameters for all motifs (which we will probably change in later versions), so we recommend to use similarly behaving motifs together (motifs that have the same or similar length of repetition). The training is usually done under a minute, and we recommend to use maximal possible number of samples to achieve best results. 

## Creating motif reports

To collect results from multiple reports and group them by motif, we provide a script `create_motif_report.py`. It has one required and one optional argument: `input dir` and `output dir`. Input dir is a path to folder with reports from samples. There can be other files as well, the script filters out all files that don't match `*.html`. Output dir is a path to directory, where the motif reports will be generated. If the path isn't specified, the script generates reports to `example/motif_report`.

After running Dante on example dataset, you may run this script as:

```bash
python create_motif_report.py example/report [example/motif_report]
```

Script generates separate report for each unique motif in Dante reports. Script uses the file name of the report to differentiate the source of result in table, so it is recommended to rename the reports before starting this script, so no two reports have the same name.

