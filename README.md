# Dante
 Dante ("Da Amazing NucleoTide Exposer") is an alignment-free algorithm for genotyping STR alleles based on NGS reads originating from the STR locus of interest. The method accounts for natural deviations from the expected sequence, such as variation in the repeat count, sequencing errors, ambiguous bases, and complex loci containing several different motifs.

Dante profiles sequenced DNA fragments against STR motifs defined in a user-friendly configuration file. For each of the target motifs, [the final report](http://158.195.68.48/dante/example/report.html) contains the following: 
- pair of the most probable genotypes with confidence scores
- observed number of repetitions, probabilities of possible genotypes in a color plot
- sequence logos corresponding to called alleles

Reported figures provide an evidence for expanded allele, which is too long to be captured by a single NGS read as well as for allelic single point mutations, small insertions, and deletions that may be relevant for a diagnostic evaluation.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine. 

### Prerequisites

Dante is developed and tested on Python 2.7 with libraries determined in the requirements.txt file: 

```
numpy==1.11.0
multiprocess==0.70.5
pysam==0.11.2.2
regex==2017.9.23
numba==0.35.0
pandas==0.19.1
scipy==0.18.1
matplotlib==1.5.1
PyYAML==3.12
```

All libraries may be easily installed using pip.

```
pip install <dependency>
```

### Installation

Dante requires only Python 2.7 interpreter with libraries listed in section Prerequisites. After cloning the repository, you may check out by: 
```
python dante.py --help
```

### Example dataset

Cloned repository contains set of example reads in Fastq format. You may try Dante with: 
```
cd <dante_dir>/example
python ../dante.py config.yaml
```

Dante will annotate and genotype 7 selected STR motifs using reads in files example_R1.fastq.gz and example_R2.fastq.gz. Results will be stored in the <dante_dir>/example/report directory. Summary web page with genotyping calls and supporting figures would be generated at <dante_dir>/example/report.html. 

## Config file

Dante requires (and accepts) only one input: configuration file in YAML format. All of the configuration options are described directly in the default configuration file ("input_data/default.yaml") and in the example config file ("example/config.yaml").

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

Dante uses the same stutter parameters for all motifs (which we will probably change in later versions), so we recommend to use similarly behaving motifs together (motifs that have the same or similar length of repetition). The training is usually done under a minute and we recommend to use maximal possible number of samples to achieve best results. 

