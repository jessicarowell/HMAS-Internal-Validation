# HMAS Internal Validation Pipeline

A pipeline for computing validation based on the positive control primers.
It is intended to be run after the HMAS QC Pipeline.

Started by [@jessicarowell](https://github.com/jessicarowell)
Active development by [@jessicarowell](https://github.com/jessicarowell) 

## TOC
* [Description](#description)
* [Requirements](#requirements)
* [INSTALL](#install)
* [USAGE](#usage)
* [Contributing](#contributing)
* [Future Plans](#future-plans)
* [Resources](#resources)

## Description

This pipeline takes 4 input files. The first 3 are outputs from HMAS QC Pipeline.
* A fasta file: this file contains the high-quality consensus reads that passed filtering steps in the HMAS QC Pipeline 
* A names file: this text file matches each consensus reads with all the identical reads
* A group file: this text file matches each read that passed QC to it's associated primer pair
* A primer design file: this is a fasta-formatted file with primer pair name as the header, and its associated amplicon as the sequence.

The design file contains amplicons specifically designed as positive controls for our 823 AMD primers.  Those
823 primers were designed to walk a set of AMD genes of interest to enterics.  For each primer pair, we design
ed an amplicon of the same length and GC content using Coliphage phi-X17, Wolbachia pipientis, or Streptomyces
 coelicolor (depending on target GC content).

In a perfect world, we expect the primers to bind and amplify the amplicons they were designed to target. So
in our final, high-quality fasta file, we hope the reads and their associated primer pair will match the desig
n file.  A blast search against this **reference database** built from the amplicons in our design file helps
us evaluate how well our primers hit their expected targets.


The pipeline runs the following steps:
1. It creates a blast database from the primer design file. (this is the **reference database**)
2. It blasts the final, post-QC fasta file against this reference database. It allows a maximum of 20 hits.
3. Create a SQL database from the group file {read: primer pair name}.  We will search against this database.
4. Get and store the number of consensus reads for each unique reads in the final, post-QC fasta file.
5. Collate all of this information so we have the primer pair and number of consensus reads for each unique read in our final fasta file. 
6. Merge our final reads with the blast results to find out the following information:
   - How many reads perfectly match with their intended primer pairs in the design file?
   - How many reads don't match anything in the design file?
   - How many reads match more than just their intended primer pairs?


## Requirements

1. Python 3 or higher. Download python [here](https://www.python.org/downloads/). 

2. You must have pandas installed. 

To be completed...


## INSTALL

(Note: I haven't elaborated here because these instructions will change when we containerize the pipeline.)

1. Copy the Github repository to a folder  
`git clone https://github.com/jessicarowell/HMAS-Internal-Validation.git` 

2. Add `hmas_validation.py` to your $PATH


## USAGE

To be elaborated...

`hmas_validation.py -h`

`hmas_validation.py -f twist.final.fasta -r amr_design_primers.fasta -g twist.final.groups -n twist.final.names -o amr2020`


## Contributing

(Note: this section might also change depending on how we package this.)

Please feel free to fork this repo, make improvements, and share them with me.

Please also post any issues you encounter to Github and I'll be sure to look into them as soon as I can.


## Future Plans

We plan to containerize this pipeline in the future.

## Resources

