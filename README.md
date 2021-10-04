# in-silico-typing
## Purpose
This script is an *in-silico* approach to type WGS datasets or to perform an *in-silico* PCR. Input files are expected to be (multi)FASTA files, one file per strain. Main purpose is to do a simple MLVA16 of *Brucella* sequences (that's why 'Brucella melitensis bv 1 str 16M.fasta' is bundled for testing), however other *in-silico* PCRs are also possible.

## Installation

The python script requires a [Python 3 release](https://www.python.org/downloads/) and some modules (Bio, glob, os, sys, argparse, configparser) to run.

1. Download the repository to local file (Code > Download ZIP) and unzip data.
1. Start script for the first time. A config.ini file with default values is created.
```sh
python3 In-silico-PCR_v3.py
```
1. *Optional*: Install missing dependencies as indicated by the error message, e.g. Bio
```sh
pip3 install biopython
```
1. *Optional*: Modify values in the config.ini to your needs.
1. Put your FASTA files into the input folder as defined in config.ini file. For testing of you can use the bundled 'Brucella melitensis bv 1 str 16M.fasta' file.
1. The file containing all the primer/probe sequneces should not be in the seq input folder. Usually it is located in the same folder as the main py-script. As a start use the bundled 'Primer_In_silico_PCR_Brucella.fasta' file. However, you can use your own of course.
1. Restart the script. Output is saved in output folder as defined in config.ini file.
```sh
python3 In-silico-PCR_v3.py
```

## Mismatches
Up to the current version the script is not intended to handle mismatches within conserved primer sequences. So, if you get no results for a specific locus, it may be due to mismatches in your input sequence. One approach to resolve this is to append loci with SNPs to the primer file (see 'Primer_In_silico_PCR_Brucella.fasta' for an example).

## Help
If you need additional help, please feel free to send an informal message to enrico1georgi@bundeswehr.org.