#!/usr/bin/python
# Script simulates in-silico PCRs for Brucella sp.
# @author enrico1georgi@bundeswehr.org
# @version 1.1 (06.06.2021, _v2) --- comments added to improve usablity
# @version 1.2 (01.10.2021, _v3) --- added argparser to improve usablity

from Bio import SeqIO
import glob
import os
import sys
import argparse
import configparser

# **********
#    Help
# **********

parser=argparse.ArgumentParser(
    description='''This script simulates in-silico PCRs for given input files. Please refer to README.md and adapt the config.ini (in the same directory as the .py-file!) as desired.''',
    epilog="""If you need additional help, please feel free to send an informal message to enrico1georgi@bundeswehr.org.""")
args=parser.parse_args()



# ****************
#    Parameters
# ****************

# Parameters are changed in a config.ini file
config_error = False
config = configparser.ConfigParser()
if not os.path.exists('config.ini'):
    # create config.ini with default values if not exists
    config['PATHS'] = {'input_seqs_dir': './seqs/',
                       'input_primer_dir': './',
                       'output_dir': './output/'}
    config['OPTIONS'] = {'primer_filename': 'Primer_In_silico_PCR_Brucella.fasta',
                         'pcr_designation': 'MLVA16',
                         'MLVA16_verbose': 'no'}
    print('--> No config.ini file was detected. A new config.ini with default values is created.')
    print()
    with open('config.ini', 'w') as configfile:
        config.write(configfile)
else:
    # Read File
    config.read('config.ini')

# Check if file has distinct parameters
try:
    seqpath = config.get('PATHS', 'input_seqs_dir')
    print('--> Input seqfiles (FASTA format) are expected at:', seqpath)
    print()
except configparser.NoOptionError:
    print('ERROR: Input directory for seqfiles was not found in config.ini file.')
    config_error = True

try:
    print('--> File with primers is expected at:', config.get('PATHS', 'input_primer_dir'))
    print()
except configparser.NoOptionError:
    print('ERROR: Input directory for primer file was not found in config.ini file.')
    config_error = True

try:
    print('--> Output files will be saved at:', config.get('PATHS', 'output_dir'))
    print()
except configparser.NoOptionError:
    print('ERROR: Output directory was not found in config.ini file.')
    config_error = True

try:
    primerfile = config.get('PATHS', 'input_primer_dir') + config.get('OPTIONS', 'primer_filename')
    print('--> The following primer input file is considered:', primerfile)
    print()
except configparser.NoOptionError:
    print('ERROR: Primer filename was not found in config.ini file.')
    config_error = True

try:
    # Use OPTIONS > pcr_designation to select the desired analysis. For in-silico MLVA16 use 'MLVA16'.
    # Values may be: BCSP31, IS711, BruceLadder, MLVA16, rpoB, recA1, recA2, recA3
    pcr = config.get('OPTIONS', 'pcr_designation')
    print('--> The following in-silico PCR is performed:', pcr)
    output = pcr + '_results.fasta' if (pcr == 'BCSP31' or pcr == 'rpoB' or pcr == 'recA1' or pcr == 'recA2' or pcr == 'recA3') else pcr + '_results.txt'
    output = config.get('PATHS', 'output_dir') + output
    print('--> The following file will be created or updated:', output)
    print()
except configparser.NoOptionError:
    print('ERROR: A pcr designation like MLVA-16 was not found in config.ini file.')
    config_error = True

try:
    if pcr=='MLVA16':
        mlva16_verbose = config.getboolean('OPTIONS', 'MLVA16_verbose')
except configparser.ValueError:
    print('ERROR: A MLVA-16 verbose option was not found in config.ini file.')
    config_error = True

# All fasta files in the variable 'seqpath' are analyzed.
filelist = []
for filename in glob.glob(seqpath + '*.fasta'):
    filelist.append(filename)

if filelist:
    print('--> The following files are considered for analysis:')
    for filename in filelist:
       print(filename)
    print()
else:
    print('ERROR: No fasta seqfiles detected in given input folder.')
    config_error = True

# terminate program if some errors due to config file occured
if config_error == True:
    print()
    sys.exit('Program terminated for some trouble with config.ini file or if no input seq files are provided ...')



# ***********************************
#    Declarations of the functions
# ***********************************

# Returns the reverse complement sequence.
def reverseComplement(sequence):
    complement = {'a':'t','c':'g','g':'c','t':'a','n':'n'}
    return "".join([complement.get(nt.lower(), '') for nt in sequence[::-1]])

# Returns the amplicon size and writes the sequence to ampliconfile.
def getFragments(inputDNA, forward, reverse, plus, ampliconfile, pcr):
    intFragments = []
    split1 = inputDNA.split(forward)
    if len(split1)>1:
        for i in range (1,len(split1)):
            split2 = split1[i].split(reverseComplement(reverse).upper())
            if len(split2)>1:
                for j in range(0, len(split2)-1):
                    intFragments.append(len(forward)+len(split2[j])+len(reverse))
                    if (pcr == 'BCSP31' or pcr == 'rpoB' or pcr == 'recA1' or pcr == 'recA2' or pcr == 'recA3'):
                        if plus:
                            ampliconfile.write(str(forward + split2[j] + reverseComplement(reverse).upper() + '\n'))
                        else:
                            ampliconfile.write(str(reverse + reverseComplement(split2[j]).upper() + reverseComplement(forward).upper() + '\n'))
    return intFragments

# Ensures that both directions (3'->5' and 5'>3') are checked for a match.
def runPCR(inputDNA, forward, reverse, resultfile, pcr):
    return getFragments(inputDNA, forward, reverse, True, resultfile, pcr) + getFragments(inputDNA, reverse, forward, False, resultfile, pcr)

# Returns the primers from input file in a specific data format.
def extractPrimer(filename, designation):
    handlePrimer = open(filename, "rU")
    locus = []
    forward = []
    reverse = []
    # Starting value for counter: chr(65) returns A.
    k = 65
    for primer in SeqIO.parse(handlePrimer, "fasta"):
        if primer.id.split('_')[1] == 'Primer':
            if primer.id.split('_')[0] == designation:
                if (primer.id.split('_')[2][0] == 'F'):
                    if designation == 'MLVA16':
                        locus.append(primer.id.split('_')[3].split('-')[1])
                    elif designation == 'BruceLadder':
                        locus.append(chr(k))
                        k += 1
                    else:
                        locus.append(primer.id.split('_')[0])
                    forward.append(primer.seq)
                else:
                    reverse.append(primer.seq)
    handlePrimer.close()
    return list(zip(locus, forward, reverse))

# Calculates results for MLVA16.
def extractMLVAPattern(filename):
    handlePrimer = open(filename, "rU")
    locus = []
    for primer in SeqIO.parse(handlePrimer, "fasta"):
        p = primer.id.split('_')
        if p[0] == 'MLVA16':
            if (p[2][0] == 'F'):
                repeatsize = int(p[4][:-2])
                refsize = int(p[5][:-2])
                refunit = int(p[6][:-1])
                offset = refsize - (repeatsize * refunit)
                locus.append((p[3].split('-')[1], offset, repeatsize))
    handlePrimer.close()
    return locus

# Performs in-silico electrophoresis.
def runGel(myDict):
    fragments = []
    for key in myDict:
        fragments += myDict[key]
    return sorted(set(fragments))

# ******************
#    Main program
# ******************

myPrimers = extractPrimer(primerfile, pcr)
myMLVAPattern = extractMLVAPattern(primerfile)
if mlva16_verbose:
    print('Following loci with offset and repeatsize are given:')
    for locus in myMLVAPattern:
        print(locus)
    print()

output_handle = open(output, "w")

# Writes header to MLVA16 output file. For each locus there will be three columns: A=ampliconSize, C=repeatCount, M=modulo.
if pcr=='MLVA16':
    output_handle.write('StrainDesignation' + '\t')
    for locus in myMLVAPattern:
        output_handle.write(locus[0] + 'A' + '\t')
        output_handle.write(locus[0] + 'C' + '\t')
        output_handle.write(locus[0] + 'M' + '\t')
    print()
    output_handle.write('\n')

for file in filelist:
    myBands = {}
    print('***')
    print(os.path.basename(file).split('.')[0])
    lineend = '\n' if (pcr == 'BCSP31' or pcr == 'rpoB' or pcr == 'recA1' or pcr == 'recA2' or pcr == 'recA3') else '\t'
    output_handle.write('>' + os.path.basename(file).split('.')[0] + lineend)
    print('***')
    handle = open(file, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        if (len(record)>10000):
            #print(record.id, "(%i bp)"  % len(record))
            for i in range (0, len(myPrimers)):
                if myPrimers[i][0] in myBands:
                    myBands[myPrimers[i][0]] += runPCR(record.seq, myPrimers[i][1], myPrimers[i][2], output_handle, pcr)
                else:
                    myBands[myPrimers[i][0]] = runPCR(record.seq, myPrimers[i][1], myPrimers[i][2], output_handle, pcr)
    handle.close()

    print()
    if pcr!='MLVA16':
        print('Amplicons per Locus:')
        print('********************')
        for key in sorted(myBands):
            if pcr == 'IS711':
                output_handle.write(str(myBands[key]) + '\n')
            else:
                print(key, myBands[key])
    else:
        print('MLVA16: Repeats per Locus:')
        print('**************************')
        for locus in myMLVAPattern:
            if not myBands[locus[0]]:
                ampliconSize = 0
                repeatCount = 0
                modulo = 0
            else:
                ampliconSize = int(myBands[locus[0]][0])
                repeatCount = (ampliconSize - locus[1]) // locus[2]
                modulo = (ampliconSize - locus[1]) % locus[2]
            output_handle.write(str(myBands[locus[0]]) + '\t' + str(repeatCount) + '\t' + str(modulo) + '\t')
            print(locus[0], '...', str(repeatCount))
        output_handle.write('\n')


    print()

    if pcr != 'MLVA16':
        print('In-silico-Gel:')
        print('*************')
        print(runGel(myBands))
        if pcr == 'BruceLadder':
            output_handle.write(str(runGel(myBands)) + '\n')
        print()

output_handle.close()
