#!/usr/bin/python
# Script simulates in-silico PCRs for Brucella sp.
# @author enrico1georgi@bundeswehr.org
# @version 1.0 (13.07.2016)

from Bio import SeqIO
import glob
import os

# ****************
#    Parameters
# ****************

mode = 'Own'    # Ref, Own
pcr = 'MLVA16'  # BCSP31, IS711, BruceLadder, MLVA16, rpoB, recA1, recA2, recA3

# ***********
#    Files
# ***********

primerfile = 'Primer_In_silico_PCR_Brucella.fasta'
output = pcr + '_results.fasta' if (pcr == 'BCSP31' or pcr == 'rpoB' or pcr == 'recA1' or pcr == 'recA2' or pcr == 'recA3') else pcr + '_results.txt'
seqpath = '..\\..\\..\\+++ NGS +++\\RefSeq_Brucella\\' if mode == 'Ref' else '.\\seq\\'

filelist = []
for filename in glob.glob(seqpath + '*.fasta'):
    filelist.append(filename)



# ****************************
#    Funktionsdeklarationen
# ****************************

def reverseComplement(sequence):
    complement = {'a':'t','c':'g','g':'c','t':'a','n':'n'}
    return "".join([complement.get(nt.lower(), '') for nt in sequence[::-1]])

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

def runPCR(inputDNA, forward, reverse, resultfile, pcr):
    return getFragments(inputDNA, forward, reverse, True, resultfile, pcr) + getFragments(inputDNA, reverse, forward, False, resultfile, pcr)

def extractPrimer(filename, designation):
    handlePrimer = open(filename, "rU")
    locus = []
    forward = []
    reverse = []
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
    
def runGel(myDict):
    fragments = []
    for key in myDict:
        fragments += myDict[key]
    return sorted(set(fragments))

# *******************
#    Hauptprogramm
# *******************

myPrimers = extractPrimer(primerfile, pcr)
myMLVAPattern = extractMLVAPattern(primerfile)
for locus in myMLVAPattern:
    print(locus)
print()

# ***** Check for file io errors

output_handle = open(output, "w")

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
    print('Amplicons per Locus:')
    print('********************')
    if pcr!='MLVA16':
        for key in sorted(myBands):
            if pcr == 'IS711':
                output_handle.write(str(myBands[key]) + '\n')
            else:
                print(key, myBands[key])
    else:
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


