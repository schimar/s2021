#! /usr/bin/python
#
# This script filters a vcf file based on overall sequence coverage, number of
# non-reference reads, number of alleles, and reverse orientation reads.
# See below for default values, and to change them, if necessary. Additionally,
# note that currently, the number of retained loci is not being written at the end
# of the file.
# Usage: ./vcfFilter.py <variant file>.vcf > outfile.vcf

from sys import argv
import re
import shutil
#import tempfile

### stringency variables, edit as desired
minCoverage = 128 # minimum number of seqs; DP
minAltRds = 4 # minimum number of sequences with the alternative allele; AD
fixed = [0.0, 1.0] # removes loci fixed for alt; RAF
mapQual = 50 # minimum mapping quality
typls = ["DEL", "INS", "SUB"]


# added 10/02/17
#minBqrs = -8 # minimum absolute value of the base quality rank sum test; BaseQRankSum
#minMqrs = -12.5 # minimum absolute value of the mapping quality rank sum test; MQRankSum
#minRprs = -8 # minimum absolute value of the read position rank sum test; ReadPosRankSum
#minQd = 2 # minimum ratio of variant confidenct to non reference read depth; QD
#maxFish = 60 # maximum phred-scaled p-value using Fisher's Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls; FS

n_seqs_retained = int()
with open(argv[1], 'rb') as file:
    for line in file:
        if line[0] == '#':
            #continue
            print line.strip('\n')
        else:
            lspl = line.strip('\n').split('\t')
            test = lspl[6]
            altAll = lspl[4]
            if test == 'PASS':
                dp = int(re.findall('DP=[0-9]+', line)[0].split('=')[1])
                ac = int(re.findall('AD=[0-9]+', line)[0].split('=')[1])
                raf = float(re.findall('RAF=[0.0-9.0]+', line)[0].split('=')[1])
                mqs = int(re.findall('MQS=[0-9]+', line)[0].split('=')[1])
                mqm = int(re.findall('MQM=[0-9]+', line)[0].split('=')[1])
                ls = int(re.findall('LS=[0-9]+', line)[0].split('=')[1])
                typ = str(re.findall('TYP=[A-Z]+', line)[0].split('=')[1])
                ppc = int(re.findall('PPC=[0-9]+', line)[0].split('=')[1])
                #
                if (raf not in fixed and dp >= minCoverage and typ in typls and mqs >= mapQual and ac >= minAltRds):
                    print line.strip('\n')


    file.close()

#print '#Retained %i variable loci' % n_seqs_retained
