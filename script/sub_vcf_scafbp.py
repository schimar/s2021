#! /usr/bin/python
#
# This file subsets a vcf file to only contain variants from a tab-separated file of scaffold bp positions to include.

# Usage: ~/sub_vcf_scafbp.py in.vcf vars.scafbp.txt > out.vcf


from sys import argv

with open(argv[2], 'rb') as scafbp_file:
    scafbp = list()
    for line in scafbp_file:
        line = bytes.decode(line).strip('\n')
        tmp = ' '.join(line.split('\t')[:2])
        scafbp.append(tmp)#':'.join(line.split('\t')[:2]))

with open(argv[1], 'rb') as vcfile:
    for line in vcfile:
        line = bytes.decode(line).strip('\n')
        if line[0] == '#':
            print(line.split('\n')[0])
            #continue
        else:
            scaf, bp = line.split('\t')[:2]
            if ' '.join([scaf, bp]) in scafbp:
                print(line)
            else:
                continue


