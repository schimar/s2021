#! /usr/bin/python

# This script reads a vcf file and converts it into zarr format, for faster access from scikit-allel.

# Usage: ./vcf2zarr.py <file.vcf>

import pandas as pd
import numpy as np
import sys
import allel as al
import zarr

vcfPath = sys.argv[1]

variants = al.read_vcf(vcfPath, numbers= {'GT': 2, 'ALT': 1}, fields= '*')	# with all (*) fields read



zarrPath = str(vcfPath.split('.')[0] + '.zarr')
al.vcf_to_zarr(vcfPath, zarrPath, fields='*', log=sys.stdout, overwrite=True)



