#! /usr/bin/python

# This script reads a vcf file and converts it into zarr format, for faster access.

# Usage: ./vcf2zarr.py <file.vcf>

import pandas as pd
import numpy as np
import sys
import zarr

vcf = sys.argv[1]

print(vcf)

