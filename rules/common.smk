import os
import pandas as pd

###### Config file and sample sheets #####


configfile: 
  "config/config.yaml"

# load sample info 
samples_information = pd.read_csv("samples.txt", sep='\t', index_col=False)

sample_names = list(samples_information['sample'])
sample_locations = list(samples_information['location'])
samples_set = zip(sample_names, sample_locations)
samples_dict = dict(zip(sample_names, sample_locations))


###### helper functions ######

def getFqHome(sample):
  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))

def getTrmHome(sample):
  return(list(os.path.join('trm', "{0}_{1}_trm.fq.gz".format(sample,pair)) for pair in ['R1','R2']))




