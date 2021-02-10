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

# load sample info for all ant samples (the 85 + 23 = 108) 
samples108 = pd.read_csv("samples108.txt", sep='\t', index_col=False)

sample_names108 = list(samples108['sample'])
sample_locations108 = list(samples108['location'])
samples_set108 = zip(sample_names108, sample_locations108)
samples_dict108 = dict(zip(sample_names108, sample_locations108))


###### helper functions ######

def getFqHome(sample):
  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))

def getTrmHome(sample):
  return(list(os.path.join('trm', "{0}_{1}_trm.fq.gz".format(sample,pair)) for pair in ['R1','R2']))


#def getSamHome(sample):
#  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))

