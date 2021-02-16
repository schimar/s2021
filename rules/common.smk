import os
import pandas as pd
from glob import glob

###### Config file and sample sheets #####


configfile: 
  "config/config.yaml"

# load sample info 
samples23 = pd.read_csv(config['samplesSub'], sep='\t', index_col=False)

sample_names = list(samples23['sample'])
sample_locations = list(samples23['location'])
samples_set = zip(sample_names, sample_locations)
samples_dict = dict(zip(sample_names, sample_locations))

# load sample info for all ant samples (the 85 + 23 = 108) 
ids = pd.read_table(config['ids'], sep='\t', index_col=False)
ids['id_nest'] = ids['id'] + '_' + ids['nest']
iddict = ids[['sample', 'id_nest']].set_index('sample').T.to_dict('list')



###### helper functions ######

def getFqHome(sample):
  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))

def getTrmHome(sample23):
  return(list(os.path.join('trm', "{0}_{1}_trm.fq.gz".format(sample23,pair)) for pair in ['R1','R2']))

#def getBamHome(ids):
#    return(os.path.join('bbmap', "{ids}") #.format(ids=wildcards.ids)))
 
#def getBamHome(ids, ext):
#  return(list(os.path.join('bbmap', '.'.join([ids[ids.id_nest == wildcards.id].sample, ext]))))
  
#def getSamHome(sample):
#  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))




def renameBamBai(path, namedict):
  bamls = glob(os.path.join(path, '*.bam'))
  bamls.sort()
  bails = glob(os.path.join(path, '*.bam.bai'))
  bails.sort()
  for i, fname in enumerate(bamls):
    os.rename(fname, os.path.join(path, '.'.join([namedict[fname], 'bam'])))
    os.rename(bails[i], os.path.join(path, '.'.join([namedict[fname], 'bam.bai'])))

 
# renameBamBai('bbmap', iddict)
