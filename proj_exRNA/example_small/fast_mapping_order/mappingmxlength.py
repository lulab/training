#! /usr/bin/env python
from tqdm import tqdm
import argparse
import numpy as np
import time
import sys,os,errno,gc
import numba
import itertools
import h5py
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='ind')  
parser.add_argument('-o', dest='output_file')
args = parser.parse_args()




originalnames = np.array(['length','lncRNA', 'miRNA', 'mRNA', 'piRNA', 'snoRNA', 'snRNA', 'srpRNA',
       'tRNA', 'vaultRNA', 'Y_RNA', 'tucpRNA'])
rnanames = np.array(['length', 'miRNA', 
            'Y_RNA', 'piRNA', 'srpRNA', 
            'snoRNA', 'snRNA', 'tRNA','vaultRNA',
            'tucpRNA', 'lncRNA','mRNA' ])
xsorted = np.argsort(originalnames)
ypos = np.searchsorted(originalnames[xsorted], rnanames)
transind  = xsorted[ypos]

lengthind = np.arange(16,151).astype('int')

filename = np.loadtxt('nocutlengthfilename.txt'
                      ,dtype='str')[int(args.ind)]  
samplearrmx = np.array(pd.read_table('/Share/home/caojingyi/exRNA/process/16.one_map/no_cut/02.hg38/'+filename)
                       .iloc[:,:])[:,transind]


print ('sample matrix loaded')
print (samplearrmx.shape)


mirnacountafter = np.where(samplearrmx[:,1]==0) 
mirnacountneed = np.where(samplearrmx[:,1]!=0) 
countind1,counts1 =np.unique(samplearrmx[mirnacountneed,0],return_counts=True)
countind1 = countind1.astype('int')
positionind = np.where(np.in1d(lengthind,countind1))
mirnalength= np.zeros([135])
print (counts1.shape)
print (positionind[0].shape)
mirnalength[positionind] = counts1


group1needarr = samplearrmx[mirnacountafter[0]][:,2:5]#np.concatenate((np.array([0]),np.arange(2,5)))]
group1countafter = np.where(np.sum(group1needarr,axis=1)==0)
group1needind = np.where(np.sum(group1needarr,axis=1)!=0)
group1needarr = group1needarr[group1needind]
lengthgroup1 = samplearrmx[:,0][mirnacountafter[0]][group1needind]

group2needarr = samplearrmx[mirnacountafter[0]][group1countafter][:,5:9]#np.concatenate((np.array([0]),np.arange(5,9)))]
group2countafter = np.where(np.sum(group2needarr,axis=1)==0)
group2needind = np.where(np.sum(group2needarr,axis=1)!=0)
group2needarr = group2needarr[group2needind]
lengthgroup2 = samplearrmx[:,0][mirnacountafter[0]][group1countafter[0]][group2needind]

group3needarr = samplearrmx[mirnacountafter[0]][group1countafter][group2countafter][:,9:12]#np.concatenate((np.array([0]),np.arange(9,12)))]
group3needind = np.where(np.sum(group3needarr,axis=1)!=0)
group3needarr = group3needarr[group3needind]
lengthgroup3 = samplearrmx[:,0][mirnacountafter[0]][group1countafter[0]][group2countafter[0]][group3needind]


def get_single_column_length(lengthgroup,groupneedarr,i):
    groupcol = (groupneedarr !=0).argmax(axis=1)
    groupcountind = np.unique( lengthgroup[np.where(groupcol ==i)],return_counts=True)
    groupcountind_ = groupcountind[0].astype('int')
    returnlength = np.zeros([135])
    for i in range(135):
        returnlength[i] = groupcountind[1][np.where(groupcountind_ ==(i+16))] if np.where(groupcountind_ ==(i+16))[0].shape[0]!=0 else 0
    #positionind = np.where(np.in1d(lengthind,groupcountind_))[0]
    #returnlength = np.zeros([35])
    #returnlength[positionind] = groupcountind[1]
    return returnlength

print (group1needarr.shape,group2needarr.shape,group3needarr.shape)

counts1 = np.ndarray([1,3,135]).astype('int')

group1countlength1 = get_single_column_length(lengthgroup1,group1needarr,0)
group1countlength2 = get_single_column_length(lengthgroup1,group1needarr,1)
group1countlength3 = get_single_column_length(lengthgroup1,group1needarr,2)
counts1 = np.concatenate((group1countlength1,group1countlength2
                            ,group1countlength3)).reshape(3,-1)

counts2 = np.ndarray([1,4,135]).astype('int')

group2countlength1 = get_single_column_length(lengthgroup2,group2needarr,0)
group2countlength2 = get_single_column_length(lengthgroup2,group2needarr,1)
group2countlength3 = get_single_column_length(lengthgroup2,group2needarr,2)
group2countlength4 = get_single_column_length(lengthgroup2,group2needarr,3)
counts2 = np.concatenate((group2countlength1,group2countlength2
                            ,group2countlength3,group2countlength4)).reshape(4,-1)
counts3 = np.ndarray([1,3,135]).astype('int')

group3countlength1 = get_single_column_length(lengthgroup3,group3needarr,0)
group3countlength2 = get_single_column_length(lengthgroup3,group3needarr,1)
group3countlength3 = get_single_column_length(lengthgroup3,group3needarr,2)
counts3= np.concatenate((group3countlength1,group3countlength2
                            ,group3countlength3)).reshape(3,-1)

wholecounts = np.zeros([1,11,135]).astype('int')


wholecounts =np.concatenate((mirnalength.reshape(1,-1),counts1,counts2,counts3),axis=0)


print ('whole array dim: ')
print (wholecounts.shape)



with h5py.File('02.matrix/new/length/testlength/'+filename.split('/')[-1]+'lengthbyseq') as f:
    f.create_dataset('counts',data = wholecounts)


'''

{
inds=$(seq 0 1 60)
for ind in $inds; do
echo bin/mappingmxlength_.py \
-i ${ind} 
done
}> Jobs/mappingwhole.txt
qsubgen -n mappingwhole -q Z-LU -a 1-60 -j 1 --bsub --task-file Jobs/mappingwhole.txt
bsub < Jobs/mappingwhole.sh

'''