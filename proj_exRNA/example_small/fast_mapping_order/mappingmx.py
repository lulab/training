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
parser.add_argument('-i', dest='ind')  #0~9
parser.add_argument('-o', dest='output_file')
#parser.add_argument('-r', dest='indind')
args = parser.parse_args()


originalnames = np.array(['lncRNA', 'miRNA', 'mRNA', 'piRNA', 'snoRNA', 'snRNA', 'srpRNA',
       'tRNA', 'vaultRNA', 'Y_RNA', 'tucpRNA'])
rnanames = np.array([ 'miRNA', 
            'Y_RNA', 'piRNA', 'srpRNA', 
            'snoRNA', 'snRNA', 'tRNA','vaultRNA',
            'tucpRNA', 'lncRNA','mRNA' ])
xsorted = np.argsort(originalnames)
ypos = np.searchsorted(originalnames[xsorted], rnanames)
transind  = xsorted[ypos]


filename = np.loadtxt('filename.txt',dtype='str')[int(args.ind)]  #'17402567-B.all.mx'
samplearrmx = np.array(pd.read_table(filename)
                       .iloc[:,1:])[:,transind]
#/Share/home/caojingyi/exRNA/process/16.one_map/02.matrix
print ('sample matrix loaded')
print (samplearrmx.shape)

mirnacountafter = np.where(samplearrmx[:,0]==0)  #filter after mirna
mirnacount = samplearrmx.shape[0] - mirnacountafter[0].shape[0]

group1needarr = samplearrmx[mirnacountafter[0]][:,1:4]
group1countafter = np.where(np.sum(group1needarr,axis=1)==0)
group1needind = np.where(np.sum(group1needarr,axis=1)!=0)
group1needarr = group1needarr[group1needind]

group2needarr = samplearrmx[mirnacountafter[0]][group1countafter][:,4:8]
group2countafter = np.where(np.sum(group2needarr,axis=1)==0)
group2needind = np.where(np.sum(group2needarr,axis=1)!=0)
group2needarr = group2needarr[group2needind]

group3needarr = samplearrmx[mirnacountafter[0]][group1countafter][group2countafter][:,8:11]
group3needind = np.where(np.sum(group3needarr,axis=1)!=0)
group3needarr = group3needarr[group3needind]

print (group1needarr.shape,group2needarr.shape,group3needarr.shape)

counts1 = np.ndarray([1,3]).astype('int')
uniind1 = np.unique(( group1needarr !=0).argmax(axis=1),return_counts=True)
counts1[0][uniind1[0]]= uniind1[1]

counts2 = np.ndarray([1,4]).astype('int')
uniind2 = np.unique(( group2needarr !=0).argmax(axis=1),return_counts=True)
counts2[0][uniind2[0]]= uniind2[1]

counts3 = np.ndarray([1,3]).astype('int')
uniind3 = np.unique(( group3needarr !=0).argmax(axis=1),return_counts=True)
counts3[0][uniind3[0]]= uniind3[1]       
    
wholecounts = np.zeros([1,10]).astype('int')

wholecounts =np.concatenate((counts1,counts2,counts3),axis=1)
wholecounts = np.concatenate((np.repeat(mirnacount,1).reshape(-1,1),wholecounts),axis=1)


with h5py.File('02.matrix/new/counts/5.31/'+filename.split('/')[-1]+'countsbyseq') as f:
    f.create_dataset('counts',data = wholecounts)


'''
{
inds=$(seq 0 1 60)
for ind in $inds; do
echo bin/mappingmx_.py \
-i ${ind} 
done
}> Jobs/mappingwhole.txt
qsubgen -n mappingwhole -q Z-LU -a 1-60 -j 1 --bsub --task-file Jobs/mappingwhole.txt
bsub < Jobs/mappingwhole.sh

#bpeek -J "mappingwhole[2]"
'''