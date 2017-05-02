### This script will take your 2D fastq file and clean the read names. It will only keep everything before the first white space.
### It will also output the reads into 2 new files named 2D.fastq and 2D.fasta in the respective formats.

import sys
import numpy as np

#This script takes 3 arguments:

path=sys.argv[1]   ### Point at the folder your 2D fastq file is located in
file_name=sys.argv[2]   ### Name of 2D fastq file
quality_cutoff=int(sys.argv[3])  ### We usually go for '9'
    
infile=path+'/'+file_name
length=0
for line in open(infile):
    length+=1

out1=open(path+'/2D.fastq','w')
out2=open(path+'/2D.fasta','w')


in1=open(infile,'r')

x=4

readlength=dict()

while x<length:     

    a=in1.readline()
    b=in1.readline()
    c=in1.readline()
    d=in1.readline()
    
    a_new=a.strip().split()[0]+'\n'
    b_new=b.strip()+'\n'
    c_new=c.strip()+'\n'
    d_new=d.strip()+'\n'

    a_fasta='>'+a_new[1:]

    quals=[]
    for character in d.strip():
        number = ord(character)-33
        
        quals.append(number) 
    if np.average(quals)>=quality_cutoff:
        out1.write(a_new+b_new+c_new+d_new)
        out2.write(a_fasta+b_new)

    x+=4

















