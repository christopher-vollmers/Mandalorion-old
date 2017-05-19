### This script demultiplexes a fastq file. You will need BLAT in your PATH to run it or add a path to your BLAT binary.
### The script will create a subfolder in your path for each of your samples

import sys
import os
import numpy as np

### This script requires three inputs

path=sys.argv[1]  ### path to where your fastq file is located
adapter_fasta=sys.argv[2]  ### fasta file containing adapter sequences
adapter_comb=sys.argv[3]  ### file containin sample names and index combinations

#os.system('blat -noHead -stepSize=1 -minScore=20 -minIdentity=20 '+adapter_fasta +' '+ path+'/2D.fasta ' +path+'/2D_adapter_alignment.psl')

readlist=[]
adapter_dict={}

length=0
for line in open(path+'2D.fastq'):
    length+=1


iterator=0
infile=open(path+'2D.fastq','r')
while iterator<length:
    line=infile.readline()
    sequence=infile.readline()
    plus=infile.readline()
    quality=infile.readline()
    name=line[1:].strip()
    readlist.append(name)
    adapter_dict[name]={}

    adapter_dict[name]['+']=[]
    adapter_dict[name]['-']=[]

    adapter_dict[name]['+'].append('-')
    adapter_dict[name]['-'].append('-')

    adapter_dict[name]['+_pos']=[]
    adapter_dict[name]['-_pos']=[]

    adapter_dict[name]['+_pos'].append(0)
    adapter_dict[name]['-_pos'].append(len(sequence))
    iterator+=4   


burn_dict={}
for line in open(path+'2D_adapter_alignment.psl'):
    a=line.strip().split('\t')
    if len(a)>8:
        read_name=a[9]
        adapter=a[13]
        strand=a[8]
        position='-'
        if strand=='+':
            if int(a[12])<300:
                position=int(a[12])+(int(a[14])-int(a[16]))
            else:
                burn_dict[read_name]=1

        if strand=='-':
            if int(a[11])>int(a[10])-300:
                position=int(a[11])-(int(a[14])-int(a[16]))
            else:
                burn_dict[read_name]=1

        if position!='-':
            if int(a[0])>20:
                adapter_dict[read_name][strand].append(adapter)
                adapter_dict[read_name][strand+'_pos'].append(position)

print('Reads with Internal Adapters', len(burn_dict))

adapters=[]
for line in open(adapter_fasta):
    a=line.strip()
    if a[0]=='>':

        adapters.append(a[1:])
 

count_dict={}
count_dict['-']={}
count_dict['-']['-']=0
for adapter1 in adapters:
    count_dict[adapter1]={}
    count_dict[adapter1]['-']=0
    count_dict['-'][adapter1]=0
    for adapter2 in adapters:
        count_dict[adapter1][adapter2]=0



for read in readlist:
    try:
        bla=burn_dict[read]
    except:
        if len(adapter_dict[read]['+'])<=2 or len(adapter_dict[read]['-'])<=2:
            count_dict[adapter_dict[read]['+'][-1]][adapter_dict[read]['-'][-1]]+=1



out_dict={}
count_dict={}
for line in open(adapter_comb):
    a=line.strip().split('\t')
    Sample_Name=a[0]
    index1=a[1]
    index2=a[2]
    print(index1,index2)
    out_dict[index1+'_'+index2]=path+'/'+Sample_Name+'/2D.fast'
    out_dict[index1+'_'+index1]=path+'/'+Sample_Name+'/2D.fast'
    out_dict[index2+'_'+index2]=path+'/'+Sample_Name+'/2D.fast'
    out_dict[index2+'_'+index1]=path+'/'+Sample_Name+'/2D.fast'
    out_dict[index1+'_-']=path+'/'+Sample_Name+'/2D.fast'
    out_dict['-_'+index2]=path+'/'+Sample_Name+'/2D.fast'
    out_dict[index2+'_-']=path+'/'+Sample_Name+'/2D.fast'
    out_dict['-_'+index1]=path+'/'+Sample_Name+'/2D.fast'
    count_dict[path+'/'+Sample_Name+'/2D.fast']=0

    os.system('mkdir '+path+'/'+Sample_Name)
    out=open(path+'/'+Sample_Name+'/2D.fasta','w')
    out.close()
    out=open(path+'/'+Sample_Name+'/2D.fastq','w')
    out.close()




total=0
good=0

x=0
infile=open(path+'2D.fastq')
while x<length:
  total+=1
  a=infile.readline().strip()
  b=infile.readline().strip()
  c=infile.readline().strip()
  d=infile.readline().strip()
  name=a[1:].strip()
  try:
    bla=burn_dict[name]
  except:
    try:
      if len(adapter_dict[name]['+'])<=2 and len(adapter_dict[name]['-'])<=2:
        outfile=out_dict[adapter_dict[name]['+'][-1]+'_'+adapter_dict[name]['-'][-1]]
        add_left='p'
        try:  
            position_plus= adapter_dict[name]['+_pos'][1]
        except:
            add_left='n'
            position_plus=0

        add_right='p'
        try:  
            position_minus= adapter_dict[name]['-_pos'][1]
        except:
            add_right='n'
            position_minus=len(b)-1

        outa=open(outfile+'a','a')
        outq=open(outfile+'q','a')
        if len(b[position_plus:position_minus])>0:
            count_dict[outfile]+=1
            good+=1
            outa.write('>'+a[1:]+'_'+add_left+'_'+add_right+'\n'+b[position_plus:position_minus]+'\n')
            outq.write(a+'_'+add_left+'_'+add_right+'\n'+b[position_plus:position_minus]+'\n'+c+'\n'+d[position_plus:position_minus]+'\n')

    except:
      pass
  x+=4

for outfile in count_dict:
    print(outfile,count_dict[outfile])

print(total, good)

