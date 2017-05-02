import numpy as np
import sys
import os


def loop(input_dict):
    
  output_dict={}
  for entry1 in input_dict:
    try:
        bla=output_dict[entry1]
    except:
        output_dict[entry1]=set()
    length2=1
    length1=0 
    while length2>length1:    
        length1=len(output_dict[entry1]) 
        for entry2 in input_dict[entry1]:
            output_dict[entry1].add(entry2)
            try:
                output_dict[entry2].add(entry1)
            except:
                output_dict[entry2]=set()
                output_dict[entry2].add(entry1)
 
            for entry3 in input_dict[entry2]:
                output_dict[entry1].add(entry3)
                try:
                    output_dict[entry3].add(entry1)
                    output_dict[entry3].add(entry2)
                except:
                    output_dict[entry3]=set()
                    output_dict[entry3].add(entry1)
                    output_dict[entry3].add(entry2)



        length2=len(output_dict[entry1]) 


  return output_dict


feature_dict={}
feature_dict['Sl']={}
feature_dict['El']={}
feature_dict['5l']={}
feature_dict['3l']={}
feature_dict['Sr']={}
feature_dict['Er']={}
feature_dict['5r']={}
feature_dict['3r']={}

path=sys.argv[2]
content_file=sys.argv[1]


chromosome_list=set()
for line in open(path+'/SS.bed'):
    a=line.strip().split('\t')
    chromosome_list.add(a[0])

for line in open(path+'/TESS.bed'):
    a=line.strip().split('\t')
    chromosome_list.add(a[0])

for chromosome in chromosome_list:
    feature_dict['Sl'][chromosome]={}
    feature_dict['El'][chromosome]={}
    feature_dict['5l'][chromosome]={}
    feature_dict['3l'][chromosome]={}
    feature_dict['Sr'][chromosome]={}
    feature_dict['Er'][chromosome]={}
    feature_dict['5r'][chromosome]={}
    feature_dict['3r'][chromosome]={}


for line in open(path+'/SS.bed'):
    a=line.strip().split('\t')
    chromosome=a[0]
    start=int(a[1])
    end=int(a[2])
    type1=a[3]
    for position in np.arange(start,end,1):
        feature_dict[type1[:2]][chromosome][position]=type1

for line in open(path+'/TESS.bed'):
    a=line.strip().split('\t')
    chromosome=a[0]
    start=int(a[1])
    end=int(a[2])
    type1=a[3]
    for position in np.arange(start,end,1):
        feature_dict[type1[:2]][chromosome][position]=type1


chromosome_dict={}

master_dict={}
for line in open(content_file):
  print(line)
  bad=0
  b=line.strip().split('\t')
  infile=b[0]
  for line in open(infile):
    big_list=set()
    lefts=[]
    rights=[]
    stop=0
    a=line.strip().split('\t')
    chromosome=a[13]
    name=a[9]
    begin=int(a[15])
    span=int(a[16])
    try:
        big_list.add(feature_dict['Sl'][chromosome][begin])
    except:
        pass
    try:
        big_list.add(feature_dict['El'][chromosome][begin])
    except:
        pass
    try:
        big_list.add(feature_dict['Sr'][chromosome][span])
    except:
        pass
    try:
        big_list.add(feature_dict['Er'][chromosome][span])
    except:
        pass
                 


    blocksizes=a[18].split(',')[:-1]
    blockstarts=a[20].split(',')[:-1]
    readstarts=a[19].split(',')[:-1]
    bases=[]
    indels=[]
    covered_list=[]
    alt_splice_set=set()
    for x in range(0,len(blocksizes)-1,1):
        blockstart=int(blockstarts[x])
        blocksize=int(blocksizes[x])
        left_splice=blockstart+blocksize
        right_splice=int(blockstarts[x+1])
        if right_splice-left_splice>50:
            try:
                big_list.add(feature_dict['5l'][chromosome][left_splice])
            except:
                pass
            try:
                big_list.add(feature_dict['3l'][chromosome][left_splice])
            except:
                pass
            try:
                big_list.add(feature_dict['5r'][chromosome][right_splice])
            except:
                pass
            try:
                big_list.add(feature_dict['3r'][chromosome][right_splice])
            except:
                pass
            
    for entry1 in big_list:
        chromosome_dict[entry1]=chromosome
        try:
            for entry2 in big_list:
                master_dict[entry1].add(entry2)
        except:
            master_dict[entry1]=big_list

print(1)
master_dict2=loop(master_dict)
print(2)
master_dict3=loop(master_dict2)
print(3)
master_dict4=loop(master_dict3)
#master_dict5=loop(master_dict4)

master_set=set()
for entry in master_dict4:
    master_set.add(tuple(sorted(list(master_dict4[entry]))))   


out=open(path+'/Combined_TESS_SS.txt','w')
counter=0
for sets in master_set:
    counter+=1
    left_ends=[]
    right_ends=[]
    left_splice=[]
    right_splice=[]
    for entry in sets:
        if entry[1]=='l':
             if entry[0] in ['S','E']:
                 if entry[0]=='S':
                     direction='+'
                 else:
                     direction='-'
                 left_ends.append(entry)
             else:
                 left_splice.append(entry)
        if entry[1]=='r':
             if entry[0] in ['S','E']:
                 if entry[0]=='S':
                     direction='-'
                 else:
                     direction='+'
                 right_ends.append(entry)
             else:
                 right_splice.append(entry)
    if len(left_ends)>0 or len(right_ends)>0:
        out.write(str(counter)+'_'+chromosome_dict[entry]+'\t'+chromosome_dict[entry]+'\t-\t'+str(direction)+'\t')
        for entry in left_ends:
            out.write(entry+',')
        out.write('\t')  
        for entry in right_ends:
            out.write(entry+',')
        out.write('\t')
        for entry in left_splice:
            out.write(entry+',')
        out.write('\t')
        for entry in right_splice:
            out.write(entry+',')
        out.write('\n')    


















          
   




