from __future__ import division
import numpy
import numpy as np
import sys


def myround(x, base=10):
    return int(round(x,-1))

def load_gene_annotation_data(infile):
    all_splices={}
    combined_splices={}
    read_length_dict={}
    data_dict={}

    for line in open(infile):
     if line[0]!='#':
      a=line.strip().split('\t')
 

      if a[2]=='exon':
       transcript_name=a[8].split('transcript_id "')[1].split('"')[0]
       transcript=transcript_name
       chromosome=a[0]
       try:
           bla=all_splices[chromosome]
       except:
           all_splices[chromosome]={}
           combined_splices[chromosome]={}


       start=int(a[3])
       end=int(a[4])
       direction=a[6] 
       try:
           bla=data_dict[transcript]
       except:
           data_dict[transcript]={}
           data_dict[transcript]['starts']=[]
           data_dict[transcript]['ends']=[]

       data_dict[transcript]['chromosome']=chromosome
       data_dict[transcript]['direction']=direction
       data_dict[transcript]['starts'].append(start) 
       data_dict[transcript]['ends'].append(end)





    for transcript in data_dict:
      chromosome=data_dict[transcript]['chromosome']
      direction=data_dict[transcript]['direction']  
      starts=sorted(data_dict[transcript]['starts'],key=int)
      ends=sorted(data_dict[transcript]['ends'],key=int)

      
      combined=sorted(starts+ends,key=int)

      for x in range(0,len(combined)-1,1):
         left=combined[x]
         right=combined[x+1]
         if x%2==0:
             for y in range(left,right,10): 
               try:        
                    combined_splices[chromosome][myround(y)].add(transcript) 
               except:        
                    combined_splices[chromosome][myround(y)]=set() 
                    combined_splices[chromosome][myround(y)].add(transcript) 

             for yy in range(y,right,1):
                try:        
                    combined_splices[chromosome][myround(yy)].add(transcript) 
                except:        
                    combined_splices[chromosome][myround(yy)]=set() 
                    combined_splices[chromosome][myround(yy)].add(transcript) 

         else:
            if right-left>50:
                  try:
                     all_splices[chromosome][left].append(right)

                  except: 
                     all_splices[chromosome][left]=[]
                     all_splices[chromosome][left].append(right)
                  try:
                     all_splices[chromosome][right].append(left)
                  except:
                     all_splices[chromosome][right]=[]
                     all_splices[chromosome][right].append(left)



    return all_splices,combined_splices,read_length_dict


def count_matches_retention(Type1, ONT,cutoff):
  out=open(Type1,'w')
  Alt_Splice=0
  Alt_Splice_Match=0
  for key in ONT:
    if ONT[key]>0:
      key_split=key.split('~')
      out.write(key_split[2]+'\t'+key_split[0]+'\t'+key_split[1]+'\t'+str(ONT[key])+'\n')
      Alt_Splice+=1


  out.close() 
  print(Type1,Alt_Splice,Alt_Splice_Match)



def count_matches_alt_splicing(Type1,ONT):
  out=open(Type1,'w')
  Alt_Splice=0
  Alt_Splice_Match=0
  for key in ONT:
    total=0
    for second in ONT[key]:
        total+=ONT[key][second]
    print(key)
    first_splice=key.split('_')[-3]+'_'+key.split('_')[-2]+'_'+key.split('_')[-1]
    name=key.split('_')[0]+'_'+key.split('_')[1]
    ONT_2nd=set()
    for second in ONT[key]:
      if ONT[key][second]>1 and ONT[key][second]/total>0.01:

        ONT_2nd.add(second)

    if len(ONT_2nd)>1:
      out.write(name+'\t'+first_splice+'\t')
      for entry in ONT_2nd:
          out.write(entry+',')
      out.write(str(total)+'\n')
      Alt_Splice+=1

  print(Type1,Alt_Splice,Alt_Splice_Match)



def collect_splices(infile):
  total=0
  total1=0

  for line in open(infile):
    a=line.replace('\n','').split('\t')

    chromosome=a[1]
    gene_name=a[0]
    expression=np.array(a[2].split(',')[:-1],dtype=int)

    Direction=a[3]
    Lefts=a[4].split(',')[:-1]
    Rights=a[5].split(',')[:-1]
    try:
        bla=data_dict1[chromosome]
    except:
        data_dict1[chromosome]={}
        data_dict2[chromosome]={}


    Left_Splices=a[6].split(',')[:-1]
    Right_Splices=a[7].split(',')[:-1]

    splice_list=[]

    for Splice in Left_Splices:
        total+=1
    for Splice in Right_Splices:
        total+=1


    for Splice in Left_Splices:
        splice_list.append(Splice)
    for Splice in Right_Splices:
        splice_list.append(Splice)

    splice_list_sorted=sorted(splice_list, key=lambda x: int(x.split('_')[1]))


    for x in range(0,len(splice_list_sorted),1):
        type_string='' 
        total1+=1
        type_string+=splice_list_sorted[x][:2]
        type1=splice_list_sorted[x][:2]
        for y in range(0,len(splice_list_sorted),1):
         if x<y:
           if y<x+10:
            type_string+=splice_list_sorted[y][:2]
            type2=splice_list_sorted[y][:2]
            bin1=range(int(splice_list_sorted[x].split('_')[1])-1,int(splice_list_sorted[x].split('_')[2]),1) 
            bin2=range(int(splice_list_sorted[y].split('_')[1])-1,int(splice_list_sorted[y].split('_')[2]),1) 
            binbetween=range(myround(min(int(splice_list_sorted[x].split('_')[2]),int(splice_list_sorted[y].split('_')[2]))),myround(max(int(splice_list_sorted[x].split('_')[1]),int(splice_list_sorted[y].split('_')[1]))),10)
            match=0
            if type1=='5l' and type2=='3r':
              match=1
            if type1=='3l' and type2=='5r':
              match=1
            if match==1:

              try:
                    bla=data_dict1[chromosome][type_string]
              except:
                    data_dict1[chromosome][type_string]={}
                    data_dict2[chromosome][type_string]={}

              for base1 in bin1:
                try: 
                    bla=data_dict1[chromosome][type_string][base1]
                except: 
                    data_dict1[chromosome][type_string][base1]={}
                for base2 in bin2:
                    data_dict1[chromosome][type_string][base1][base2]=splice_list_sorted[x]+'~'+splice_list_sorted[y]+'~'+gene_name
              if y==x+1:
               data_dict2[chromosome][type_string][splice_list_sorted[x]+'~'+splice_list_sorted[y]+'~'+gene_name]={}
               for base1 in binbetween:
                data_dict2[chromosome][type_string][splice_list_sorted[x]+'~'+splice_list_sorted[y]+'~'+gene_name][base1]=1

           else:
              break 

  return data_dict1, data_dict2


def sort_reads_into_splices(data_dict,splice_dict):
 Splices_53={}
 Splices_35={}
 for chromosome in data_dict:
  if chromosome in splice_dict:
   for key in data_dict[chromosome]:
     data_test=data_dict[chromosome][key]
     for splices in splice_dict[chromosome]:
        try:
            test_bin=data_test[splices]

            for second_splices in splice_dict[chromosome][splices]:

                try:
                    info=test_bin[second_splices]

                    try: 
                        Splices_53[info.split('~')[2]+'_'+info.split('~')[0]][info.split('~')[1]]+=1 
                    except: 
                        try:
                            Splices_53[info.split('~')[2]+'_'+info.split('~')[0]][info.split('~')[1]]=1 
                        except:
                            Splices_53[info.split('~')[2]+'_'+info.split('~')[0]]={}
                            Splices_53[info.split('~')[2]+'_'+info.split('~')[0]][info.split('~')[1]]=1 

                    try: 
                        Splices_35[info.split('~')[2]+'_'+info.split('~')[1]][info.split('~')[0]]+=1 
                    except: 
                        try:
                            Splices_35[info.split('~')[2]+'_'+info.split('~')[1]][info.split('~')[0]]=1 
                        except:
                            Splices_35[info.split('~')[2]+'_'+info.split('~')[1]]={}
                            Splices_35[info.split('~')[2]+'_'+info.split('~')[1]][info.split('~')[0]]=1 


                except:
                    pass
        except:
           pass 
 return Splices_53,Splices_35

def sort_reads_into_introns(data_dict,splice_dict,read_length_dict,cat):
 Retentions={}
 Complete_Retentions={}
 for chromosome in data_dict:
  if chromosome in splice_dict:
   for key in data_dict[chromosome]:
    for intron in data_dict[chromosome][key]:

     data_test=data_dict[chromosome][key][intron]
     matched=0
     to_erase=[]
     for splices in splice_dict[chromosome]:
        try:
            bla=data_test[splices]
            for read in splice_dict[chromosome][splices]:
                try:
                    Retentions[intron].append(read)
                except: 
                    Retentions[intron]=[]
                    Retentions[intron].append(read) 
                matched=1
            to_erase.append(splices) 
        except:
            pass
     for splices in to_erase:
         splice_dict[chromosome].pop(splices, None)
   
     if matched==1:
       intron_length=np.absolute((int(intron.split('~')[1].split('_')[1])-int(intron.split('~')[0].split('_')[1])))/10  


       for splice in set(Retentions[intron]):
         coverage=Retentions[intron].count(splice)

         if cat==1:
          if coverage/intron_length>0.7:
           try:   
             Complete_Retentions[intron]+=1
           except:
             Complete_Retentions[intron]=1
         if cat==2:
          read_length=read_length_dict[splice]/10
          if coverage/read_length>0.7:
           try:   
             Complete_Retentions[intron]+=1
           except:
             Complete_Retentions[intron]=1


  
 return Retentions,Complete_Retentions



def read_data(content_file):
  all_splices={}
  combined_splices={}
  read_length_dict={}
  for file1 in open(content_file):
    b=file1.split('\t')[0]
    for line in open(b):
        a=line.strip().split('\t')
        chromosome=a[13]

        try:
           bla=all_splices[chromosome]
        except:
           all_splices[chromosome]={}
           combined_splices[chromosome]={}

        read=a[9]
        read_length=int(a[10])
        read_length_dict[read]=read_length
        begin=int(a[15])
        span=int(a[16])
        previous_blocksize=-1
        previous_start=-1
        previous_blockend=100000000000000
        blocksizes=numpy.array(a[18].split(',')[:-1],dtype=int)
        blockstarts=numpy.array(a[20].split(',')[:-1],dtype=int)
        readstarts=numpy.array(a[19].split(',')[:-1],dtype=int)
        aligned_bases=sum(numpy.array(blocksizes,dtype=int))

        for x in range(0,len(blocksizes),1):
            left=int(blockstarts[x])

            if blocksizes[x]>10:           
              for y in range(0,int(blocksizes[x]),10):
                try:        
                    combined_splices[chromosome][myround(left+y)].add(read) 
                except:        
                    combined_splices[chromosome][myround(left+y)]=set() 
                    combined_splices[chromosome][myround(left+y)].add(read) 

              for yy in range(y,int(blocksizes[x]),1):
                try:        
                    combined_splices[chromosome][myround(left+yy)].add(read) 
                except:        
                    combined_splices[chromosome][myround(left+yy)]=set() 
                    combined_splices[chromosome][myround(left+yy)].add(read) 




        intron=0
        if aligned_bases/read_length>0.85:


            for x in range(0,len(blocksizes),1):
                blockstart=int(blockstarts[x])
                blocksize=int(blocksizes[x])
                readstart=int(readstarts[x])
                blockend=blockstart+blocksize



                if previous_start==-1:
                    previous_start=blockstart
                    min_length=10
                else:
                    min_length=10

                if blockstart-previous_blockend>20:
                    previous_start=blockstart

                if blockend-previous_start>min_length:
                    if intron==1:
                        try:
                            all_splices[chromosome][remember_blockend].append(previous_start)

                        except: 
                            all_splices[chromosome][remember_blockend]=[]
                            all_splices[chromosome][remember_blockend].append(previous_start)
                        try:
                            all_splices[chromosome][previous_start].append(remember_blockend)
                        except:
                            all_splices[chromosome][previous_start]=[]
                            all_splices[chromosome][previous_start].append(remember_blockend)


                        intron=0
                    else:
                        try:
                            next_blockstart=int(blockstarts[x+1])
                            next_blocksize=int(blocksizes[x+1])
                            next_readstart=int(readstarts[x+1])

                            insert=next_blockstart-blockend
                            if insert>50:
                                indel1=next_readstart-(readstart+blocksize) 
                                if indel1<=0:   
                                    remember_blockend=blockend
                                    remember_indel1=indel1
                                    remember_start=previous_start
  
                                    intron=1
                                    previous_start=next_blockstart 
                                    blockend=next_blockstart  
                        except:
                            pass

                previous_blockend=blockend 




  return all_splices,combined_splices,read_length_dict





intron_dict={}
counter=0


conversion={}
conversion['+']={}
conversion['-']={}
conversion['+']['R']='3'
conversion['+']['L']='5'
conversion['-']['R']='5'
conversion['-']['L']='3'

data_dict1={}
data_dict2={}


content_file=sys.argv[1]
path=sys.argv[2]

data_dict1,data_dict2=collect_splices(path+'/Matched_Combined_TESS_SS.txt')


print(1)
all_ONT_splices, combined_ONT_splices,read_length_dict_ONT=read_data(content_file)
print(2)
Retention_ONT,Complete_Retention_ONT=sort_reads_into_introns(data_dict2,combined_ONT_splices,read_length_dict_ONT,1)
count_matches_retention(path+'Retention.txt',Complete_Retention_ONT,1)


Matches_ONT_53,Matches_ONT_35=sort_reads_into_splices(data_dict1,all_ONT_splices)
count_matches_alt_splicing(path+'Alt_3_SS.txt',Matches_ONT_53)
count_matches_alt_splicing(path+'Alt_5_SS.txt',Matches_ONT_35)


