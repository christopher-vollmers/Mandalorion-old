import numpy as np
import os
import sys


path=sys.argv[2]
content_file=sys.argv[1]
upstream_buffer=int(sys.argv[4])
downstream_buffer=int(sys.argv[3])


def collect_alt_events(path):
    alt_dict={}
    chromosome_list=set()

    for line in open(path+'/Sequential_SS.txt'):
      a=line.strip().split('\t')
      parts=a[0].split('-')

      chromosome=parts[-1]
      
      gene_name=a[0]
      try:
          bla=alt_dict[gene_name]
      except:
          alt_dict[gene_name]={}
          alt_dict[gene_name]['starts']=set() 
          alt_dict[gene_name]['ends']=set() 
          alt_dict[gene_name]['alt']={}
          alt_dict[gene_name]['ret']={}
          alt_dict[gene_name]['seq']=[]


      seqs=a[1].split(',')[:-1]
      alt_dict[gene_name]['seq'].append(seqs)




    for line in open(path+'/Retention.txt'):
      a=line.strip().split('\t')
      parts=a[0].split('-')

      chromosome=parts[-1]
      
      gene_name=a[0]
      try:
          bla=alt_dict[gene_name]
      except:
          alt_dict[gene_name]={}
          alt_dict[gene_name]['starts']=set() 
          alt_dict[gene_name]['ends']=set() 
          alt_dict[gene_name]['alt']={}
          alt_dict[gene_name]['ret']={}
          alt_dict[gene_name]['seq']=[]


      retention5=a[1]
      retention3=a[2]
      retention5_coord=(int(retention5.split('_')[1]),int(retention5.split('_')[2]))
      retention3_coord=(int(retention3.split('_')[1]),int(retention3.split('_')[2]))

      left_start=min(retention5_coord[1],retention3_coord[1])
      right_start=max(retention5_coord[0],retention3_coord[0])
      try:
          alt_dict[gene_name]['ret'][left_start].add((retention5_coord,retention3_coord,'retention'))
          alt_dict[gene_name]['ret'][left_start].add((retention5_coord,retention3_coord,'retention'))
      except:
          alt_dict[gene_name]['ret'][left_start]=set()
          alt_dict[gene_name]['ret'][left_start].add((retention5_coord,retention3_coord,'retention'))




    file_list=[path+'/Alt_5_SS.txt',path+'/Alt_3_SS.txt']
    for infile in file_list:
        for line in open(infile):
            a=line.strip().split('\t')
            chromosome=a[0].split('_')[-1]        
            gene_name=a[0]

            starts5=(int(a[1].split('_')[1]),int(a[1].split('_')[2]))
            starts3=[]
            for entry in a[2].split(',')[:-1]:
                starts3.append((int(entry.split('_')[1]),int(entry.split('_')[2])))
       

            for entry in starts3:
                try:
                    alt_dict[gene_name]['alt'][starts5].add((entry,'splice'))
                except:
                    try:
                        alt_dict[gene_name]['alt'][starts5]=set()
                        alt_dict[gene_name]['alt'][starts5].add((entry,'splice'))
                    except:
                        alt_dict[gene_name]={}
                        alt_dict[gene_name]['alt']={}
                        alt_dict[gene_name]['ret']={}
                        alt_dict[gene_name]['starts']=set() 
                        alt_dict[gene_name]['ends']=set() 
                        alt_dict[gene_name]['seq']=[]
                        alt_dict[gene_name]['alt'][starts5]=set()
                        alt_dict[gene_name]['alt'][starts5].add((entry,'splice'))



    Isoform_counter={}

    for line in open(path+'/Matched_Combined_TESS_SS.txt'):
        a=line.strip().split('\t')
        if len(a)>5:
            gene_name=a[0]
            chromosome=a[1]
            chromosome_list.add(chromosome)
            Isoform_counter[gene_name]=0
            Direction=a[3]
            TSS=a[4].split(',')[:-1]
            TES=a[5].split(',')[:-1]
            starts=TSS
            ends=TES

            if len(starts)>0 and len(ends)>0:
               for entry in starts:
                  position=(entry.split('_')[0][0],int(entry.split('_')[1])-upstream_buffer,int(entry.split('_')[2])+downstream_buffer)
                  try:
                      alt_dict[gene_name]['starts'].add(position) 

                  except:
                      try:
                          alt_dict[gene_name]['starts']=set() 
                          alt_dict[gene_name]['starts'].add(position) 
                      except:
                          alt_dict[gene_name]={}
                          alt_dict[gene_name]['alt']={}
                          alt_dict[gene_name]['ret']={}
                          alt_dict[gene_name]['seq']=[]
                          alt_dict[gene_name]['starts']=set() 
                          alt_dict[gene_name]['starts'].add(position) 



               for entry in ends:
                  position=(entry.split('_')[0][0],int(entry.split('_')[1])-downstream_buffer,int(entry.split('_')[2])+upstream_buffer)
                  try:
                      alt_dict[gene_name]['ends'].add(position) 

                  except:
                      try:
                          alt_dict[gene_name]['ends']=set() 
                          alt_dict[gene_name]['ends'].add(position) 
                      except:
                          alt_dict[gene_name]={}
                          alt_dict[gene_name]['alt']={}
                          alt_dict[gene_name]['ret']={}
                          alt_dict[gene_name]['seq']=[]
                          alt_dict[gene_name]['ends']=set() 
                          alt_dict[gene_name]['ends'].add(position) 

    return alt_dict,chromosome_list,Isoform_counter

def generate_map_and_detail_dicts(alt_dict,chromosome_list):
    map_dict={}
    for chromosome in chromosome_list:
        map_dict[chromosome]={}
    detail_dict={}
    gene_set=set()
    for gene,data in alt_dict.items():
        chromosome='_'.join(gene.split('_')[1:])
        starts=sorted(list(data['starts']),key=lambda x: x[1], reverse=False)
        ends=sorted(list(data['ends']),key=lambda x: x[1], reverse=True)
        alts=data['alt']
        rets=data['ret']
        seqs=data['seq']

        detail_dict[gene]={}
        detail_dict[gene]['ret']={}
        detail_dict[gene]['alt']={}
        detail_dict[gene]['seqs']={}
        start_and_end=0
        alt_starts=0
        alt_splicing=0

        if len(starts)>0 and len(ends)>0:
            start_and_end=1

        if len(starts)>1 or len(ends)>1:
            alt_starts=1
        if len(alts)>0 or len(rets)>0:
            alt_splicing=1
        if start_and_end==1:
#          if alt_splicing==1 or alt_starts==1:
            gene_set.add(gene)
            for start in starts: 

                for x in range(start[1],start[2],1):
                 try:
                    bla=map_dict[chromosome][x]
                 except:
                    map_dict[chromosome][x]={}
                 for end in ends:
                        for y in range(end[1],end[2],1):
                            map_dict[chromosome][x][y]=(gene,start[0]+'_'+str(start[1])+'_'+str(start[2]),end[0]+'_'+str(end[1])+'_'+str(end[2]))




            if len(alts)>0:
             for alt in alts:
              associated_junctions=alts[alt]
              for associated_junction in associated_junctions:
                position=associated_junction

                for x in range(alt[0],alt[1],1):
                    try:
                        bla=detail_dict[gene]['alt'][x]
                    except:
                        detail_dict[gene]['alt'][x]={} 
                    for y in range(position[0][0],position[0][1],1):
                        detail_dict[gene]['alt'][x][y]=('alt',alt[0],position[0][0])


            if len(rets)>0:
             for ret in rets:
              associated_junctions=rets[ret]

              for associated_junction in associated_junctions:
                first=associated_junction[0]
                second=associated_junction[1]
                mixed=sorted(first+second)
                type1=associated_junction[2]
                for x in range(mixed[1],mixed[2],1):
                     detail_dict[gene]['ret'][x]=('retention',tuple(mixed))
            group=0
            if len(seqs)>0:

               for seq in seqs:
                 group+=1
                 for seq1 in seq:
                  seq_list=[]
                  for seq2 in seq:
                      if seq2!=seq1:
                          seq_list.append(seq2)
                      for x in range(int(seq1.split('_')[1]),int(seq1.split('_')[2]),1):
                          detail_dict[gene]['seqs'][x]=(seq_list,group,seq1)



    return map_dict,detail_dict

                
alt_dict,chromosome_list,Isoform_counter=collect_alt_events(path)
map_dict,detail_dict=generate_map_and_detail_dicts(alt_dict,chromosome_list)




number=-1
Isoform_Ident={}
isoform_dict={}
read_numbers={}
isoform_set=set()
for line in open(content_file):

  bad=0

  b=line.strip().split('\t')
  infile=b[0]
  fasta_file=b[1]
  individual_path=b[2]     
  os.system('mkdir '+individual_path+'/parsed_reads')
  os.system('rm '+individual_path+'/parsed_reads/*')

  length=0
  for line in open(fasta_file):
      length+=1
  infile1=open(fasta_file,'r')
  read_numbers[infile]=length/2
  x=0
  read_dict={}
  while x<length:
      first=infile1.readline()
      second=infile1.readline()
      name=first.strip()[1:]
      read_dict[name]=[first.strip(),second]
      x+=2





  for line in open(infile):

    stop=0
    a=line.strip().split('\t')
    chromosome=a[13]
    name=a[9]
    begin=int(a[15])
    span=int(a[16])

    


    first_number=int(a[11])
    last_number=(int(a[10])-int(a[12]))
    if a[8]=='+':
        extra_start_bases=int(first_number)
        extra_end_bases=int(last_number)
        left_ispcr=name.split('_')[-2]
        right_ispcr=name.split('_')[-1]

    if a[8]=='-':
        extra_start_bases=int(last_number)
        extra_end_bases=int(first_number)
        left_ispcr=name.split('_')[-1]
        right_ispcr=name.split('_')[-2]

    matches=name.split('_')
    ratio=sum(np.array(a[18].split(',')[:-1],dtype=int))/int(a[10])

    matched=0
    try:
        gene_data=map_dict[chromosome][begin][span]
        
        gene_chromosome='_'.join(gene_data[0].split('_')[1:])
        matched=1 
    except:
        pass


    if matched==1:# and extra_start_bases<15 and extra_end_bases<15 and left_ispcr==right_ispcr=='p':

        if gene_chromosome==chromosome:
             
          start_info=gene_data[1]
          end_info=gene_data[2]
          start=int(start_info.split('_')[1])
          end=int(end_info.split('_')[2]) 

 
          gene_name=gene_data[0]



          left_bin=gene_data[1][0]
          right_bin=gene_data[2][0]
          if left_bin!=right_bin:
            gene_info=detail_dict[gene_name]

            blocksizes=a[18].split(',')[:-1]
            blockstarts=a[20].split(',')[:-1]
            readstarts=a[19].split(',')[:-1]
            bases=[]
            indels=[]
            covered_list=[]
            alt_splice_set=set()
            seq_splice_set=set()
            seq_splice_present=set()

            for x in range(0,len(blocksizes)-1,1):
                blockstart=int(blockstarts[x])
                blocksize=int(blocksizes[x])
                left_splice=blockstart+blocksize
                right_splice=int(blockstarts[x+1])
                if right_splice-left_splice>50:# and int(readstarts[x])+blocksize==int(readstarts[x+1]):

                  try:
                    match=detail_dict[gene_name]['alt'][left_splice][right_splice]
                    alt_splice_set.add(match)
                  except:
                    pass
                  try:
                    match=detail_dict[gene_name]['alt'][right_splice][left_splice]
                    alt_splice_set.add(match)
                  except:
                    pass

                  try:
                    matches=detail_dict[gene_name]['seqs'][left_splice]
                    for match in matches[0]:
                        seq_splice_set.add(int(match.split('_')[1]))
                        seq_splice_present.add(matches[1])
                  except:
                    pass
                  try:
                    matches=detail_dict[gene_name]['seqs'][right_splice]
                    for match in matches[0]:
                        seq_splice_set.add(int(match.split('_')[1]))
                        seq_splice_present.add(matches[1])
                  except:
                    pass


                for y in range(blockstart,blockstart+blocksize,1):

                    try:
                        match=detail_dict[gene_name]['ret'][y]
                        covered_list.append(match)
                    except:
                        pass


            

            retention_list=[]
            for retention in set(covered_list):

                if covered_list.count(retention)/(retention[1][2]-retention[1][1])>0.3:
                     retention_list.append(retention)


            identity=''
            identity+=gene_name+'_'+str(start_info)+'_'+str(end_info)+'_'


            retention_present=set()
            for entry in sorted(retention_list, key=lambda x: x[1]):
                identity+=str(entry[1][0])[-4:]+'-'+str(entry[1][2])[-4:]+','
                retention_present.add(entry[1][0])
                retention_present.add(entry[1][2])
            identity+='_'   

            splice_out={}
            alt_splices_present=set() 
            rights=set()
            for entry in sorted(list(alt_splice_set), key=lambda x: x[1]):

                identity+=str(entry[1])[-4:]+'-'+str(entry[2])[-4:]+','
                rights.add(entry[2])
                alt_splices_present.add((entry[1]))
                for intro_spliced in range(min(entry[1],entry[2])+1,max(entry[1],entry[2]),1):
                    splice_out[intro_spliced]=1

            alt_splices_required=set()
            for key in alt_dict[gene_name]['alt']:
                

                if start<key[0]<end:
                  for second in alt_dict[gene_name]['alt'][key]:

                    if start<int(second[0][0])<end:  
                      try: 
                        bla=splice_out[key[0]]
                      except:
                        alt_splices_required.add(key[0]) 


            seq_splice_required=set()
            seq_splice_required_dict={}
            for key in detail_dict[gene_name]['seqs']:
                group=detail_dict[gene_name]['seqs'][key][1]
                site=detail_dict[gene_name]['seqs'][key][2]
                if start<key<end:
                    try: 
                        bla=splice_out[key]
                    except:
                        try:
                            seq_splice_required_dict[group].add(site) 
                        except:
                            seq_splice_required_dict[group]=set()
                            seq_splice_required_dict[group].add(site)

            seq_retention_set=set()
            for key in retention_present:
                try:
                    group=detail_dict[gene_name]['seqs'][key][1]
                    site=detail_dict[gene_name]['seqs'][key][2]
                    seq_retention_set.add(group)
                except:
                    pass

            for group in seq_splice_required_dict:
                if len(seq_splice_required_dict[group])>1:
                    seq_splice_required.add(group)



#            if chromosome=='SIRV3' and start==1935 and end==8949 and 6323 in alt_splices_present:
#                print(start,end,alt_splices_required,alt_splices_present, retention_present,seq_splice_set,seq_splice_present,seq_splice_required)
#                print(begin,start,end,alt_splices_present,((alt_splices_required-retention_present)-seq_splice_set),retention_present,seq_retention_set,(seq_splice_required-seq_retention_set),seq_splice_present)
            if alt_splices_present==((alt_splices_required-retention_present)-seq_splice_set) and (seq_splice_required-seq_retention_set) <= seq_splice_present:



              try:
                  Isoform_Identity=Isoform_Ident[identity]
              except:
                  Isoform_counter[gene_name]+=1
                  Isoform_Ident[identity]='Isoform'+str(Isoform_counter[gene_name])
                  Isoform_Identity=Isoform_Ident[identity]

#              if chromosome=='SIRV3' and start==1935 and end==8949 and 6323 in alt_splices_present:

#                  print('success',Isoform_Ident[identity],identity)

              detailed_identity=identity
              filename=gene_name[:200]+'_'+Isoform_Identity+'_'+left_bin+'_'+right_bin

              identity=Isoform_Identity
              isoform_set.add(Isoform_Identity+'\t'+detailed_identity)

              try:
                isoform_dict[gene_name][identity][infile]+=1
              except:
                try:
                    isoform_dict[gene_name][identity][infile]=1
                except:
                    try:
                        isoform_dict[gene_name][identity]={}
                        isoform_dict[gene_name][identity][infile]=1
                    except:
                        isoform_dict[gene_name]={}
                        isoform_dict[gene_name][identity]={}
                        isoform_dict[gene_name][identity][infile]=1
              


              out_reads=open(individual_path+'/parsed_reads/'+filename,'a')
              out_reads.write('>'+read_dict[name][0][1:]+'_'+str(extra_start_bases)+'_'+str(extra_end_bases)+'\n'+read_dict[name][1])
              out_reads.close() 

out_iso=open(path+'/Isoforms.txt','w')
for isoform in isoform_set:
    out_iso.write(isoform+'\n')
out=open(path+'/Isoform_expression.txt','w')

for gene in sorted(list(isoform_dict.keys())):
    gene_total={}
    for line in open(content_file):
        b=line.strip().split('\t')
        infile=b[0]
        gene_total[infile]=0
    out.write(gene+'\t')

    for identity in sorted(list(isoform_dict[gene].keys())):
        for line in open(content_file):
            b=line.strip().split('\t')
            infile=b[0] 

            try: 
                expression=(isoform_dict[gene][identity][infile]/read_numbers[infile])*10000
            except: 
                expression=0
                isoform_dict[gene][identity][infile]=0

            gene_total[infile]+=expression

    for line in open(content_file):
        b=line.strip().split('\t')
        infile=b[0]
        out.write(str(gene_total[infile])+',')

    out.write('\t')

    for identity in sorted(list(isoform_dict[gene].keys())):
        out.write(identity+'_')
        for line in open(content_file):
            b=line.strip().split('\t')
            infile=b[0]
            try: 
                expression=(isoform_dict[gene][identity][infile]/read_numbers[infile])*10000
            except: 
                expression=0
            out.write(str(expression)+',')
        out.write('\t')
    out.write('\n')






















          
   




