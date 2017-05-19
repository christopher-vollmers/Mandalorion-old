import numpy
import sys

illumina_content_file=sys.argv[1]
genome_annotation=sys.argv[2]
out_path=sys.argv[3]
mode=sys.argv[4]




def parse_genome(input_file,left_bounds,right_bounds):
    gene_dict={} 

    for line in open(input_file):
        a=line.strip().split('\t')
        if len(a)>7:
             if a[2]=='exon':
                 try: 
                     gene_dict[a[8].split('; transcript_id "')[1].split('"')[0]].append((a[0],a[3],a[4]))
                 except:
                     gene_dict[a[8].split('; transcript_id "')[1].split('"')[0]]=[]
                     gene_dict[a[8].split('; transcript_id "')[1].split('"')[0]].append((a[0],a[3],a[4]))

    read_list=[]
    for transcript_id in gene_dict:
        transcript_data=gene_dict[transcript_id]

        chromosome=transcript_data[0][0]
        try:
            bla=right_bounds[chromosome]
        except:
            left_bounds[chromosome]=[]
            right_bounds[chromosome]=[]
        start=sorted(transcript_data,key=lambda x: int(x[1]))[0][1]
        end=sorted(transcript_data,key=lambda x: int(x[2]),reverse=True)[0][2]
        print(start,end)
        for entry in transcript_data:
              if entry[1]!=start:
                  right_bounds[chromosome].append(int(entry[1])-1)  
                  right_bounds[chromosome].append(int(entry[1])-1)
              if entry[2]!=end:
                  left_bounds[chromosome].append(int(entry[2]))  
                  left_bounds[chromosome].append(int(entry[2]))
     
    return left_bounds,right_bounds
    




def collect_Illumina_read_alignments(illumina_content_file,left_bounds,right_bounds):

    for line in open(illumina_content_file):
        print(line.strip())
        total=0
        b=line.strip().split('\t')
        infile=b[0]
     
        for line in open(infile):
            total+=1
            a=line.strip().split('\t')
            chromosome=a[13]
            try:
                bla=right_bounds[chromosome]
            except:
                left_bounds[chromosome]=[]
                right_bounds[chromosome]=[]
 
            blocksizes=a[18].split(',')[:-1]
            blockstarts=a[20].split(',')[:-1]

            for x in range(0,len(blocksizes),1):
                blockstart=int(blockstarts[x])
                blocksize=int(blocksizes[x])

                blockend=blocksize+blockstart 
                try: 
                    nextblockstart=int(blockstarts[x+1])
                    right_bounds[chromosome].append(nextblockstart)
                    left_bounds[chromosome].append(blockend)
                except: 
                    pass

    return left_bounds,right_bounds





def load_splice_sites(infile):
    splice_dict={}
    match_dict={}
    for line in open(infile):
        a=line.strip().split('\t')
        chromosome=a[0]
        try:
            bla=splice_dict[chromosome]
        except:
            splice_dict[chromosome]={}
            splice_dict[chromosome]['l']={}
            splice_dict[chromosome]['r']={}

        start=int(a[1])
        end=int(a[2])
        name=a[3]
        splice_type=name[1]
        number=a[4]
        match_dict[line]={}
        for base in range(start,end,1):  
            splice_dict[chromosome][splice_type][base]=line

    return splice_dict,match_dict

def match_illumina_to_splice_sites(splice_dict,left_bounds,right_bounds,match_dict):
    
    for chromosome in left_bounds:
        for splice in left_bounds[chromosome]:
            try:
                line=splice_dict[chromosome]['l'][splice]
                
                try: 
                    match_dict[line][splice]+=1
                except:
                    match_dict[line][splice]=1
            except:
                pass
        for splice in right_bounds[chromosome]:
            try:
                line=splice_dict[chromosome]['r'][splice]

                try:
                    match_dict[line][splice]+=1
                except:
                    match_dict[line][splice]=1
            except:
                pass
    return match_dict
            
def split_splice_sites(infile,match_dict):
    counter=0
    for line in open(infile):
        split=0
        a=line.strip().split('\t')
        matched_splices=match_dict[line]
        if len(matched_splices)>1:
            
            splice_list=[]
            for splice in matched_splices:
                if matched_splices[splice]>1:
                    splice_list.append(splice)
            if len(splice_list)>1:

                sorted_splice_list=sorted(splice_list,key=int)
                splice_distances=[]
                for splice_pos in range(0,len(sorted_splice_list)-1,1):
                    splice_distances.append(int(sorted_splice_list[splice_pos+1])-int(sorted_splice_list[splice_pos]))
                if min(splice_distances)>4:

                  split=1
                  counter+=1

            
                  for x in range(0,len(sorted_splice_list),1):
                    if x!=0:
                        start=int(sorted_splice_list[x]-((sorted_splice_list[x]-sorted_splice_list[x-1])/2))
                        s=1
                    else:
                        start=int(sorted_splice_list[x])-10
                        s=0
                    if x!=len(sorted_splice_list)-1:
                        end=int(sorted_splice_list[x]+((sorted_splice_list[x+1]-sorted_splice_list[x])/2)) 
                        e=1
                    else:
                        end=int(sorted_splice_list[x])+10
                        e=0


#                    print(a[0]+'\t'+str(start)+'\t'+str(end)+'\t'+a[3].split('_')[0]+'.'+str(x+1)+'_'+str(start)+'_'+str(end)+'_I'+'\t'+a[4]+'\n')


                    out.write(a[0]+'\t'+str(start)+'\t'+str(end)+'\t'+a[3].split('_')[0]+'.'+str(x+1)+'_'+str(start)+'_'+str(end)+'_I'+'\t'+a[4]+'\n')
        if split==0:
             out.write(line)        
    print(counter)        
                       

out=open(out_path+'/SS.bed','w')
   
left_bounds={}
right_bounds={}

splice_dict,match_dict=load_splice_sites(out_path+'/SS_raw.bed')

if 'i' in mode:
    if illumina_content_file!='-':
        left_bounds,right_bounds=collect_Illumina_read_alignments(illumina_content_file,left_bounds,right_bounds)

match_dict=match_illumina_to_splice_sites(splice_dict,left_bounds,right_bounds,match_dict)
split_splice_sites(out_path+'/SS_raw.bed',match_dict)











    
