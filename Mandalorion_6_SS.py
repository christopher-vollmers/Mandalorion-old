import numpy
import sys




def scan_for_best_bin(entry,distance_range,iterator_shift,density_dict,extra_base_cutoff_bottom,extra_base_cutoff_top,peak_areas):
    best_extra_list=[]
    peak_center=0
    bases=[]
    coverage_area=[]
    best_direction_l=[] 
    best_direction_r=[]


    for x in distance_range:
        extra_list_bases=[]
        extra_list_expression=[]
        extra_list_starts=[]
        extra_list_ends=[]
        direction_l={}
        direction_r={} 
        coverage_set=[]
        bla=0
        for y in distance_range:
            try:
               bla=peak_areas[entry+x+y]  
            except:
               pass

        if bla==0:  
             highest_y=0
             highest_y_pos=0
             for y in distance_range:
                try: 
 

                    for item in density_dict[entry+x+y]:
                        
                        extra_list_bases.append(item[0])
                        extra_list_expression.append(1)
                        extra_list_starts.append(item[1])
                        extra_list_ends.append(item[2])
                        try:
                            direction_l[item[4]]+=1
                        except:
                            direction_l[item[4]]=1 
                        try:
                            direction_r[item[5]]+=1
                        except:
                            direction_r[item[5]]=1
                        for covered_position in item[3]:

                            coverage_set.append(covered_position)
                except:
                        pass

        if extra_base_cutoff_bottom<=numpy.median(extra_list_bases)<=extra_base_cutoff_top:
            if sum(extra_list_expression)>sum(best_extra_list):
                best_extra_list=extra_list_expression
                peak_center=entry+x

                bases=extra_list_bases
                coverage_area=coverage_set
                best_direction_l=direction_l 
                best_direction_r=direction_r


    return best_extra_list,peak_center,bases,coverage_area,best_direction_l,best_direction_r

def determine_coverage(coverage_area,chromosome,reverse,peak_center,histo_coverage):
    coverage=[]
    coverage.append(0)
    forward=0
    reverse1=0
    forward_missed=0
    reverse_missed=0
    coverage_area2=[]
    coverage_area3=[]
    for covered_position in set(coverage_area):
         coverage_area3.append(covered_position)
         if coverage_area.count(covered_position)>1:
               coverage_area2.append(covered_position)

 
    coverage_area=sorted(coverage_area2,reverse=reverse)
    
    counter=0
    for base_f in coverage_area:
         count=0 
         if reverse==False:
             if base_f>peak_center:
                    count=1
         elif reverse==True:
             if base_f<peak_center:
                    count=1

         if count==1:
             if counter<=3:
                  counter+=1
                  base_f=myround(base_f)
                  try:
                      bla=histo_coverage[chromosome][base_f]

                      coverage.append(histo_coverage[chromosome][base_f])
                  except:
                      pass
             else:
                   break
    coverage=max(coverage)
    return coverage, coverage_area


def read_seq_file(seq_file):

    read_seq={}
        
    length=0
    for line2 in open(seq_file):
        length+=1
    seq_file_open=open(seq_file,'r')
    counter=0
    while counter<length:
        fasta_name=seq_file_open.readline().strip()
        fasta_seq=seq_file_open.readline().strip()
        fasta_name=fasta_name[1:] 
        read_seq[fasta_name]=fasta_seq
        counter+=2
    return read_seq


def reverse_complement(sequence):
  Seq=''
  complement = {'a':'T','c':'G','g':'C','t':'A','n':'N','A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
  for item in sequence[::-1]:
    Seq=Seq+complement[item]
  return Seq



def myround(x, base=10):
    return int(base * round(float(x)/base))


def find_peaks(density_dict,out,Peaks,reverse,cutoff,extra_base_cutoff_bottom,extra_base_cutoff_top,histo_coverage,side):

    if reverse==False:
        distance_range=range(-10,10,1)
        iterator_shift=1
    if reverse==True:
        distance_range=range(10,-10,-1)
        iterator_shift=-1

    peak_areas={}
    entry_list=[]
    for entry in density_dict:
      entry_list.append([entry,density_dict[entry]])

    for entry,density in sorted(entry_list,key=lambda x: sum(numpy.array(x[1])[:,2]),reverse=True):
        if len(density)>=2:

      
            try:
                bla=peak_areas[entry]

            except:
                best_extra_list,peak_center,bases,coverage_area,best_direction_l,best_direction_r=scan_for_best_bin(entry,distance_range,iterator_shift,density_dict,extra_base_cutoff_bottom,extra_base_cutoff_top,peak_areas)

                coverage,coverage_area=determine_coverage(coverage_area,chromosome,reverse,peak_center,histo_coverage)


                if coverage>0:
                    if sum(best_extra_list)/coverage>cutoff:
                         try:
                             Left_TSS=best_direction_l['TSS']
                         except:
                             Left_TSS=0
                         try:
                             Left_TES=best_direction_l['TES']
                         except:
                             Left_TES=0
                         try:
                             Right_TSS=best_direction_r['TSS']
                         except:
                             Right_TSS=0
                         try:
                             Right_TES=best_direction_r['TES']
                         except:
                             Right_TES=0
                         Left_to_Right=Left_TSS+Right_TES
                         Right_to_Left=Left_TES+Right_TSS
                         Type='-'
          
                         if Left_to_Right>Right_to_Left:
                             if reverse==True:
                                 Type='5'
                             elif reverse==False:
                                 Type='3'
                         if Left_to_Right<Right_to_Left:
                             if reverse==True:
                                 Type='3'
                             elif reverse==False:
                                 Type='5'

                         if Type!='-':
                             Peaks+=1

                             out.write(chromosome+'\t'+str(peak_center-10)+'\t'+str(peak_center+10)+'\t'+str(Type)+side+str(Peaks)+'_'+str(peak_center-10)+'_'+str(peak_center+10)+'\t'+str(Peaks)+'\n')
                             


                             for base in range(peak_center-10,peak_center+10,1):                             
                                 peak_areas[base]=1

        else:
            break                    

    return Peaks





def collect_reads(content_file):
    histo_left_bases={}
    histo_right_bases={}
    chromosome_list_left=set()
    chromosome_list_right=set()
    histo_coverage={}

    extra_base_cutoff_top=5
    extra_base_cutoff_bottom=0
    for line in open(content_file):
        print(line.strip())
        read_seq={}
        total=0
        b=line.strip().split('\t')
        infile=b[0]
        seq_file=b[1]
        length=0
        read_seq=read_seq_file(b[1])
        for line in open(infile):
            total+=1
            a=line.strip().split('\t')
            chromosome=a[13]
            try:
                bla=histo_coverage[chromosome]
            except:
                histo_coverage[chromosome]={}


            score=int(a[0])
            direction=a[8]
            name=a[9]
            begin=int(a[15])
            span=int(a[16])
            seq=read_seq[name]
            length=len(seq)

            blocksizes=a[18].split(',')[:-1]
            blockstarts=a[20].split(',')[:-1]
            readstarts=a[19].split(',')[:-1]

            if direction=='+':
                start_seq=seq[2:6]
                end_seq=reverse_complement(seq)[2:6]

            if direction=='-':
                start_seq=reverse_complement(seq)[2:6]
                end_seq=seq[2:6]

            left_match='-'
            right_match='-'
            if start_seq=='ATGG':
                left_match='TSS'
            elif start_seq=='TTTT':
                left_match='TES'
            if end_seq=='ATGG':
                right_match='TSS'
            elif end_seq=='TTTT':
                right_match='TES'


            coverage_set=set()
            previous_blocksize=-1
            previous_start=-1
            previous_blockend=100000000000000
            intron=0
            indel=0
            low_bounds=[]
            up_bounds=[]
            indel1=0

            aligned_bases=0

            for x in range(0,len(blocksizes),1):
                blockstart=int(blockstarts[x])
                blocksize=int(blocksizes[x])
                readstart=int(readstarts[x])
                aligned_bases+=blocksize
                blockend=blockstart+blocksize
                if blocksize>10:


                   for y in range(0,blocksize,10):

                       rounded=myround(blockstart+y)
                       coverage_set.add(rounded)

                   for yy in range(y,blocksize,1):
                       rounded=myround(blockstart+yy)
                       coverage_set.add(rounded)
    





                   if previous_start==-1:
                       previous_start=blockstart
                       min_length=10
                   else:
                       min_length=10
 
                   if blockstart-previous_blockend>20:
                       previous_start=blockstart
 
                   if blockend-previous_start>min_length:
                       if intron==1:
                          up_bounds.append([previous_start,indel1,blockend])
                          low_bounds.append([remember_blockend,remember_indel1,remember_start])
                          intron=0
                       else:
                          try:
                             next_blockstart=int(blockstarts[x+1])
                             next_blocksize=int(blocksizes[x+1])
                             next_readstart=int(readstarts[x+1])

                             insert=next_blockstart-blockend
                             if insert>50:
                                 indel1=next_readstart-(readstart+blocksize) 
                                 remember_blockend=blockend
                                 remember_indel1=indel1
                                 remember_start=previous_start
                                 intron=1
                                 previous_start=next_blockstart 
                                 blockend=next_blockstart  
                          except:
                              pass

                   previous_blockend=blockend 

            for rounded in coverage_set:           
                try:
                    histo_coverage[chromosome][rounded]+=1
                except:
                    histo_coverage[chromosome][rounded]=1

            ratio=aligned_bases/length
            if ratio>0.90:
  
                for low_bound,indel1,blockend in low_bounds:  

                    chromosome_list_left.add(chromosome)
                    try:
                        histo_left_bases[chromosome][low_bound].append([indel1,begin,span,coverage_set,left_match,right_match])
                    except:
                        try: 
                            histo_left_bases[chromosome][low_bound]=[]
                            histo_left_bases[chromosome][low_bound].append([indel1,begin,span,coverage_set,left_match,right_match])
                        except:
                            histo_left_bases[chromosome]={}
                            histo_left_bases[chromosome][low_bound]=[]
                            histo_left_bases[chromosome][low_bound].append([indel1,begin,span,coverage_set,left_match,right_match])


                for up_bound,indel1,blockend in up_bounds:   

                    chromosome_list_right.add(chromosome)            
                    try:
                         histo_right_bases[chromosome][up_bound].append([indel1,begin,span,coverage_set,left_match,right_match])
                    except:
                         try: 
                             histo_right_bases[chromosome][up_bound]=[]
                             histo_right_bases[chromosome][up_bound].append([indel1,begin,span,coverage_set,left_match,right_match])
                         except:
                             histo_right_bases[chromosome]={}
                             histo_right_bases[chromosome][up_bound]=[]
                             histo_right_bases[chromosome][up_bound].append([indel1,begin,span,coverage_set,left_match,right_match])

    chromosome_list=chromosome_list_left & chromosome_list_right
    return histo_left_bases, histo_right_bases,chromosome_list,histo_coverage




content_file=sys.argv[1]  ### This file should simply contain the paths to to all the psl files you want to include in the analysis. The aligned reads have to be parsed using Isopore_3.
out_path=sys.argv[2]
cutoff=float(sys.argv[3])
histo_left_bases, histo_right_bases,chromosome_list,histo_coverage=collect_reads(content_file)



Left_Peaks=0
Right_Peaks=0

out=open(out_path+'/SS.bed','w')

for chromosome in chromosome_list:
    Left_Peaks=find_peaks(histo_left_bases[chromosome],out,Left_Peaks,True,cutoff,0,5,histo_coverage,'l')
    Right_Peaks=find_peaks(histo_right_bases[chromosome],out,Right_Peaks,False,cutoff,0,5,histo_coverage,'r')


print(Left_Peaks)
print(Right_Peaks)












    
