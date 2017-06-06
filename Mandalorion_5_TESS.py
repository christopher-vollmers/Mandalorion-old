import sys
import numpy

content_file=sys.argv[1]  ### This file should simply contain the paths to to all the psl files you want to include in the analysis. The aligned reads have to be parsed using Isopore_3.
out_path=sys.argv[2]
cutoff=float(sys.argv[3])
genome_file=sys.argv[4]
refine=sys.argv[5]

TSS='ATGG'
TES='TTTT'

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

def add_to_coverage(begin,span,blocksizes,blockstarts,histo_coverage,chromosome):
    coverage_set=set()
    coverage_set.add(myround(begin))
    coverage_set.add(myround(span))

    for x in range(0,len(blocksizes),1):
        blockstart=int(blockstarts[x])
        blocksize=int(blocksizes[x])
        if blocksize>10:
            for y in range(0,blocksize,10):
                rounded=myround(blockstart+y)
                coverage_set.add(rounded)

            for yy in range(y,blocksize,1):
                rounded=myround(blockstart+yy)
                coverage_set.add(rounded)
    
    for rounded in coverage_set:           
        try:
            histo_coverage[chromosome][rounded]+=1
        except:
            histo_coverage[chromosome][rounded]=1
    return coverage_set,histo_coverage


def collect_reads(content_file):
    histo_left_bases={}
    histo_right_bases={}
    chromosome_list_left=set()
    chromosome_list_right=set()
    histo_coverage={}
    start_bases=[]
    end_bases=[]
    for line in open(content_file):
        print(line.strip())
        b=line.strip().split('\t')
        infile=b[0]
        read_seq=read_seq_file(b[1])
        for lines in open(infile):
            a=lines.strip().split('\t')
            
            name=a[9]

            direction=a[8]
            seq=read_seq[name]
            if direction=='+':
                start_seq=seq[2:6]
                end_seq=reverse_complement(seq)[2:6]
                pick=[-2,-1]
            if direction=='-':
                start_seq=reverse_complement(seq)[2:6]
                end_seq=seq[2:6]
                pick=[-1,-2] 

            matches=name.split('_')
            chromosome=a[13]
            try:
                bla=histo_coverage[chromosome]
            except:
                histo_coverage[chromosome]={}


            begin=int(a[15])
            span=int(a[16])

            if direction=='+':
                extra_start_bases=int(a[11])
                extra_end_bases=int(a[10])-int(a[12])

            if direction=='-':
                extra_end_bases=int(a[11])
                extra_start_bases=int(a[10])-int(a[12])

            start_bases.append(extra_start_bases)
            end_bases.append(extra_end_bases)
            blocksizes=a[18].split(',')[:-1]
            blockstarts=a[20].split(',')[:-1]

            coverage_set,histo_coverage=add_to_coverage(begin,span,blocksizes,blockstarts,histo_coverage,chromosome)

            
            if start_seq!=end_seq:

                low_bound=int(begin)
                up_bound=int(span)
                if matches[pick[0]]=='p':
                    chromosome_list_left.add(chromosome)
                    try:
                        histo_left_bases[chromosome][low_bound].append([extra_start_bases,begin,span,coverage_set,start_seq,matches[pick[0]]])
                    except:
                         try: 
                             histo_left_bases[chromosome][low_bound]=[]
                             histo_left_bases[chromosome][low_bound].append([extra_start_bases,begin,span,coverage_set,start_seq,matches[pick[0]]])
                         except:
                             histo_left_bases[chromosome]={}
                             histo_left_bases[chromosome][low_bound]=[]
                             histo_left_bases[chromosome][low_bound].append([extra_start_bases,begin,span,coverage_set,start_seq,matches[pick[0]]])

                if matches[pick[1]]=='p':
                    chromosome_list_right.add(chromosome)
                    try:
                        histo_right_bases[chromosome][up_bound].append([extra_end_bases,begin,span,coverage_set,end_seq,matches[pick[1]]])
                    except:
                        try: 
                            histo_right_bases[chromosome][up_bound]=[]
                            histo_right_bases[chromosome][up_bound].append([extra_end_bases,begin,span,coverage_set,end_seq,matches[pick[1]]])
                        except:
                            histo_right_bases[chromosome]={}
                            histo_right_bases[chromosome][up_bound]=[]
                            histo_right_bases[chromosome][up_bound].append([extra_end_bases,begin,span,coverage_set,end_seq,matches[pick[1]]])

    chromosome_list=chromosome_list_left & chromosome_list_right
    return histo_left_bases,histo_right_bases,chromosome_list,histo_coverage



def scan_for_best_bin(entry,density_dict,distance_range,iterator_shift,bottom,top):
    best_extra_list=[]
    best_extra_list_bases=[]
    base_match=[]
    best_adapter=[]
    peak_center=0
    coverage_area=[]
    for x in distance_range:
        extra_list_bases=[]
        extra_list_expression=[]
        extra_list_starts=[]
        extra_list_ends=[]
        extra_adapter=[]
        coverage_set=[]
        base_seqs=[]

        bla=0

        for y in distance_range:
            try:

               bla=peak_areas[chromosome][side][entry+x+y]  
            except:
               pass

        if bla==0:  
            for y in distance_range:

                try: 
                    for item in density_dict[entry+x+y]:
                        extra_list_bases.append(item[0])
                        extra_list_expression.append(1)
                        
                        base_seqs.append(item[4])
                        for covered_position in item[3]:
                            coverage_set.append(covered_position)
                except:
                        pass
        if bottom<=numpy.median(extra_list_bases)<top:
             if len(extra_list_expression)>len(best_extra_list):
                  best_extra_list=extra_list_expression
                  best_extra_list_bases=extra_list_bases

                  coverage_area=coverage_set
                  base_match=base_seqs
                  peak_center=entry+x  

    return best_extra_list,best_extra_list_bases,coverage_area,base_match,peak_center 


def block_out_peak_area(peak_center,peak_areas,iterator_shift,Peaks,coverage_area,reverse,chromosome,side):
    for base in range(peak_center-(10*iterator_shift),peak_center+(10*iterator_shift),iterator_shift):
        peak_areas[chromosome][side][base]=Peaks+1

    iterator=0
    count=0
    base_minus=peak_center
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
                base_minus=base_f
                for base in range(peak_center,base_minus+(iterator_shift*5),iterator_shift):
                    peak_areas[chromosome][side][base]=str(Peaks+1)+'_'+str(counter)
    return peak_areas


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
    coverage_area1=sorted(coverage_area3,reverse=reverse)     
    counter=0
    for base_f in coverage_area1:
         count=0 
         if reverse==False:
             if base_f>peak_center:
                    count=1
         elif reverse==True:
             if base_f<peak_center:
                    count=1

         if count==1:
             if counter<=4:
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


def reverse_complement(sequence):
  Seq=''
  complement = {'a':'T','c':'G','g':'C','t':'A','n':'N','A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
  for item in sequence[::-1]:
    Seq=Seq+complement[item]
  return Seq

def myround(x, base=10):
    return int(round(x,-1))


def find_peaks(density_dict,Peaks,reverse,bottom,top,bottom25,top75,cutoff,chromosome,histo_coverage,side,peak_areas):

    if reverse==False:
        distance_range=range(-10,10,1)
        iterator_shift=1
    if reverse==True:
        distance_range=range(10,-10,-1)
        iterator_shift=-1


    entry_list=[]
    for entry in density_dict:
        entry_list.append([entry,density_dict[entry]])
 
    for entry,density in sorted(entry_list,key=lambda x: int(x[0]),reverse=reverse):

        if len(density_dict[entry])>=2: 
            best_extra_list_bases=[]
            for item in density_dict[entry]:
                best_extra_list_bases.append(item[0])
            if bottom<=numpy.median(best_extra_list_bases)<=top: 
                 if bottom25<=numpy.percentile(best_extra_list_bases,25): 
                     if numpy.percentile(best_extra_list_bases,75)<=top75:
                         try:
                             bla=peak_areas[chromosome][side][entry]
                             if bla=='annotation':
                                 best_extra_list,best_extra_list_bases,coverage_area,base_match,peak_center=scan_for_best_bin(entry,density_dict,distance_range,iterator_shift,bottom,top)
                                 coverage,coverage_area=determine_coverage(coverage_area,chromosome,reverse,peak_center,histo_coverage)
                                 peak_areas=block_out_peak_area(peak_center,peak_areas,iterator_shift,Peaks,coverage_area,reverse,chromosome,side)

                         except:
                             best_extra_list,best_extra_list_bases,coverage_area,base_match,peak_center=scan_for_best_bin(entry,density_dict,distance_range,iterator_shift,bottom,top)
                             coverage,coverage_area=determine_coverage(coverage_area,chromosome,reverse,peak_center,histo_coverage)

                             if coverage>0:
                                if round(sum(best_extra_list)/coverage,3)>cutoff:
                                   TSS_distance=[]
                                   TES_distance=[]
                                   for seq in base_match:
            
                                       if seq==TSS:
                                           TSS_distance.append(0)
                                           TES_distance.append(4) 
                                       if seq==TES:
                                           TSS_distance.append(4)
                                           TES_distance.append(0) 
                                   matched=0
                                   if reverse==False:
                                       if numpy.median(TES_distance)>numpy.median(TSS_distance) and numpy.median(TSS_distance)<=2:
                                           out.write(chromosome+'\t'+str(peak_center-10)+'\t'+str(peak_center+10)+'\t'+'S'+side+str(Peaks)+'_'+str(peak_center-10)+'_'+str(peak_center+10)+'_'+str(sum(best_extra_list)/coverage)+'\n')
                                           matched=1

                                       if numpy.median(TSS_distance)>numpy.median(TES_distance) and numpy.median(TES_distance)<=2:
                                           out.write(chromosome+'\t'+str(peak_center-10)+'\t'+str(peak_center+10)+'\t'+'E'+side+str(Peaks)+'_'+str(peak_center-10)+'_'+str(peak_center+10)+'_'+str(sum(best_extra_list)/coverage)+'\n')
                                           matched=1 

                                   if reverse==True:
                                       if numpy.median(TES_distance)>numpy.median(TSS_distance) and numpy.median(TSS_distance)<=2:
                                           out.write(chromosome+'\t'+str(peak_center-10)+'\t'+str(peak_center+10)+'\t'+'S'+side+str(Peaks)+'_'+str(peak_center-10)+'_'+str(peak_center+10)+'_'+str(sum(best_extra_list)/coverage)+'\n')
                                           matched=1

                                       if numpy.median(TSS_distance)>numpy.median(TES_distance) and numpy.median(TES_distance)<=2:
                                           out.write(chromosome+'\t'+str(peak_center-10)+'\t'+str(peak_center+10)+'\t'+'E'+side+str(Peaks)+'_'+str(peak_center-10)+'_'+str(peak_center+10)+'_'+str(sum(best_extra_list)/coverage)+'\n')
                                           matched=1 

                                   if matched==1:
                                       Peaks+=1
                                       peak_areas=block_out_peak_area(peak_center,peak_areas,iterator_shift,Peaks,coverage_area,reverse,chromosome,side)

                    
                                        

    return Peaks,peak_areas



def parse_genome(input_file,left_bounds,right_bounds):

    gene_dict={} 
    for line in open(input_file):
        a=line.strip().split('\t')
        if len(a)>7:
             if a[2]=='exon':
                 try: 
                     gene_dict[a[8].split('; transcript_id "')[1].split('"')[0]].append((a[0],a[3],a[4],a[6]))
                 except:
                     gene_dict[a[8].split('; transcript_id "')[1].split('"')[0]]=[]
                     gene_dict[a[8].split('; transcript_id "')[1].split('"')[0]].append((a[0],a[3],a[4],a[6]))

    read_list=[]
    for transcript_id in gene_dict:
        transcript_data=gene_dict[transcript_id]

        chromosome=transcript_data[0][0]
        try:
            bla=right_bounds[chromosome]
        except:
            left_bounds[chromosome]={}
            right_bounds[chromosome]={}
            left_bounds[chromosome]['S']=[]
            right_bounds[chromosome]['S']=[]
            left_bounds[chromosome]['E']=[]
            right_bounds[chromosome]['E']=[]

        start=sorted(transcript_data,key=lambda x: int(x[1]))[0][1]
        end=sorted(transcript_data,key=lambda x: int(x[2]),reverse=True)[0][2]

        if transcript_data[0][3]=='+':
             left_bounds[chromosome]['S'].append(int(start)) 
             right_bounds[chromosome]['E'].append(int(end)) 
        elif transcript_data[0][3]=='-':
             left_bounds[chromosome]['E'].append(int(start))
             right_bounds[chromosome]['S'].append(int(end))  

    return left_bounds,right_bounds

def make_genome_bins(bounds,side,Peaks,chromosome,peak_areas,reverse):
    for type1 in ['S','E']:
        covered={}
        position_list=sorted(bounds[type1],key=int,reverse=reverse)
        for index1 in range(0,len(position_list),1):
            try:
                bla=peak_areas[chromosome][side][position_list[index1]]
            except:
                start=position_list[index1]-10
                end=position_list[index1]+10
                out.write(chromosome+'\t'+str(start)+'\t'+str(end)+'\t'+type1+side+str(Peaks)+'_'+str(start)+'_'+str(end)+'_A'+'\t'+str(Peaks)+'\n')
                if side =='l':
                    for base in range(start,end+40,1):                             
                        peak_areas[chromosome][side][base]='annotation'
                if side =='r':
                    for base in range(start-40,end,1):                             
                        peak_areas[chromosome][side][base]='annotation'

                Peaks+=1                    
                    
    return Peaks,peak_areas



   
histo_left_bases, histo_right_bases,chromosome_list,histo_coverage=collect_reads(content_file)
left_bounds={}
right_bounds={}
left_bounds,right_bounds=parse_genome(genome_file,left_bounds,right_bounds)


out=open(out_path+'TESS.bed','w')

Left_Peaks=0
Right_Peaks=0



peak_areas={}
for chromosome in chromosome_list:
    peak_areas[chromosome]={}
    peak_areas[chromosome]['l']={}
    peak_areas[chromosome]['r']={}
    print(chromosome)
    if 'g' in refine:
        Left_Peaks_old=Left_Peaks
        Right_Peaks_old=Right_Peaks
        Left_Peaks,peak_areas=make_genome_bins(left_bounds[chromosome],'l',Left_Peaks,chromosome,peak_areas,False)
        Right_Peaks,peak_areas=make_genome_bins(right_bounds[chromosome],'r',Right_Peaks,chromosome,peak_areas,True)
        print('Annotation-Based', Left_Peaks-Left_Peaks_old,Right_Peaks-Right_Peaks_old)
    Left_Peaks_old=Left_Peaks
    Right_Peaks_old=Right_Peaks
    Left_Peaks,peak_areas=find_peaks(histo_left_bases[chromosome],Left_Peaks,False,6,15,0,20,cutoff,chromosome,histo_coverage,'l',peak_areas)
    Right_Peaks,peak_areas=find_peaks(histo_right_bases[chromosome],Right_Peaks,True,6,15,0,20,cutoff,chromosome,histo_coverage,'r',peak_areas)
    print('Read-Based', Left_Peaks-Left_Peaks_old,Right_Peaks-Right_Peaks_old)



print(Left_Peaks)
print(Right_Peaks)












    
