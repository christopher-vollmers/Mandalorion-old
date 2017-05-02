import sys
import os
import editdistance

path=sys.argv[1]
gene_file=sys.argv[2]
poa = 'poa'
score_matrix = 'NUC.4.4.mat'
subsample=100  # Change this if you want to subsample more or less reads per isoform
progressive_read_cutoff=50 ## Change this if you want to use the -do_progressive option for poa for isoform with more or less than 50 reads.

def read_fasta(infile):
    reads={}
    sequence=''
    
    for line in open(infile):
        a=line.strip()
        if a[0]=='>':
            if sequence!='':

                reads[name]=sequence
            name=a[1:].split()[0]
            sequence=''
        else:
            sequence+=a
    reads[name]=sequence
    return reads

def reverse_complement(sequence):
  Seq=''
  complement = {'A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
  for item in sequence[::-1]:
    Seq=Seq+complement[item]
  return Seq


def collect_file_paths(path,gene_file):
   genes_of_interest=[]
   for line in open(gene_file):
       genes_of_interest.append(line.strip())

   isoform_list=[]
   for file1 in sorted(os.listdir(path+'/parsed_reads')):
       for gene in genes_of_interest: 
           if gene in file1:
               file2=file1+'_sub'
               out_sub=open(path+'/parsed_reads/'+file2,'w') 
               counter=0
               isoform_reads=read_fasta(path+'/parsed_reads/'+file1)
               isoform_read_list=list(isoform_reads.keys())
               read1 = isoform_read_list[0]
               out_sub.write('>'+read1+'\n'+isoform_reads[read1]+'\n')
               for read2 in isoform_read_list[1::]:
                   if counter<subsample:
                       out_sub.write('>'+read2+'\n')
                       dist_1 = editdistance.eval(isoform_reads[read1],isoform_reads[read2])**2/float(len(isoform_reads[read1])*len(isoform_reads[read2]))
                       dist_2 = editdistance.eval(isoform_reads[read1],reverse_complement(isoform_reads[read2]))**2/float(len(isoform_reads[read1])*len(isoform_reads[read2]))
                       if dist_1 < dist_2:
                           out_sub.write(isoform_reads[read2]+'\n')
                       else:
                           out_sub.write(reverse_complement(isoform_reads[read2])+'\n')
                   counter+=1


               isoform_list.append(path+'/parsed_reads/'+file2)

   return isoform_list

def run_poa(isoform_list,out,cutoff):
    for isoform in isoform_list:
        print(isoform)
        FASTA = isoform
        number_of_reads=0
        for line in open(FASTA):
             number_of_reads+=1
        number_of_reads=int(number_of_reads/2)

        PIR = FASTA+'.pir'
        if number_of_reads<cutoff:
            os.system('%s -read_fasta %s -hb -pir %s -best -do_progressive %s' %(poa,FASTA,PIR,score_matrix))
        else:
            os.system('%s -read_fasta %s -hb -pir %s -best %s' %(poa,FASTA,PIR,score_matrix))        
        i=0
        for line in open(PIR):
            i+=1
        print(i)
        reads=read_fasta(PIR)
        print(len(reads))
        os.system('rm '+FASTA)
        for read in reads:
            if read == 'CONSENS0':
                 
                 out.write('>'+FASTA.split('/')[-1].split('_sub')[0]+'\n'+reads[read].replace('-','').replace('.','')+'\n')



combined_consensus_file=open(path+'/Isoform_Consensi.fasta','w')
isoform_list=collect_file_paths(path,gene_file)
run_poa(isoform_list,combined_consensus_file,progressive_read_cutoff)

    
