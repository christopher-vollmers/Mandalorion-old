import sys
import os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-s','--sample_sheet',type=str)
parser.add_argument('-f','--fastq',type=str)
parser.add_argument('-g','--gmap_genome',type=str)
parser.add_argument('-a','--adapter_fasta',type=str)
parser.add_argument('-q','--quality_cutoff',type=str)
args=parser.parse_args()

fastq=args.fastq
sample_sheet=args.sample_sheet
gmap_genome=args.gmap_genome
adapter_fasta=args.adapter_fasta
quality_cutoff=args.quality_cutoff

fastq_name=fastq.split('/')[-1]
path='/'.join(fastq.split('/')[:-1])+'/'
print(path)

content_file=open(path+'/content_file','w')
os.system('python3 Mandalorion_1_Clean_Read_and_File_Names.py %s %s %s ' %(path, fastq_name,quality_cutoff))
os.system('python3 Mandalorion_2_Demultiplex.py %s %s %s ' %(path,adapter_fasta,sample_sheet))
for line in open(sample_sheet):
    a=line.strip().split('\t')
    sub_path=path+a[0]
    os.system('python3 Mandalorion_3_Remove_ISPCR_Sequences.py %s ' %(sub_path))
    os.system('python3 Mandalorion_4_Align_Reads_With_Gmap.py %s %s ' %(sub_path, gmap_genome))
    content_file.write('%s\t%s\t%s\n' %(sub_path+'/2D_trimmed_l_gmapoutput_filtered.psl',sub_path+'/2D_trimmed_l_filtered.fasta',sub_path+'/'))


