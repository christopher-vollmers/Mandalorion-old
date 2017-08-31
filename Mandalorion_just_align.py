import sys
import os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-f','--fastq',type=str)
parser.add_argument('-g','--gmap_genome',type=str)
parser.add_argument('-G','--gmap_path',type=str)
parser.add_argument('-q','--quality_cutoff',type=str, default='9')
parser.add_argument('-t','--gmap_threads',type=str, default='1')
args=parser.parse_args()

fastq=args.fastq
gmap_genome=args.gmap_genome
gmap_path=os.path.abspath(args.gmap_path)
quality_cutoff=args.quality_cutoff
gmap_threads=args.gmap_threads

fastq=os.path.abspath(fastq)
fastq_name=fastq.split('/')[-1]
path='/'.join(fastq.split('/')[:-1])+'/'
print(path)

scripts=os.path.dirname(os.path.realpath(__file__))
print(scripts)

os.system('python3 %s/Mandalorion_1_Clean_Read_and_File_Names.py %s %s %s ' %(scripts,path, fastq_name,quality_cutoff))
os.system('python3 %s/Mandalorion_3_Remove_ISPCR_Sequences.py %s ' %(scripts,path))
os.system('python3 %s/Mandalorion_4_Align_Reads_With_Gmap.py %s %s %s %s' %(scripts,path, gmap_path,gmap_genome,gmap_threads))
content_file=open(path+'/content_file','w')
content_file.write('%s\t%s\t%s\n' %(path+'/2D_trimmed_l_gmapoutput_filtered.psl',path+'/2D_trimmed_l_filtered.fasta',path+'/'))



