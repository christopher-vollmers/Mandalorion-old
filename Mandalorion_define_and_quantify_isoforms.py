import sys
import os
import argparse


parser=argparse.ArgumentParser()
parser.add_argument('-c','--content_file',type=str)
parser.add_argument('-p','--path',type=str)
parser.add_argument('-a','--genome_annotation',type=str)
parser.add_argument('-g','--gmap_genome',type=str)
parser.add_argument('-l','--gene_list',type=str)
parser.add_argument('-i','--illumina_content_file',type=str, default='-')
parser.add_argument('-r','--refine',type=str,default=None,choices=['g','gi','i'])  # if set to 'i' or 'gi' and no illumina_content_file is provided, the i will be ignored
args=parser.parse_args()


content_file=args.content_file         # file containing paths to psl and fasta/q files you want analyzed 
path=args.path+'/'	         #path where you want your output files to go
genome_annotation=args.genome_annotation    #path to your gtf file		
gmap_genome=args.gmap_genome	         #name of the gmap genome. Just whatever you name it when you use gmap-build to make it
gene_list=args.gene_list            #list of gene name for which you want to create consensus   
illumina_content_file=args.illumina_content_file  #If you have Illumina_data it will be used to improve Splice Sites
refine=args.refine

os.system('python3 Mandalorion_5_TESS.py %s %s %s %s %s' %(content_file, path,'0.02', genome_annotation, refine)) # This script identifies TSS and TES. Output file 'TESS.bed' is written to 'path'. The last argument determines how stringent the script filters potential TSS/TES. The lower the more stringent. 0.03 is our default
os.system('python3 Mandalorion_6_SS.py %s %s %s %s %s' %(content_file, path,'0.05', genome_annotation, refine)) # This script identifier splice sites. Output file 'SS.bed' is written to 'path'.
if 'i' in refine:
    os.system('python3 Mandalorion_6.5_Refine_SS.py %s %s %s %s' %(illumina_content_file,genome_annotation, path,refine)) # This script refines splice sites using genome annotation and/or Illumina data. Output file 'SS_refined.bed' is written to 'path'.
else:
    os.system('cp %s %s' %(path+'SS_raw.bed',path+'SS.bed')) # If SSs are not refined SS.bed is simply copied to SS_refined to keep the workflow consistent

os.system('python3 Mandalorion_7_Gene_Expression.py %s %s %s' %(content_file,path,genome_annotation)) # This script determine gene-level expression. Requires >8G of RAM for mouse/human genome. All other outputs are independent of this script so you can just comment it out if you are not interested in gene level data.
os.system('python3 Mandalorion_8_Combined_TESS_SS.py %s %s ' %(content_file,path)) # This script groups TSS/TES and SS if they are linked by raw reads
os.system('python3 Mandalorion_9_Match_TESS_SS_Combinations_to_genes.py %s %s ' %(genome_annotation,path)) # This script matches the groups wit annotated genes
os.system('python3 Mandalorion_10_Determine_Alternative_Splicing.py %s %s' %(content_file,path)) # This script looks at raw read alignments to determine alt. splicing and intron retentions
os.system('python3 Mandalorion_11_Define_and_Quantify_Isoforms.py %s %s %s %s' %(content_file, path,40,0)) # This script sort raw reads into isoform bins. The two number variables determine the window around TSS and TES in which read ends can fall and still be matched to the site.

for line in open(content_file):
    subpath=line.strip().split('\t')[2]
    os.system('python3 Mandalorion_12_Create_Consensi.py %s %s %s %s %s' %(subpath, gene_list,'0.009',100,50))  # This script uses poaV2 to create isoform consensus reads for each gene specified in the gene list. The number at the end is the ratio an isoform has to represent of all isoforms of each respective gene to be processed 
    os.system('python3 Mandalorion_14_Align_Consensi_With_Gmap.py %s %s' %(subpath,gmap_genome))  # This script aligns those isoform consensus reads to the genome.


