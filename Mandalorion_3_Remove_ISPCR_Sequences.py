### This script removes ISPCR in 2D reads after they are already demultiplexed.
### It assumes your 2D reads are in a file called 2D.fastq and will return two files called 2D_trimmed_l.fastq 2D_trimmed_l.fasta
### Read names will be appended with 'p's if adapters are found and removed on either end, 
### i.e. _n_p would mean that an adapter was found on the right but not the left end of the read
### You'll need to install the "editdistance" module for python



import editdistance
import sys

path=sys.argv[1]
input1='2D.fastq'
output1='2D_trimmed_l'

ISPCR='AAGCAGTGGTATCAACGCAGAGTAC'

Counter_dict={}
Counter_dict['ATGG']=0
Counter_dict['TTTT']=0
def reverse_complement(sequence):
  Seq=''
  complement = {'a':'T','c':'G','g':'C','t':'A','n':'N','A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
  for item in sequence[::-1]:
    Seq=Seq+complement[item]
  return Seq

def find_ispcr(infile):
  outq=open(path+'/'+output1+'.fastq','w')
  outa=open(path+'/'+output1+'.fasta','w')
  
  length=0
  for line in open(infile):
    length+=1



  x=0

  file1=open(infile,'r')

  total=0
  forward=0
  reverse=0
  total_double=0
  difference=0
  while x<=length:
    a=file1.readline()
    b=file1.readline().strip()  
    c=file1.readline()
    d=file1.readline().strip()   
    reverse_b=reverse_complement(b)
    length_seq=len(b)
    forward_trimm='-'
    reverse_trimm='-'
    total+=1
    mindist1=12
    for xx in range(0,100,1):
        forward_seq=b[xx:xx+25]        
        dist=editdistance.eval(forward_seq,ISPCR)
        if dist<mindist1:
            for y in range((xx+25)-2,xx+35,1):
                if b[y:y+4]=='ATGG' or b[y:y+4]=='TTTT':

                    match_left=b[y:y+4]
                    mindist1=dist
                    forward_trimm=y-2
                    break

    mindist2=12
    for xx in range(0,100,1):
        reverse_seq=reverse_b[xx:xx+25]        
        dist=editdistance.eval(reverse_seq,ISPCR)
        if dist<mindist2:
            for y in range((xx+25)-2,xx+35,1):
                if reverse_b[y:y+4]=='ATGG' or reverse_b[y:y+4]=='TTTT':
                    Counter_dict[reverse_b[y:y+4]]+=1
                    match_right=reverse_b[y:y+4]
                    mindist2=dist
                    reverse_trimm=length_seq-(y-2)
                    break 

    x+=4

    add_left='n'
    add_right='n'
    if forward_trimm!='-':
        forward+=1
        add_left='p' 
    else:
        forward_trimm=0
    if reverse_trimm!='-':
        reverse+=1
        add_right='p'
    else:
        reverse_trimm=length_seq
    if add_right!='n' and add_left!='n':
        total_double+=1
        if match_left!=match_right:
            difference+=1

    outq.write(a.strip()+'_'+add_left+'_'+add_right+'\n'+b[forward_trimm:reverse_trimm]+'\n'+c+d[forward_trimm:reverse_trimm]+'\n')
    outa.write('>'+a.strip()[1:]+'_'+add_left+'_'+add_right+'\n'+b[forward_trimm:reverse_trimm]+'\n')
  print(total_double,total, (total_double/total)*100,(difference/total_double)*100, forward,reverse, (forward/total)*100, (reverse/total)*100,difference, (difference/total_double)*100) 
  print(Counter_dict)



find_ispcr(path+'/'+input1)


















