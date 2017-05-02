import sys


content_file=sys.argv[1]
path=sys.argv[2]
genome_annotation=sys.argv[3]
out=open(path+'/Gene_Expression.txt','w')


def read_data(file1):
    print(file1)
    data={}
    counter=0
    for line in open(file1):
        counter+=1
        a=line.split('\t')
        if len(a)>6:
            begin=int(a[15])
            span=int(a[16])

            read=a[9]
            chromosome=a[13]

            blocksizes=a[18].split(',')[:-1]
            blockstarts=a[20].split(',')[:-1]

            for x in range(0,len(blocksizes),1):
                blockstart=int(blockstarts[x])
                blocksize=int(blocksizes[x])
                for y in range(0,blocksize,1):
                    rounded=int((blockstart+y)/5)*5
                    try:
                                                  
                        data[chromosome][rounded].add(read+'_'+blocksizes[0])
                    except:
                        try:
                            data[chromosome][rounded]=set()
                            data[chromosome][rounded].add(read+'_'+blocksizes[0])
                        except:
                            data[chromosome]={} 
                            data[chromosome][rounded]=set()
                            data[chromosome][rounded].add(read+'_'+blocksizes[0])
                            


    print(counter)
    return data,counter


def parse_genome_annotation(genome_annotation):
    gene_data={}
    for line in open(genome_annotation):
        a=line.strip().split('\t')
        if len(a)>6:
            if a[2]=='exon':
                start=int(a[3])
                end=int(a[4])
                chromosome=a[0]
                gene_name=a[8].split('gene_id "')[1].split('"')[0]+'_'+chromosome

                try:
                    gene_start=gene_data[gene_name][0][1]
                    gene_end=gene_data[gene_name][0][2]
                    gene_data[gene_name][0]=(chromosome,min(start,gene_start),max(end,gene_end))
                except:
                    gene_data[gene_name]=[(chromosome,start,end),set()]

                for x in range(start,end,1):
                    gene_data[gene_name][1].add(int(x/5)*5)

    return gene_data



def quantify_gene_expression(gene_data,content_file):
    data_list=[]
    quant_dict={}
    for line in open(content_file):
        a=line.strip().split('\t')
        data_list.append(a[0])

    for x_count in range(0,len(data_list),1):

        data_dict={}
        data_dict[x_count],counter=read_data(data_list[x_count])
        for gene in gene_data:


            list1=gene_data[gene]
            gene_chromosome=list1[0][0]
            gene_start=list1[0][1]
            gene_end=list1[0][2]
            exon_length=len(list1[1])
            exon_bases=list1[1]


            matched=set()
            for exon_base in exon_bases:
                try: 
                    matched_reads=data_dict[x_count][gene_chromosome][exon_base]
                    for match in matched_reads:
                        matched.add(match)


                except:
                    pass
 
            try:
                bla=quant_dict[gene]
                quant_dict[gene][x_count]=(len(matched)/counter)*10000
            except:
                quant_dict[gene]={}
                quant_dict[gene][x_count]=(len(matched)/counter)*10000

    return quant_dict,data_list

def write_output_file(gene_data,quant_dict,out,data_list):
    for gene in gene_data:
        list1=gene_data[gene]
        gene_chromosome=list1[0][0]
        gene_start=list1[0][1]
        gene_end=list1[0][2]
        exon_length=len(list1[1])
        exon_bases=list1[1]


        out.write(gene+'\t'+gene_chromosome+'\t'+str(gene_start)+'\t'+str(gene_end)+'\t'+str(exon_length)+'\t')
        for x_count in range(0,len(data_list),1):
            out.write(str(round(quant_dict[gene][x_count],2))+'\t')
        out.write('\n')




gene_data=parse_genome_annotation(genome_annotation)
quant_dict,data_list=quantify_gene_expression(gene_data,content_file)
write_output_file(gene_data,quant_dict,out,data_list)
  
                



        

