import sys



genome_annotation=sys.argv[1]
path=sys.argv[2]
out=open(path+'/Matched_Combined_TESS_SS.txt','w')


def parse_genome_annotation(genome_annotation):
    gene_data={}
    for line in open(genome_annotation):
        a=line.strip().split('\t')
        if len(a)>6:
            if a[2]=='exon':
                start=int(a[3])
                end=int(a[4])
                chromosome=a[0]
                gene_name=a[8].split('gene_id "')[1].split('"')[0]

                try:
                    gene_start=gene_data[gene_name][0][1]
                    gene_end=gene_data[gene_name][0][2]
                    gene_data[gene_name][0]=(chromosome,min(start,gene_start),max(end,gene_end))
                except:
                    gene_data[gene_name]=(chromosome,start,end)


    return gene_data


def collect_gene_locations(gene_data):
    gene_locations={}
    for gene_name,gene_info in gene_data.items():
        gene_chromosome=gene_info[0]
        gene_start=gene_info[1]
        gene_end=gene_info[2]
        for location in range(round(gene_start,-2),round(gene_end+51,-2),50):
            try:
                gene_locations[gene_chromosome][location].append(gene_name)
            except:
                try:
                    gene_locations[gene_chromosome][location]=[]
                    gene_locations[gene_chromosome][location].append(gene_name)

                except:
                    gene_locations[gene_chromosome]={}
                    gene_locations[gene_chromosome][location]=[]
                    gene_locations[gene_chromosome][location].append(gene_name)

    return gene_locations


def match_feature_groups_to_genes(gene_locations,path,out):
    for line in open(path+'/Combined_TESS_SS.txt'):
        a=line.strip().split('\t')
        chromosome=a[1]
        edges=[]
        if len(a)>4:
            if a[4]!='':
                for item in a[4].split(',')[:-1]:
                    edges.append(int(item.split('_')[1]))
        if len(a)>5:
            if a[5]!='':
                for item in a[5].split(',')[:-1]:
                    edges.append(int(item.split('_')[1]))
        start=round(min(edges),-2)
        end=round(max(edges)+50,-2)
        matched_genes=[]
        for location in range(start,end,50):
            try:
                for gene in gene_locations[chromosome][location]:
                    matched_genes.append(gene)
            except:
                pass

        for gene in set(matched_genes):
            out.write(gene+'-')
        out.write(line)


gene_data=parse_genome_annotation(genome_annotation)
gene_locations=collect_gene_locations(gene_data)
match_feature_groups_to_genes(gene_locations,path,out)


