# Changes

-Mandalorion_3_Remove_ISPCR_Sequences.py
    Removed bug that appended empty read at the end of output fastq/a files that would stall gmap.
-Mandalorion_demultiplex_and_align.py
    Added '-t' flag option to run multiple gmap threads. Default is 1
-Mandalorion_just_align.py
    Added '-t' flag option to run multiple gmap threads. Default is 1
-Mandalorion_define_and_quantify_isoforms.py
    Added '-r' flag (options 'g' or 'gi', or 'i') to use use genome annation when defining TESS and SS.
    Added command to run Mandalorion_6.5_Refine_SS.py
-Mandalorion_5_TESS.py
    Modified to accomodate '-r' flag. If '-r' flag contains 'g' it uses TESS defined in the genome annotation file in addition to ONT read data
-Mandalorion_6_SS.py
    Modified to accomodate '-r' flag. If '-r' flag contains 'g' it uses SS defined in the genome annotation file in addition to ONT read data
-Mandalorion_6.5_Refine_SS.py
    New script that is run if -r flag contains 'i'. 
    Looks for SS bins that contains more than 1 Illumina splice junction. 
    It will split those bins if splice junctions are more than 4 bp apart

-Mandalorion_12_Create_Consensi.py
    Added a minimum isoform ratio for which consensus sequences are generated. It's set to 0.01 by default. 
    Also exposed values for subsampling (default: 50) and progressive cutoff (default:20) for poa. Increase numbers for more accuracy but longer run time.
    To change these default values open Mandalorion_define_and_quantify_isoforms.py and look for the third, fourth, and fifth input options to this script.

-Mandalorion_2_Demultiplex.py 
    Fixed bug that only allowed demultiplexing of reads with Q values > 9.
    Exclude reads with adapter recognized in the middle of the read (internal)


# Mandalorion
Analysis Pipeline to analysis Nanopore RNAseq data

### Mandalorion v0.1 ###

This collection of scripts is intended to provide a complete workflow to process Oxford Nanopore reads to identify and quantify isoforms. 
Starting with basecalled (Metrichor or Albacore) 1D or 2D reads (better), 
Mandalorion outputs isoforms, their consensus sequences, and approximate expression levels.  
Mandalorion is still highly experimental so definitely take a good look at all the outputs to make sure they make sense.
Most of the code was written by Christopher Vollmers (UCSC) with help from Charles Cole (UCSC) for the wrapper script to poaV2 as well as quality control.
Mandalorion relies on gmap,blat, and poaV2 all of which come with their own citation/licenses.
So far, we tested the code on mouse single cell RNAseq (Biorxiv).

### Requirements ###

Mandalorion consists of several independent python scripts and 3 wrapper scripts to call them in the right way and order.  

### Hardware: ###

Mandalorion was tested on both Linux Mint, Ubuntu, and Mac but I’m sure you’ll run into creative problems on your own machine. 
The amount of RAM you will need to run this pipeline will depend on the size of your genome of interest. 
16GB of RAM should be sufficient for mouse/human transcriptomes.

### Software: ### 

You will need:

-python3

with the following packages installed:

-numpy
-scipy
-editdistance

Further, you will need

-gmap (http://research-pub.gene.com/gmap/)
-blat (https://genome.ucsc.edu/FAQ/FAQblat.html)

installed to your PATH. And you will need to build or download the gmap genome index of interest. 
If you download, make sure you put them wherever gmap looks for them by default.

Finally, to build consensus isoforms you will need 

-poaV2 (https://github.com/tanghaibao/bio-pipeline)

with the poa executable in your PATH as well.

The NUC.4.4.mat distance matrix used by poaV2 is provided.

### Usage: ### 

All Mandalorion scripts have to be in the same folder and you should be in that folder when you call the wrapper scripts.
For the first part of the Mandalorion pipeline, you can call wrapper scripts to either 'demultiplex and align' reads or 'just align' reads.


    This part of the pipeline uses blat to demultiplex and gmap to align the reads to the genome.

    Demultiplex and Align:
        The output will be a separate folder in the path specified under 2 for each sample containing:
            2D_trimmed_l_gmapoutput_filtered.psl    -> contains filtered alignments to the genome specified in 4 below
            2D_trimmed_l_filtered.fasta    -> contain all fasta entry present in the psl file

        This will process a single MinION run that contains multiple indexed samples.

        The wrapper script for this part of the pipeline is called like this:

        $python3 Mandalorion_demultiplex_and_align.py [Options...]

        Options must include

        -s, --sample_sheet	path/to/sample_sheet (example provided) 
        -f, --fastq	path/to/fastq_file (make sure it is not called 2D.fastq, or it will be overwritten)
        -g, --gmap_genome	gmap genome name (e.g. mm10 or hg38, depending on what you named it when you built it) 
        -a, --adapter_fasta	path/to/sample_index.fasta (Provided)
        -q, --quality_cutoff	sets the minimum average q-score for a read to be kept (We use '9')
        -t, --gmap_threads	sets the number of threads to run gmap with, default is 1 

    Just Align:
        The output will be 2 files in the path specified under 2 below:
            2D_trimmed_l_gmapoutput_filtered.psl    -> contains filtered alignments to the gene specified in 3 below
            2D_trimmed_l_filtered.fasta    -> contain all fasta entry present in the psl file

        This assumes that all reads in your fastq file come from the same sample.

        The wrapper script is called like this:

        $python3 Mandalorion_just_align.py [Options...]

        Options must include

        -f, --fastq	path/to/fastq_file (make sure it is not called 2D.fastq, or it will be overwritten)
        -g, --gmap_genome	gmap genome name (e.g. mm10 or hg38, depending on what you named it when you built it) 
        -q, --quality_cutoff	sets the minimum average q-score for a read to be kept (We use '9')
        -t, --gmap_threads	sets the number of threads to run gmap with, default is 1 

    Both wrapper scripts generate a 'content_file' in the directory the fastq file was located in. 
    You will need this file to run the next part of the pipeline. 
    If you want to combine several multiplexed runs for isoform analysis, simply combine their content_files using 'cat' for example

After you aligned your reads either way you can go ahead and 'Define and Quantify Isoforms':

    This part of the pipeline assumes that the reads were processed with the previous part of the pipeline. 
    Read names have to be free of spaces and contain information on whether a read contained ISPCR adapters (e.g. end with '_p_p' if they did on both ends) 
    The content_file links this part of the pipeline to the previous part. 
    It contains in tab delimited format one line for each sample you want to analyze.
    The scripts will use all lines in your content_file to create isoform features, then create and quantify isoforms for each line separately.

    The output of this part of the pipeline should be the following files:

        In the path specified under 2 below:

            TESS.bed    -> Contains all putative Transcription Start and End Sites (TSS/TES) in BED format. 
                       Column 4 contains Site information required for downstream analysis.
            SS.bed    -> Contains all putative Splice Sites (SS)  in BED format. 
                     Column 4 contains Site information required for downstream analysis.
            Combined_TESS_SS.txt    -> Contains groups of TSS/TES and SS that appear in the same reads. 
                                   One line per group.    
            Matched_Combined_TESS_SS.txt    -> Contains TSS/TES/SS matched to genes. 
                                           Same format as Combined_TESS_SS, just Group name in first columnn is modified if group matches gene.
            Gene_Expression.txt    -> Gene level expression values (RP10K) for all samples in 1 below. 
                                  Format(tab separated): 
                                  Gene_name+'_'+chromosome,chromosome,start,end,combined_exon_length,expression_sample_1,expression_sample_2,...
            Isoform_expression.txt    -> Contains expression values (RP10K) for all isoforms of a gene in a single line for each sample in 1 below. 
                                         Format(tab separated):
                                         1.column    gene_names separated by '-' followed by group_numbers+'_'+chromosome
                                         2.column    total expression in all isoform of gene in all samples separated by commas
                                         3+.columns  Isoform number+'_'+expression of isoform in all samples separated by commas

        In each sample_specific_output_folder in the content:

            Isoform_Consensi_aligned.fasta    -> Contains consensus sequences for all genes specified in 5 below
            Isoform_Consensi_gmapoutput_filtered.psl    -> Contains alignment of consensus sequences to the genome specified in 4 below
 
    $python3 Mandalorion_define_and_quantify_isoforms.py [Options...]

    Options must include

    -c, --content_file	path/to/content_file (generated by first part of pipeline) 
    -p, --path	path/to/where_you_want_your_outputfiles 
    -a, --genome_annotation	path/to/genome_annotation_file (Has to be a gtf file and contain "exon" features and "gene_id" fields) 
    -g, --gmap_genome	gmap genome name (e.g. mm10 or hg38, depending on what you named it when you built it) 
    -l, --gene_list	path/to/list_of_genes_for_consensi (Consensus sequences are generated for genes in this file. Example provided)
    -i, --illumina_content_file	file containing paths to .psl files containing illumina read alignments. One path per line
    -r, --refine {g,gi,i}	if 'g' is specified, the genome annotation file is used to populate TESS and SS bins. 
		If 'i' is specified, illumina reads are used to refine splice junctions. 

