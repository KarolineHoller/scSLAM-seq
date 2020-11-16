# This repository contains scripts used to analyze scSLAM-seq data of zebrafish embryos.

## How to run the pipeline:  

Instructions for all scripts can be obtained by calling the script without parameters.
Pysam (https://github.com/pysam-developers/pysam) is required to run the pipeline.
Samtools (http://www.htslib.org/) is also required.


1. map demultiplexed data with cellranger using defauls parameters (tested with release 3.0.2)

2. unzip ```outs/filtered_feature_bc_matrix/barcodes.tsv.gz```

3. run ___script 1___

4. zip ```outs/filtered_feature_bc_matrix/barcodes.tsv``` again for later use with Seurat if desired

5. compress, sort and index the .sam file produced by script 1 using samtools

6. run ___script 2___ on the file obtained in step 5

7. compress, sort and index the .sam file produced by script 2 using samtools

8. if desired, bulk labeling rates can be analyzed at this point using ___script 4___ on the result of step 7. A sequencing quality cutoff for a mutation to be considered valid can be set. For analysis purposes nucleotides at the beginning or the end of reads can also be disregarded in counting ('clipped'). This has no effect on the source file or any step downstream. If any clipping of the reads turns out to be desired, please do so using the software of your choice and start the pipeline again.  
This step produces 2 output files that can be utilized to plot mutation efficiencies. **"mutation_occurences"** has the number of occurences for each possible mutation and **"nucleotide_counts"** has the total counts for each nucleotide.  
Additionally a **"mutation_locations"** file is written that has information in which position on the read a mutation occured and if the read was forward (F) or reverse (R). This is a large file, but can be utilized to visualize distribution of mutations along the reads.

9. run ___script 3___ on the file obtained in step 7. Minimum sequencing quality of T to C mutations as well as minimum mutations per UMI can be set. Reads will be separated on a UMI basis.  
Three output files will be created:  
**"labeled"** has the reads from all UMIs that passed the minimum number of qualifying T to C mutations,  
**"maybe-labeled"** has the reads from all UMIs that did not pass the threshold, but still had more that 0 qualifying mutations. This file is empty when separationg with 1 qualifying mutation.  
Finally **"unlabeled"** contains the reads from all UMIs that have no qualitying T to C mutations.  
This step also creates two metric files:  
**"labeled_unlabeled_stats"** contains information on the number of labeled and unlabeled reads for each gene  
**"umi_stats"** contains extensive statistics on gene, cell barcode, number of reads, number of T to C mutations and number of unlabeled reads for each UMI. This is a large file, but useful when analyzing labeling information on a per-cell basis.

10. compress, sort and index the .sam files produced by script 3 using samtools

11. prepare fastq files for separating reads into labeled and unlabeled:
```
cat {read1}.fastq | paste -d ' ' - - - - | sort -k1,1 > R1_cat.fastq
cat {read2}.fastq | paste -d ' ' - - - - | sort -k1,1 > R2_cat.fastq
```
This needs to be done as the fastq format uses 4 lines for each read. To work on a line-by-line basis those 4 lines need to be pasted onto one.

12. obtain readnames of labeled and unlabeled reads
```
cut -f1 {_labeled/_unlabeled}.sam | sort -k1,1 | awk '$0="@"$0' > read_names{_labeled/_unlabeled}.txt
```

13. separate reads  
```
join read_names{_labeled/_unlabeled}.txt {R1/R2}_cat.fastq > select_{_labeled/_unlabeled}_cat_{R1/R2}.fastq
awk '{printf("%s %s\n%s\n%s\n%s\n", $1, $2, $3, $4, $5)}' select_{_labeled/_unlabeled}_cat_{R1/R2}.fastq > select_{_labeled/_unlabeled}_{R1/R2}.fastq
```

14. remap the resulting fastq files with STARsolo to obtain count matrices for use in Seurat. Keep in mind that the input files were filtered for real cells already in step 3, so the raw count matrices contain only valid cells.


-----

## to obtain data on labeling rates per cell:

run the script ___labeling_eff_per_cell_MTabs_final.py___ on a sorted and indexed .bam file that has the MT tag added (output of Step 7 above). A sequencing quality threshold for T to C mutations needs to be provided as a second argument.

## to obtain data in labeling rates per gene:

run the script ___mutation_eff_per_gene_v4_MTabs_final.py___ on a sorted and indexed .bam file that has the MT tag added (output of Step 7 above). A sequencing quality threshold for T to C mutations needs to be provided as a second argument.


## Single-cell analysis of labeled und labeled RNA:

For this analysis, use the count matrices from step 14. Alternatively, download the count matrices from GEO (accession number GSE158849).
The data is clustered separately with the information of old and new RNA, respectively. The code produces the figures included in the manuscript.
