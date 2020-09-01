import pysam
import re
import pdb
import sys
import os

#this separates labeled and unlabeled reads based on a T->C event over a certain quality threshold
# this is UMI aware - all reads from one BC-UMI-gene combination go into labeled, if one is considered labeled
# this also creates stats on how many umis per gene have how many labeling events
# this only considers reads that map to a gene

# this works using multiple dictionaries:
## gene_dict: numbers of labeled/unlabeled reads for each gene
## umi_stats: number of reads, labeling events, unlabeled reads for each barcode_gene_umi combination
## labeled_dict: number of labeling events for each gene_barcode_UMI combination for read separation purposes
## umi_dict: has all reads for each gene_barcode_UMI combination to be distributed to labeled/unlabeled files
# Author: Anika Neuschulz

if len(sys.argv) == 4:
    ifn = sys.argv[1]
    qual_thr = int(sys.argv[2])
    nmutations = int(sys.argv[3])
else:
    print("Please provide me with an input .bam file to extract labeled reads from, a minimum quality score for mutations and the number of required mutations for a UMI.")
    sys.exit()


scriptpath = os.path.abspath(sys.argv[0])
scriptdir = "/".join(scriptpath.split("/")[:-1])



ofn_l = ifn[:-4] + "_labeled_Q" + str(qual_thr) + "_" + str(nmutations) + "mutations.sam"
ofn_ml = ifn[:-4] + "_maybe-labeled_Q" + str(qual_thr) + "_" + str(nmutations) + "mutations.sam"
ofn_ul = ifn[:-4] + "_unlabeled_Q" + str(qual_thr) + "_" + str(nmutations) + "mutations.sam"

samfile = pysam.AlignmentFile(ifn, "rb")

gene_dict = {}
#structure: [#unlabeled, #labeled]

umi_stats = {}
#structure: [barcode_umi_gene][barcode][umi][gene][n reads][n T->C][n unlabeled]

labeledreads = pysam.AlignmentFile(ofn_l, "w", template=samfile)
unlabeledreads = pysam.AlignmentFile(ofn_ul, "w", template=samfile)
maybelabeledreads = pysam.AlignmentFile(ofn_ml, "w", template=samfile)

reads_processed = 0
sys.stdout.write(str(reads_processed) + " reads processed")

labeled_dict = {}
umi_dict = {}

for read in samfile.fetch():
    reads_processed+=1
    sys.stdout.write("\r" + str(reads_processed) + " reads processed.")
    try:
        CB = read.get_tag("CB")
        gene = read.get_tag("GN")
        UB = read.get_tag("UB")
        identifier = CB + "_" + UB + "_" + gene
        MT = read.get_tag("MT")

        #remove ID reads here
        if MT == "0-ID":
            continue #this skips the loop for all reads with insertions/deletions
            #those reads are therefore not included in the output

        #a umi is considered unlabeled until a labeled read is found
        if not identifier in labeled_dict:
            labeled_dict[identifier] = 0

        #this is where we collect all reads from a umi
        if not identifier in umi_dict:
            umi_dict[identifier] = [read]
        else:
            umi_dict[identifier].append(read)

        #this is where we count reads for each umi
        if not identifier in umi_stats:
            umi_stats[identifier] = [str(CB),str(UB),str(gene),1,0,0]
        else:
            umi_stats[identifier][3]+=1

        #here we check if a read is labeled
        if read.get_tag("NM") != 0:
            mutations = tuple(MT.split("_"))
            if read.is_reverse:
                label = "AG"
            else:
                label = "TC"
            passed = 0
            for mutation in mutations:
                split_mutation = mutation.split("-")
                #pdb.set_trace()
                position = int(int(split_mutation[0].split(':')[1]) - int(read.reference_start))
                if split_mutation[1] == label:
                    if ord(read.qual[position]) > qual_thr + 33:
                        passed+=1
                        #count labeling events
                        umi_stats[identifier][4]+=1
                        #no "else" that sets passed to false, as I assume one mutation over the threshold to be sufficient to assume the read as labeled
            if passed >= 1:
                #set umi marker to labeled = True
                labeled_dict[identifier]+=passed
                try:
                    gene = read.get_tag("GX")
                    if gene in gene_dict:
                        gene_dict[gene][1]+=1
                    else:
                        gene_dict[gene] = [0,1]
                except KeyError:
                    pass
                    #catch stuff that does not map to a gene, we don't use it in the count matrix anyway

            else:
                #count unlabeled reads
                umi_stats[identifier][5]+=1
        else:
            try:
                #count unlabeled reads
                umi_stats[identifier][5]+=1
                gene = read.get_tag("GX")
                if gene in gene_dict:
                    gene_dict[gene][0]+=1
                else:
                    gene_dict[gene] = [1,0]
            except KeyError:
                pass
            #catch stuff that does not map to a gene, we don't use it in the count matrix anyway
    except KeyError:
        pass
        #some reads have no corrected UMI, CB or gene (GX)

sys.stdout.write("\nwriting reads.")
for identifier in labeled_dict:
    if labeled_dict[identifier] >= nmutations:
        for read in umi_dict[identifier]:
            labeledreads.write(read)
    elif labeled_dict[identifier] == 0:
        for read in umi_dict[identifier]:
            unlabeledreads.write(read)
    else:
        for read in umi_dict[identifier]:
            maybelabeledreads.write(read)

labeledreads.close()
unlabeledreads.close()
maybelabeledreads.close()

statfile = ifn[:-4] + "_labeled_unlabeled_stats" + "_Q" + str(qual_thr) + "_" + str(nmutations) + "mutations.csv"
sys.stdout.write("\nwriting gene stats.")
with open(statfile, "w") as f:
    f.write("gene\tunlabeled\tlabeled\n")
    for gene in gene_dict:
        f.write(gene + "\t" + str(gene_dict[gene][0]) + "\t" + str(gene_dict[gene][1]) + "\n")

sys.stdout.write("\nwriting UMI stats.")
umi_statfile = ifn[:-4] + "_umi_stats" + "_Q" + str(qual_thr) + "_" + str(nmutations) + "mutations.csv"
with open(umi_statfile, "w") as g:
    #[barcode][umi][gene][n reads][n T->C][n unlabeled]
    g.write("Identifier\tCB\tUB\tgene\tnReads\tnT>C\tnUnlabeled\n")
    for identifier in umi_stats:
        g.write(identifier + "\t" + str(umi_stats[identifier][0]) + "\t" + str(umi_stats[identifier][1]) + "\t" + str(umi_stats[identifier][2]) + "\t" + str(umi_stats[identifier][3]) + "\t" + str(umi_stats[identifier][4]) + "\t" + str(umi_stats[identifier][5]) + "\n")
