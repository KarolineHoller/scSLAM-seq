import sys
import pdb


## this calculates the gene count per cell and the UMI count per cell based on the UMI-stats file


ifn = sys.argv[1]
ofn = sys.argv[1][:-4] + "UMIs+genes_per_cell.csv"


#these count UMIs and genes per cell
UMI_dict = {}
gene_dict = {}



lines_processed = 0
with open(ifn, "r") as umi_stats:
    next(umi_stats)
    for line in umi_stats:
        lines_processed+=1
        if lines_processed % 10000 == 0:
            sys.stdout.write('\r' + str(lines_processed) + ' UMIs processed.')
        lin = line.rstrip().split("\t")
        CB = lin[1]
        #this is counting UMIs per cell. Since the UMI_stats has one line per UMI we can sinply count the lines corresponding to each CB
        if not CB in UMI_dict:
            UMI_dict[CB] = 1
        else:
            UMI_dict[CB]+=1
        #here we are counting genes, which is a bit more elaborate
        gene = lin[3]
        if not CB in gene_dict:
            gene_dict[CB] = []
            gene_dict[CB].append(gene)
        else:
            if not gene in gene_dict[CB]:
                gene_dict[CB].append(gene)


print("writing...")

with open(ofn, "w") as out:
    header = "CB\tUMIs\tgenes\n"
    out.write(header)
    for CB in UMI_dict:
        l = CB + "\t" + str(UMI_dict[CB]) + "\t" + str(len(gene_dict[CB])) + "\n"
        out.write(l)
