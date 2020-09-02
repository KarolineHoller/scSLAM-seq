import pysam
import re
import pdb
import sys
import os
#import pandas as pd

## input has to be sorted and indexed
# Author: Anika Neuschulz


if len(sys.argv) == 2:
    ifn = sys.argv[1]

else:
    print("Please make my life easier by providing an input file (indexed .bam).")
    sys.exit()

ofn = ifn[:-4] + "_MTabs.sam"


scriptpath = os.path.abspath(sys.argv[0])
scriptdir = "/".join(scriptpath.split("/")[:-1])





samfile = pysam.AlignmentFile(ifn, "rb")

#this creates a .sam file, that then has to be converted to .bam, sorted and indexed


mappedreads = pysam.AlignmentFile(ofn, "w", template=samfile)




for read in samfile.fetch():

    total_content = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}
    mutations = ["0-XX"]
    sequence = read.get_reference_sequence()
    for base in total_content.keys():
        total_content[base] += sequence.count(base)
        TC_tag = ''.join([''.join(key) + str(total_content[key]) + ';' for key in total_content.keys()])[:-1]
    if int(read.get_tag("NM")) > 0 and not re.search("[ID]", read.cigarstring):
        md = re.findall('(\d+)(\^[A-Za-z]+|[A-Za-z])', read.get_tag("MD"))
                # adjust md for softclipped reads
        softclip_test = re.match("\d+S", read.cigarstring)
        if softclip_test:
            try:
                md_adjustment = softclip_test[0][0:len(softclip_test[0])-1]
                md[0] = tuple((str(int(md_adjustment) + int(md[0][0])), md[0][1]))
            except IndexError:
                pdb.set_trace()
                #lets go from relative to absolute positions for the MD
                #first MD number does not need to be altered, all others need to have 1 added for position
                #of mutated base, therefore this cheap -1 trick
                # to go for globally absolute values, let's account for the absolute mapping position on the chromosome

        contig = read.reference_name
        mapping_start = read.reference_start
        cumsum_start = -1 + mapping_start
        mutations = []
        for element in md:
            try:
                base = element[1]
                position = int(element[0]) + cumsum_start + 1
                mutations.append(str(contig) + ":" + str(position) + "-" + base + read.seq[position-mapping_start])
                cumsum_start = position
            except IndexError:
                pass
    if re.search("[ID]", read.cigarstring):
        mutations = ["0-ID"]
    read.tags = read.tags + [('MT','_'.join(mutations))]

    try:

        read.set_tag('TC',TC_tag,'Z')

        mappedreads.write(read)
    except KeyError:
        pass

mappedreads.close()
