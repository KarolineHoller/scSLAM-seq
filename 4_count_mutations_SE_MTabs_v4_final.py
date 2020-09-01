import pysam
import re
import pdb
import string
import sys
import os

# this now works with absolute MT tags, on the other hand unmutated reads and indels are
# not counted anymore
# Author: Anika Neuschulz



if len(sys.argv) == 5:
    ifn = sys.argv[1]
    qcut = int(sys.argv[2])
    crop = int(sys.argv[3])
    cropside = sys.argv[4]
else:
    print("Please give me an input file to excercise my counting abilities, a desired quality cutoff, how many nucleotides to crop and if cropping should happen at the beginning (B) and/or the end (E) of the read or not at all (N).")
    sys.exit()


scriptpath = os.path.abspath(sys.argv[0])
scriptdir = "/".join(scriptpath.split("/")[:-1])



bamfile = pysam.AlignmentFile(ifn, "rb")

mutation_location_file = ifn[:-4] + "_mutation_locations_" + str(qcut) + "_" + str(crop) + cropside +".txt"
mlf = open(mutation_location_file, "w")

mutation_dict = {}
mutation_dict["AC"] = 0
mutation_dict["AG"] = 0
mutation_dict["AT"] = 0
mutation_dict["AN"] = 0
mutation_dict["CA"] = 0
mutation_dict["CG"] = 0
mutation_dict["CT"] = 0
mutation_dict["CN"] = 0
mutation_dict["GA"] = 0
mutation_dict["GC"] = 0
mutation_dict["GT"] = 0
mutation_dict["GN"] = 0
mutation_dict["TA"] = 0
mutation_dict["TG"] = 0
mutation_dict["TC"] = 0
mutation_dict["TN"] = 0
#mutation_dict["XX"] = 0
#mutation_dict["ID"] = 0
mutation_dict["NA"] = 0
mutation_dict["NC"] = 0
mutation_dict["NG"] = 0
mutation_dict["NT"] = 0

tab = str.maketrans("ACTG", "TGAC")

countA = 0
countT = 0
countC = 0
countG = 0



for read in bamfile.fetch():

    if not read.is_secondary:

        if not read.is_reverse:

            if cropside == "B":
                sequence = read.seq[crop:]
            elif cropside == "BE":
                sequence = read.seq[crop:-crop]
            elif cropside == "E":
                sequence = read.seq[:-crop]
            elif cropside == "N":
                sequence = read.seq
            else:
                print("Please specify where/if to crop reads for counting. Beginning: B; Beginning+End: BE, End: E, not at all: N.")
                sys.exit()
            countA+= sequence.count("A")
            countT+= sequence.count("T")
            countC+= sequence.count("C")
            countG+= sequence.count("G")

        if read.is_reverse:
            if cropside == "B":
                sequence = read.seq[:-crop]
            elif cropside == "BE":
                sequence = read.seq[crop:-crop]
            elif cropside == "E":
                sequence = read.seq[crop:]
            elif cropside == "N":
                sequence = read.seq
            else:
                print("Please specify where/if to crop reads for counting. Beginning: B; Beginning+End: BE, End: E, not at all: N.")
                sys.exit()
            countT+= sequence.count("A")
            countA+= sequence.count("T")
            countG+= sequence.count("C")
            countC+= sequence.count("C")
        mutations = read.get_tag("MT").split("_")
        ## this is where we catch unmutated reads or insertions/deletions
        for mutation in mutations:
            if mutation == "0-XX" or mutation == "0-ID":
                break
            split_mutation = mutation.split("-")
            # count only mutations above certain quality +33 to include offset from IlluminaQuality to Unicode code points
            #we try to clip by establishing two boundary numbers - mutations have to fall in between in order to be considered
            if read.is_reverse:
                if cropside == "B":
                    cropnumber1 = 0
                    cropnumber2 = len(read.seq)-crop
                elif cropside == "BE":
                    cropnumber1 = crop
                    cropnumber2 = len(read.seq)-crop
                else:
                    cropnumber1 = crop
                    cropnumber2 = len(read.seq)
            else:
                if cropside == "B":
                    cropnumber1 = crop
                    cropnumber2 = len(read.seq)
                elif cropside == "BE":
                    cropnumber1 = crop
                    cropnumber2 = len(read.seq)
                else:
                    cropnumber1 = 0
                    cropnumber2 = len(read.seq)-crop
            # here we adujst for the absolute MT tags
            mapping_start = read.reference_start
            relative_position = int(split_mutation[0].split(':')[1]) - mapping_start

            if ord(read.qual[relative_position]) >= qcut + 33 and relative_position< cropnumber2 and  relative_position > cropnumber1:
                try:
                    if read.is_reverse:

                        split_mutation[1] = split_mutation[1].translate(tab)[::1]
                    mutation_dict[split_mutation[1]]+=1
                    ##write mutation locations
                    if read.is_reverse:
                        direction = "R"
                    else:
                        direction = "F"
                    mlf.write(split_mutation[0]+ "\t" + split_mutation[1] + "\t" + direction + "\n")
                except KeyError:
                    pass


mlf.close()

ofn = ifn[:-4] + "_mutation_occurences_" + str(qcut) + "_" + str(crop) + cropside + ".txt"
outfile = open(ofn, "w")
for k, v in mutation_dict.items():
    outfile.write(str(k) + "\t" + str(v) + "\n\n")
outfile.close()

ofn2 = ifn[:-4] + "_nucleotide_counts_Q" + str(qcut) + "_C" + str(crop) + cropside + ".txt"
outfile = open(ofn2, "w")
outfile.write("A\t" + str(countA) + "\n")
outfile.write("C\t" + str(countC) + "\n")
outfile.write("T\t" + str(countT) + "\n")
outfile.write("G\t" + str(countG) + "\n")
outfile.close()
