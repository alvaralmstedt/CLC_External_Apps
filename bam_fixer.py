#!/usr/bin/env python

from subprocess import call
from sys import argv
import pysam
import os

#def read_input(sam_in):
#    sam = open(sam_in, "r")
#    return sam

def sam_logic(samfile):
    with open(samfile) as sam:
        for line in sam:
            if not line.startswith("@")
                old_line = line.split("\t")
                new_line = old_line
                NH_field = old_line[12].split(":")[-1].rstrip()
                if int(old_line[4]) <= 3 and int(NH_field) > 1 :
                    new_line[1] = str(int(old_line[4]) + 256)
                new_line = "\t".join(new_line)
                print new_line
                yield str(new_line)
            else:
                continue

def write_sam(input_line, output_sam):
    with open(output_sam, "a+") as outfile:
        for element in input_line:
            outfile.write(element)

if __name__ == "__main__":

    in_bam = argv[1]
    out_sam = argv[2]
#    in_sam = "./{}".format(os.path.basename(in_bam).strip(".bam") + ".sam")
#    pysam.view("-ho {} {}".format(in_bam, in_sam))
#    test_1 = sam_logic(in_bam)
    write_sam(sam_logic(in_bam), out_sam)
