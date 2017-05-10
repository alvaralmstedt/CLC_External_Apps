#!/usr/bin/env python

from sys import argv
import subprocess
from os import path

"""
This script will take bam or sam files from CLC and try to convert them to a format that is more compatible with
third part software that wants more bwa-like bam/sam files (such as manta for example). Work In Progress.
"""

def sam_split(samfile_in, out_perfect, out_secondary):
    with open(samfile_in, "r") as sam:
        with open(out_perfect, "w+") as perfect:
            with open(out_secondary, "w+") as secondary:
                for line in sam:
                    old_line = line.split("\t")
                    #new_line = old_line
                    NH_field = old_line[12].split(":")[-1].rstrip()
                    if int(old_line[4]) <= 3 and int(NH_field) > 1:
                        secondary.write(line)
                        #new_line[1] = str(int(old_line[4]) + 256)
                    elif not line.startswith("@"):
                        perfect.write(line)
                    else:
                        print("This line will not be used:")
                        print(line)
                    #new_line = "\t".join(new_line)
                    #print new_line
                    #yield str(new_line)

if __name__ == "__main__":
    infile = argv[1]
    outfile_perfect = argv[2]
    outfile_secondary = argv[3]
    fasta_index = "/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta.fai"
    bwa_index = fasta_index.rsplit('.', 1)[0]
    directory = path.dirname(outfile_secondary)

    subprocess.call("module load samtools/1.3.1", shell=True)
    subprocess.call("module load bwa/0.7.5a", shell=True)
    if ".bam" in infile:
        intermediary = infile.replace(".bam", ".sam")
        print(str(infile))
        print(str(intermediary))
        subprocess.call("samtools view %s -o %s -@ 40" % (infile, intermediary), shell=True)
#        subprocess.call(["wait"])
        sam_split(intermediary, outfile_perfect, outfile_secondary)
    elif ".sam" in infile:
        sam_split(infile, outfile_perfect, outfile_secondary)

    #headers = str(subprocess.call("tr '\t' '\n' < %s | grep RG: | sort | uniq" % outfile_secondary, shell=True, stdout=subprocess.PIPE))
    #with open(outfile_secondary + "tmp", "w+") as sec_tmp:
    #    sec_tmp.write(headers)
    #    sec_tmp.write(outfile_secondary)
    #    subprocess.call("mv %s %s" % (sec_tmp, outfile_secondary))
    secondary_tmp = outfile_secondary + "_temp"
    subprocess.call("samtools view -ht %s %s > %s" % (fasta_index, outfile_secondary, secondary_tmp), shell=True)
    subprocess.call(["mv", secondary_tmp, outfile_secondary])
    subprocess.call("samtools fastq %s > %s/reads_interleaved.fastq" % (outfile_secondary, directory), shell=True)
    subprocess.call(["bwa", "mem", bwa_index, "-p", directory + "/reads_interleaved.fastq", "-t", "40", ">", directory + "/bwa_out.sam"])
    bwa_out = directory + "/bwa_out.sam"
    merged_bam = directory + "/merged.bam"
    sorted_bam = directory + "/merged_sorted.bam"
    subprocess.call("cat %s >> %s" % (outfile_perfect, bwa_out), shell=True)
    subprocess.call("samtools view -b -@ 40 %s -o %s" % (bwa_out, merged_bam), shell=True)
    
    with open(sorted_bam, "w+") as sorted:
        sorted.write(subprocess.check_output(["samtools", "sort", "-@," "30", "-m", "2G,", str(merged_bam)]))
    subprocess.call(["samtools", "index", sorted_bam])
    print("Location of output file: \n" + str(path.abspath(sorted_bam)))
