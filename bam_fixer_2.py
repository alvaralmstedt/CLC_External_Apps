#!/usr/bin/env python

from sys import argv
import subprocess
from os import path

"""
This script will take bam or sam files from CLC and try to convert them to a format that is more compatible with
third part software that wants more bwa-like bam/sam files (such as manta for example). Work In Progress.
"""

#Function splits input into perfectly mapped reads with NH:1 and QMAP > 3 and secondary mapped
#with NH:>1 and QMAP < 3
def sam_split(samfile_in, out_perfect, out_secondary):
    with open(samfile_in, "r") as sam:
        with open(out_perfect, "w+") as perfect:
            with open(out_secondary, "w+") as secondary:
                it1 = 0
                it2 = 0
                for line in sam:
                    old_line = line.split("\t")
                    #new_line = old_line
                    if not line.startswith("@"):    
                        try:
                            NH_field = old_line[-1].split(":")[-1].rstrip()
                            if int(old_line[4]) <= 3 and int(NH_field) > 1 or int(bin(old_line[1])[-2]) == 1:
                                secondary.write(line)
                                #new_line[1] = str(int(old_line[4]) + 256)
                            else:    
                                perfect.write(line)
                        except TypeError as te:
                            it1 += 1
                            print("type error number {}".format(str(it1)))
                            print(str(te))
                            continue
                        except ValueError as ve:
                            it2 += 1
                            print("value error number {}".format(str(it2)))
                            print(str(ve))
                            continue
                        
                            #print("Column 4 (MAPQ):" + old_line[4])
                            #print("NH" + NH_field)
                            #print(bin())
                        #print("This line will not be used:")
                        #print(line)
                    #new_line = "\t".join(new_line)
                    #print new_line
                    #yield str(new_line)

if __name__ == "__main__":
    #
    infile = argv[1]
    outfile_perfect = argv[2]
    outfile_secondary = argv[3]
    fasta_index = "/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta.fai"
    bwa_index = fasta_index.rsplit('.', 1)[0]
    directory = path.dirname(outfile_secondary)
    
    samtools_module = "module load samtools/1.3.1"
    bwa_module = "module load bwa/0.7.5a"
    
    #subprocess.call("module load samtools/1.3.1", shell=True)
    #subprocess.call("module load bwa/0.7.5a", shell=True)
    
    #Determine if the inut is bam or sam
    if ".bam" in infile:
        intermediary = infile.replace(".bam", ".sam")
        print(str(infile))
        print(str(intermediary))
        subprocess.call("{} && samtools sort {} -n -@ 112 -m 2G | samtools view - -o {} -@ 112".format(samtools_module, infile, intermediary), shell=True)
        sam_split(intermediary, outfile_perfect, outfile_secondary)
        subprocess.call(["rm", intermediary])
    elif ".sam" in infile:
        temp = infile.replace(".sam", "_tmp.sam")
        subprocess.call("{} && samtools view {} -@ 112 | samtools sort - -@ 40 -m 2G -n | samtools view - -o {} -@ 112".format(samtools_module, infile, temp), shell=True)
        sam_split(temp, outfile_perfect, outfile_secondary)
        subprocess.call(["rm", temp])
    #headers = str(subprocess.call("tr '\t' '\n' < %s | grep RG: | sort | uniq" % outfile_secondary, shell=True, stdout=subprocess.PIPE))
    #with open(outfile_secondary + "tmp", "w+") as sec_tmp:
    #    sec_tmp.write(headers)
    #    sec_tmp.write(outfile_secondary)
    #    subprocess.call("mv %s %s" % (sec_tmp, outfile_secondary))
    secondary_tmp = outfile_secondary + "_temp"
    
    #Reheader the secondary mapped sam file
    subprocess.call("%s && samtools view -ht %s %s > %s" % (samtools_module, fasta_index, outfile_secondary, secondary_tmp), shell=True)
    
    #Rename tempfile to original name
    subprocess.call(["mv", secondary_tmp, outfile_secondary])
    
    #Convert to fastq
    subprocess.call("%s && samtools fastq %s > %s/reads_interleaved.fastq" % (samtools_module, outfile_secondary, directory), shell=True)
    
    #Run bwa
    subprocess.call("%s && bwa mem %s -p %s/reads_interleaved.fastq -t 112 > %s/bwa_out.sam" % (bwa_module ,bwa_index, directory, directory), shell=True)
    
    #Variable assignment
    bwa_out = directory + "/bwa_out.sam"
    merged_bam = directory + "/merged.bam"
    sorted_bam = directory + "/merged_sorted.bam"
    original_headers = directory + "/original_headers.sam"

    #Merge perfect mapped into the bwa output
    subprocess.call("cat %s >> %s" % (outfile_perfect, bwa_out), shell=True)
    
    #Convert the merged sam file into bam
    subprocess.call("%s && samtools view -b -@ 112 %s -o %s" % (samtools_module, bwa_out, merged_bam), shell=True)
    # ---We are here---

    #Extract the header from the original bam file
    subprocess.call("%s && samtools view -H %s > %s" % (samtools_module, merged_bam, original_headers), shell=True)

    #Reheader the merged bam file
    subprocess.call("%s && samtools reheader -i %s %s" % (samtools_module, original_headers, merged_bam), shell=True)
    
    #Sort reheadered bam file
    with open(sorted_bam, "w+") as sorted:
        subprocess.check_output(["/apps/bio/apps/samtools/1.3.1/samtools", "sort", "-@", "112", "-m", "2G", merged_bam], stdout=sorted)
    subprocess.call(["/apps/bio/apps/samtools/1.3.1/samtools", "index", sorted_bam])
    
    #Give path to result
    print("Location of output file: \n" + str(path.abspath(sorted_bam)))
