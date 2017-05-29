#!/usr/bin/env python

from sys import argv
import subprocess
from os import path
import logging
#import datetime

"""
This script will take bam or sam files from CLC and try to convert them to a format that is more compatible with
third part software that wants more bwa-like bam/sam files (such as manta for example). Work In Progress.
"""


# Function splits input into perfectly mapped reads with NH:1 and QMAP > 3 and secondary mapped
# with NH:>1 and QMAP < 3


def sam_split(samfile_in, out_perfect, out_secondary):
    with open(samfile_in, "r") as sam:
        with open(out_perfect, "w+") as perfect:
            with open(out_secondary, "w+") as secondary:
                it1 = 0
                it2 = 0
                perfs = 0
                secs = 0
                for line in sam:
                    old_line = line.split("\t")
                    # new_line = old_line
                    if not line.startswith("@"):
                        NH_field = old_line[-2].split(":")[-1]
                        try:
                            if int(old_line[4]) <= 3 and int(NH_field) > 1 or bin(int(old_line[1]))[-2] == '0':
                                secondary.write(line)
                                secs += 1
                                # new_line[1] = str(int(old_line[4]) + 256)

                            else:
                                perfect.write(line)
                                perfs += 1
                        except TypeError as te:
                            it1 += 1
                            logging.warning("type error number {}".format(str(it1)))
                            print("column 4: ", old_line[4])
                            print("NH_field: ", NH_field)
                            print("Column 1: ", old_line[1])
                            print(str(te))
                            continue
                        except ValueError as ve:
                            it2 += 1
                            # logging.warning("value error number {}".format(str(it2)))
                            print(str(ve))
                            if bin(int(old_line[1]))[-2] == "b":
                                # logging.info("This read has flag 0 (mapped, unpaired). Sending it to perfect.")
                                perfect.write(line)
                                perfs += 1
                            continue
                logging.info("{} perfect. {} secondary. {} TypeErrors. {} ValueErrors".format(perfs, secs, it1, it2))
                # print("Column 4 (MAPQ):" + old_line[4])
                # print("NH" + NH_field)
                # print(bin())
                # print("This line will not be used:")
                # print(line)
                # new_line = "\t".join(new_line)
                # print new_line
                # yield str(new_line)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, filename="bam_fixer_2.log", filemode="w",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    infile = argv[1]
    outfile_perfect = argv[2]
    outfile_secondary = argv[3]
    fasta_index = "/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta.fai"
    bwa_index = fasta_index.rsplit('.', 1)[0]
    directory = path.dirname(outfile_secondary)

    samtools_module = "module load samtools/1.3.1"
    samtools_source = "source /apps/bio/apps/samtools/1.3.1/samtools"
    bwa_module = "module load bwa/0.7.5a"
    logging.info("initial variables set")
    # subprocess.call("module load samtools/1.3.1", shell=True)
    # subprocess.call("module load bwa/0.7.5a", shell=True)

    # Determine if the inut is bam or sam
    if ".bam" in infile:
        logging.info("Input file was: BAM")
        intermediary = infile.replace(".bam", ".sam")
        print(str(infile))
        print(str(intermediary))
        try:
            subprocess.check_call(
                "{} && samtools sort {} -n -@ 112 -m 2G | samtools view - -o {} -@ 112 -h".format(samtools_source, infile,
                                                                                               intermediary),
                shell=True)
        except subprocess.CalledProcessError:
            logging.warning("CALLEDPROCESSERROR in initial sort")
        except OSError:
            print("OSERROR")
            logging.warning("OSERROR in initial sort")
        logging.info("Initial conversion of BAM to SAM (name sorted) completed. Now starting SAM-file splitting")
        sam_split(intermediary, outfile_perfect, outfile_secondary)
        # subprocess.call(["rm", intermediary])
        logging.info("Splitting of the SAM-file completed")
    elif ".sam" in infile:
        logging.info("input file was: SAM")
        temp = infile.replace(".sam", "_tmp.sam")
        try:
            #subprocess.call("", shell=True)
            subprocess.check_call(
                "{} && samtools view {} -@ 112 -ht {} | samtools sort - -@ 40 -m 2G -n | samtools view - -o {} -@ 112 -h".format(
                    samtools_module, infile, fasta_index, temp), shell=True)
        except subprocess.CalledProcessError:
            logging.warning("CALLEDPROCESERROR in initial sort")
        except OSError:
            logging.warning("OSERROR in initial sort")
        logging.info("Initial name sorting completed, starting SAM-file splitting")
        sam_split(temp, outfile_perfect, outfile_secondary)
        logging.info("Splitting of the SAM-file completed")
        # subprocess.call(["rm", temp])

    perfsize = path.getsize(outfile_perfect)
    secsize = path.getsize(outfile_secondary)
    logging.info("After splitting: ByteSize of Perfect: {}. Bytesize of secondary: {}".format(perfsize, secsize))

    # old and unused
    # headers = str(subprocess.call("tr '\t' '\n' < %s | grep RG: | sort | uniq" % outfile_secondary, shell=True, stdout=subprocess.PIPE))
    # with open(outfile_secondary + "tmp", "w+") as sec_tmp:
    #    sec_tmp.write(headers)
    #    sec_tmp.write(outfile_secondary)
    #    subprocess.call("mv %s %s" % (sec_tmp, outfile_secondary))

    # Get original headers
    original_headers = directory + "/original_headers.sam"
    subprocess.call("{} && samtools view -H {} > {}".format(samtools_module, infile, original_headers), shell=True)
    logging.info("Headers extracted from original input file: {} to : {}".format(infile, original_headers))

    secondary_tmp = outfile_secondary + "_temp"

    # Reheader the secondary mapped sam file
    # Original headers used to be the .fai, change back if crash
    subprocess.call(
        "%s && samtools view -ht %s %s > %s" % (samtools_module, original_headers, outfile_secondary, secondary_tmp),
        shell=True)
    logging.info("Temporary secondary SAM-file: {} re-headered into : {}".format(outfile_secondary, secondary_tmp))

    # Rename tempfile to original name
    subprocess.call(["mv", secondary_tmp, outfile_secondary])
    logging.info("SAM-file : {} overwritten by temp file {}, which no longer exists".format(outfile_secondary,
                                                                                            secondary_tmp))

    # Convert to fastq
    subprocess.call(
        "%s && samtools fastq %s > %s/reads_interleaved.fastq" % (samtools_module, outfile_secondary, directory),
        shell=True)
    logging.info("Converted secondary to fastq")

    # Run bwa
    subprocess.call("%s && bwa mem %s -p %s/reads_interleaved.fastq -t 112 > %s/bwa_out.sam" % (bwa_module, bwa_index,
                                                                                                directory, directory),
                    shell=True)
    logging.info("Secondary reads re-mapped")

    # Variable assignment
    bwa_out = directory + "/bwa_out.sam"
    merged_bam = directory + "/merged.bam"
    sorted_bam = directory + "/merged_sorted.bam"
    new_headers = directory + "/new_headers.sam"
    logging.info("Merger-related variables set: {}, {}, {}, {} (not used)".format(bwa_out, merged_bam, sorted_bam,
                                                                                  new_headers))

    # Merge perfect mapped into the bwa output
    subprocess.call("cat %s >> %s" % (outfile_perfect, bwa_out), shell=True)
    logging.info("Perfect SAM merged into remapped secondary SAM")

    # Convert the merged sam file into bam
    subprocess.call("%s && samtools view -b -@ 112 %s -o %s" % (samtools_module, bwa_out, merged_bam), shell=True)
    logging.info("Merged SAM converted to BAM")

    # Extract the header from the new bam file. Probably unnecessary
    #subprocess.call("%s && samtools view -H %s > %s" % (samtools_module, merged_bam, original_headers), shell=True)
    #logging.info("Headers extracted")

    # Reheader the merged bam file
    rehead_tmp = directory + "/temporary_reheader.bam"
    subprocess.call("%s && samtools reheader %s %s > %s" % (samtools_module, original_headers, merged_bam, rehead_tmp),
                    shell=True)
    logging.info("BAM file: {} reheadered using: {} into tempfile: {}".format(merged_bam, original_headers, rehead_tmp))

    subprocess.call(["mv", rehead_tmp, merged_bam])
    logging.info("Non-headered bamfile: {} overwritten by headered tempfile: {}".format(merged_bam, rehead_tmp))

    # Sort reheadered bam file
    # with open(sorted_bam, "wb") as sorted:
    subprocess.call("/apps/bio/apps/samtools/1.3.1/samtools sort -@ 112 -m 2G %s > %s" % (merged_bam, sorted_bam),
                    shell=True)
    logging.info("Final BAM: {} sorted to: {}".format(merged_bam, sorted_bam))
    subprocess.call(["/apps/bio/apps/samtools/1.3.1/samtools", "index", sorted_bam])
    logging.info("Final BAM: {} indexed into: {}".format(sorted_bam, sorted_bam + ".bai"))

    # Remove all intermediary files. Fill in later when testing is done.
    # subprocess.call(["rm", ])

    # Give path to result
    print("Location of output file: \n" + str(path.abspath(sorted_bam)))
    logging.info("Everything is Completed.")
