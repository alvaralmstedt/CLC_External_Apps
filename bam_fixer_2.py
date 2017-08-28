#!/usr/bin/env python

from sys import argv
import subprocess
from os import path
import logging
import datetime
import sam_parse

"""
This script will take bam or sam files from CLC and try to convert them to a format that is more compatible with
third part software that wants more bwa-like bam/sam files (such as manta for example). Work In Progress.
Aimed to be used for manta
"""


def errchk():
    subprocess.call('date 1>&2', shell=True)


# This function determines which kind of node the analysis is running on

def determine_node():
    mem = int(subprocess.check_output('free -g | grep "Mem: " | tr -s " " | cut -f 2 -d " "', shell=True).rstrip())
    cpu = int(subprocess.check_output(["nproc"]))
    node = subprocess.check_output(["hostname"])
    if cpu > 41:
        return "coruscant"
    elif cpu == 40 and mem > 80:
        return "wgs"
    elif cpu == 40 and mem < 80:
        return "clinprod"
    elif cpu < 40:
        return "short"
    else:
        return node


# Function splits input into perfectly mapped reads with NH:1 and QMAP > 3 and secondary mapped
# with NH:>1 and QMAP < 3

def sam_split(samfile_in, out_perfect, out_secondary):
    with open(samfile_in, "r") as sam:
        with open(out_perfect, "w+") as perfect:
            with open(out_secondary, "w+") as secondary:
                it1 = 0
                it2 = 0
                it3 = 0
                perfs = 0
                secs = 0
                truvalerrors = 0
                for line in sam:
                    old_line = line.split("\t")
                    it3 += 1
                    if len(old_line) != 13:
                        logging.warning("Not 13 columns on line. Instead {} columns were found on line {}".format(str(len(old_line)),
                                                                                                                  str(it3)))
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
                            else:
                                truvalerrors += 1
                            continue
                logging.info("{} perfect. {} secondary. {} TypeErrors. {} ValueErrors of which {} are real errors".format(perfs, secs, it1, it2, truvalerrors))


if __name__ == "__main__":

    threads = 8
    memory = 2

    if determine_node() == "coruscant":
        threads = 50
        memory = 3
    elif determine_node() == "wgs":
        threads = 40
        memory = 2
    elif determine_node() == "clinprod":
        threads = 20
        memory = 2
    elif determine_node() == "short":
        threads = 8
        memory = 2
    else:
        logging.info("node: {} is being used, using default resources: 8t2G".format(determine_node()))

    errchk()
    infile = argv[1]
    outfile_perfect = argv[2]
    outfile_secondary = argv[3]
    loggloc = argv[4]
    logging.basicConfig(level=logging.DEBUG, filename=str(loggloc), filemode="w",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    fasta_index = "/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta.fai"
    bwa_index = fasta_index.rsplit('.', 1)[0]
    directory = path.dirname(outfile_secondary)
    logging.info("Inputs: Infile = {}, Outfile_prefect = {}, Outfile_secondary = {}, loggloc = {}".format(str(infile),
                                                                                                          str(outfile_perfect),
                                                                                                          str(outfile_secondary),
                                                                                                          str(loggloc)))

    # samtools_module = "module load samtools/1.3.1"
    samtools_path = "/apps/bio/apps/samtools/1.3.1/samtools"
    # bwa_module = "module load bwa/0.7.5a"
    bwa_path = "/apps/bio/local/apps/bwa/0.7.5a/bwa"
    logging.info("initial variables set")

    # Determine if the inut is bam or sam
    if ".bam" in infile:
        logging.info("Input file was: BAM")
        errchk()
        intermediary_sam = infile.replace(".bam", ".sam")
        intermediary_sorted = infile.replace(".bam", ".sorted.bam")
        print(str(infile))
        print(str(intermediary_sam))
        try:
            # used to be check_call

            subprocess.call("{0} sort {1} -n -@ {2} -m {3}G -T {4} -o {6} && {0} view {6} -o {5} -@ {2} -h".format(samtools_path,
                                                                                                                   infile,
                                                                                                                   threads,
                                                                                                                   memory, directory,
                                                                                                                   intermediary_sam,
                                                                                                                   intermediary_sorted),
                            shell=True)
            # remove intermidary bam file which is no longer needed
            subprocess.call(["rm", str(intermediary_sorted)])
            # Old way to fix the input bam. The above way is 2,7% faster than the one below.
            # subprocess.call(
            #     "{0} sort {1} -n -@ {2} -m {3}G -T {4} | {0} view - -o {5} -@ {2} -h".format(samtools_path, infile,
            #                                                                                  threads, memory,
            #                                                                                  directory,
            #                                                                                  intermediary_sam),
            #    shell=True)
            errchk()
        except subprocess.CalledProcessError:
            logging.warning("CALLEDPROCESSERROR in initial sort")
        except OSError:
            print("OSERROR (bam)")
            logging.warning("OSERROR in initial sort (bam)")
        logging.info("Initial conversion of BAM to SAM (name sorted) completed. Now starting SAM-file splitting")
        # sam_split(intermediary, outfile_perfect, outfile_secondary)
        sam_parse.sam_split_runner(intermediary_sam, directory + "/", 2500000, threads)
        subprocess.call(["rm", str(intermediary_sam)])
        logging.info("Splitting of the SAM-file completed")
    elif ".sam" in infile:
        logging.info("input file was: SAM")
        temp = infile.replace(".sam", "_tmp.sam")
        try:
            # subprocess.call("", shell=True)
            # used to be check_call
            subprocess.call(
                "{0} view {1} -@ {2} -ht {3} | samtools sort - -@ {2} -m 2G -n -T {4} | samtools view - -o {5} -@ {2} "
                "-h".format(samtools_path, infile, threads, fasta_index, directory, temp), shell=True)
            errchk()
        except subprocess.CalledProcessError:
            logging.warning("CALLEDPROCESERROR in initial sort (sam)")
        except OSError:
            logging.warning("OSERROR in initial sort (sam)")
        logging.info("Initial name sorting completed, starting SAM-file splitting at: {}".format(str(datetime.datetime.now())))
        errchk()
        # sam_split(temp, outfile_perfect, outfile_secondary)
        sam_parse.sam_split_runner(temp, directory + "/", 2500000, threads)
        errchk()
        logging.info("Splitting of the SAM-file completed at: {}".format(str(datetime.datetime.now())))
        subprocess.call(["rm", temp])

    perfsize = path.getsize(outfile_perfect)
    secsize = path.getsize(outfile_secondary)
    logging.info("After splitting: ByteSize of Perfect: {}. Bytesize of secondary: {}".format(perfsize, secsize))

    # Get original headers
    original_headers = directory + "/original_headers.sam"
    subprocess.call("{} view -H {} > {}".format(samtools_path, infile, original_headers), shell=True)
    logging.info("Headers extracted from original input file: {} to : {}".format(infile, original_headers))
    errchk()

    secondary_tmp = outfile_secondary + "_temp"

    # Remove original infile in order to save space
    subprocess.call(["rm", str(infile)])
    logging.info("{} removed".format(infile))

    # Reheader the secondary mapped sam file
    subprocess.call("{} view -ht {} {} > {}".format(samtools_path, fasta_index, outfile_secondary, secondary_tmp),
                    shell=True)
    logging.info("Temporary secondary SAM-file: {} re-headered into : {}".format(outfile_secondary, secondary_tmp))
    errchk()
    # Rename tempfile to original name
    subprocess.call(["mv", secondary_tmp, outfile_secondary])
    logging.info("SAM-file : {} overwritten by temp file {}, which no longer exists".format(outfile_secondary,
                                                                                            secondary_tmp))
    errchk()
    # Convert to fastq
    subprocess.call(
        "%s fastq %s > %s/reads_interleaved.fastq" % (samtools_path, outfile_secondary, directory),
        shell=True)
    logging.info("Converted secondary to fastq")
    errchk()
    # Run bwa
    subprocess.call("%s mem %s -p %s/reads_interleaved.fastq -t %s > %s/bwa_out.sam" % (bwa_path, bwa_index,
                                                                                        directory, threads, directory),
                    shell=True)
    logging.info("Secondary reads re-mapped")
    errchk()

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
    errchk()

    # Remove now unnecessary perfect.sam and secondary.sam
    subprocess.call(["rm", str(outfile_perfect)])
    subprocess.call(["rm", str(outfile_secondary)])
    logging.info("{} and {} removed".format(outfile_perfect, outfile_secondary))
    errchk()

    # Convert the merged sam file into bam
    subprocess.call("{} view -b -@ {} {} -o {}".format(samtools_path, threads, bwa_out, merged_bam), shell=True)
    logging.info("Merged SAM converted to BAM")
    errchk()

    # Remove large merged sam file
    subprocess.call(["rm", str(bwa_out)])
    logging.info("{} removed".format(bwa_out))

    # Reheader the merged bam file
    rehead_tmp = directory + "/temporary_reheader.bam"
    subprocess.call("%s reheader %s %s > %s" % (samtools_path, original_headers, merged_bam, rehead_tmp),
                    shell=True)
    logging.info("BAM file: {} reheadered using: {} into tempfile: {}".format(merged_bam, original_headers, rehead_tmp))
    errchk()
    subprocess.call(["mv", rehead_tmp, merged_bam])
    logging.info("Non-headered bamfile: {} overwritten by headered tempfile: {}".format(merged_bam, rehead_tmp))
    errchk()

    # Sort reheadered bam file
    # with open(sorted_bam, "wb") as sorted:
    subprocess.call("/apps/bio/apps/samtools/1.3.1/samtools sort -T %s -@ %s -m %sG %s > %s" % (directory, threads,
                                                                                                memory, merged_bam,
                                                                                                sorted_bam),
                    shell=True)
    logging.info("Final BAM: {} sorted to: {}".format(merged_bam, sorted_bam))
    errchk()

    # Removing merged.bam to save space
    subprocess.call(["rm", merged_bam])
    logging.info("{} removed".format(merged_bam))

    subprocess.call(["/apps/bio/apps/samtools/1.3.1/samtools", "index", sorted_bam])
    logging.info("Final BAM: {} indexed into: {}".format(sorted_bam, sorted_bam + ".bai"))
    errchk()
    # Remove all intermediary files. Fill in later when testing is done.
    # This is currently being done in the wrapper script so I will leave this commented
    # subprocess.call(["rm", ])

    # Give path to result
    print("Location of output file: \n" + str(path.abspath(sorted_bam)))
    logging.info("Everything is Completed at {}.".format(str(datetime.datetime.now())))
    errchk()
