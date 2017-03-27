#!/usr/bin/env python

import os
from subprocess import call
from sys import argv
from math import log
import socket

bam_file = argv[1]
mode = argv[2]
#vcf_in_1 = argv[3]
#vcf_in_2 = argv[4]
vcf_out = argv[3]
cnv_text = argv[4]
cnv_copynumber = argv[5]

error_file = open("/tmp/canvaserror.log", 'w+')
error_file.write(str(socket.gethostname()))

for i in argv:
    print(i)
    error_file.write(i + "\n")

if "Somatic" in mode:
    mode = "Somatic-WGS"
else:
    mode = "Germline-WGS"

call("mkdir /tmp/canvas", shell=True)
call("mkdir /tmp/canvas/bam", shell=True)
call("mkdir /tmp/canvas/outdir", shell=True)

indexed = 0

array = bam_file.split("/")
filename = array[-1]
bam_path = "/tmp/canvas/bam/%s" % filename

error_file.write("Bampath before looking in text: %s" % bam_path)
error_file.write("\n")

if ".txt" in str(bam_file):
    bam_temp = bam_file.split("/")
    bam_file = "/" + str("/".join(bam_temp[1:]))
    bam_text_file = open(bam_file, "r")
    bam_file = bam_text_file.readline().rstrip()
    bam_text_file.close()
    array = bam_file.split("/")
    filename = array[-1]
#    bam_path = "/tmp/canvas/bam/%s" % filename.rstrip()
    bam_path = "/tmp/canvas/bam/"
    indexed = 1
    print("bam_file: " + str(bam_file))
    error_file.write("Bamfile after looking in text: %s\n" % bam_file)

error_file.write("Bampath after looking in text: %s" % bam_path)

call(["cp", str(bam_file), "-t", str(bam_path)])
call(["cp", str(bam_file) + ".bai", "-t", str(bam_path)])
call("cp -r /medstore/External_References/Canvas_CLC_HG19_Dataset /tmp/canvas/", shell=True)
call("cp -r /medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta* /tmp/canvas/Canvas_CLC_HG19_Dataset", shell=True)

if not indexed:
    call("module load samtools/1.3.1", shell=True)
    call("/medstore/IGV_Folders/samtools index /tmp/canvas/bam/%s" % filename, shell=True)

call(["/usr/bin/mono", "/apps/CLC_ExternalApps/canvas/1.11.0/Canvas.exe", str(mode), "-b",
      "/tmp/canvas/bam/" + str(filename), "--b-allele-vcf=/tmp/canvas/Canvas_CLC_HG19_Dataset/dbsnp_common_all_20160601.vcf",
      "-o", "/tmp/canvas/outdir", "--reference=/tmp/canvas/Canvas_CLC_HG19_Dataset/kmer.fa",
      "-g", "/tmp/canvas/Canvas_CLC_HG19_Dataset/", "-f", "/tmp/canvas/Canvas_CLC_HG19_Dataset/filter13.bed",
      "-n", "WGS", "--custom-parameters=CanvasBin,-p"])


with open("/tmp/canvas/outdir/CNV.CoverageAndVariantFrequency.txt", "r") as INFILE:
    with open("/tmp/canvas/outdir/CNV_observed.seg", "w+") as OUTFILE:
        OUTFILE.write("#track graphType=points maxHeightPixels=300:300:300 color=0,220,0 altColor=220,0,0\n")
        OUTFILE.write("Sample\tChromosome\tStart\tEnd\tCNV_Observed\n")
        with open("/tmp/canvas/outdir/CNV_called.seg", "w+") as OUTFILE2:
            OUTFILE2.write("#track graphType=points maxHeightPixels=300:300:300 color=0,220,0 altColor=220,0,0\n")
            OUTFILE2.write("Sample\tChromosome\tStart\tEnd\tCNV_Called\n")
            for line in INFILE:
                array_2 = line.split("\t")
                length = len(array_2)
                cnv = 0
                ncov = 0
                try:
                    cnv = int(array_2[3])
                    ncov = int(array_2[6])
                except ValueError:
                    try:
                        cnv = float(array_2[3])
                        ncov = float(array_2[6])
                    except ValueError:
                        try:
                            print(str(array_2[3]) + " Cannot be converted to an int or float")
                            print(str(array_2[6]) + " Cannot be converted to an int or float")
                            continue
                        except IndexError:
                            continue
                if ncov > 0 and cnv > 0:
                    cnvlog = log(cnv, 10) / log(2)
                    covlog = log(ncov, 10) / log(2)
                    if not array_2[0] == "X" or not array_2[0] == "Y":
                        cnvlog = cnvlog - 1
                        covlog = covlog - 1

                    OUTFILE.write("Observed_CNVs\t%s\t%s\t%s\t%s" % (array_2[0], array_2[1], array_2[2], covlog))
                    OUTFILE2.write("Called_CNVs\t%s\t%s\t%s\t%s" % (array_2[0], array_2[1], array_2[2], cnvlog))

call("gunzip /tmp/canvas/outdir/CNV.vcf.gz", shell=True)
call("mv /tmp/canvas/outdir/CNV.vcf %s" % vcf_out, shell=True)
call("mv /tmp/canvas/outdir/CNV_observed.seg %s" % cnv_copynumber, shell=True)
call("mv /tmp/canvas/outdir/CNV_called.seg %s" % cnv_copynumber, shell=True)
call("mv /tmp/canvas/outdir/CNV.CoverageAndVariantFrequency.txt %s" % cnv_text, shell=True)
call("rm -rf /tmp/canvas", shell=True)