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

for i in str(argv):
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
bam_path = "/tmp/canvas/%s" % filename

error_file.write("Bampath before looking in text: %s" % bam_path)
error_file.write("\n")

if ".txt" in str(bam_file):
    bam_temp = bam_file.split("/")
    bam_file = "/" + str("/".join(bam_temp[1:]))
    bam_text_file = open(bam_file, "r")
    bam_file = bam_text_file.readline()
    bam_text_file.close()
    array = bam_file.split("/")
    filename = array[-1]
    bam_path = "/tmp/canvas/%s" % filename
    indexed = 1

error_file.write("Bampath after looking in text: %s" % bam_path)

call("cp %s* /tmp/canvas/bam/" % bam_file, shell=True)
call("cp -r /medstore/External_References/Canvas_CLC_HG19_Dataset /tmp/canvas/", shell=True)
call("cp -r /medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta* /tmp/canvas/Canvas_CLC_HG19_Dataset", shell=True)

if not indexed:
    call("module load samtools/1.3.1", shell=True)
    call("/medstore/IGV_Folders/samtools index /tmp/canvas/bam/%s" % filename, shell=True)

command = """/usr/bin/mono /apps/CLC_ExternalApps/canvas/1.11.0/Canvas.exe %s -b
/tmp/canvas/bam/%s
--b-allele-vcf=/tmp/canvas/Canvas_CLC_HG19_Dataset/dbsnp_common_all_20160601.vcf
--exclude-non-het-b-allele-sites
-o /tmp/canvas/outdir
--reference=/tmp/canvas/Canvas_CLC_HG19_Dataset/kmer.fa
-g /tmp/canvas/Canvas_CLC_HG19_Dataset/
-f /tmp/canvas/Canvas_CLC_HG19_Dataset/filter13.bed
-n WGS
--custom-parameters=CanvasBin,-p""" % (mode, filename)

call(command, shell=True)

with open("/tmp/canvas/outdir/CNV.CoverageAndVariantFrequency.txt", "r") as INFILE:
    with open("/tmp/canvas/outdir/CNV_log2.cn", "w+") as OUTFILE:
        for line in INFILE:
            array_2 = line.split("\t")
            length = len(array_2)
            cnv = 0
            ncov = 0
            cnv = int(array_2[3])
            ncov = int(array_2[6])

            if ncov > 0 and cnv > 0:
                cnvlog = log(cnv, 10) / log(2)
                covlog = log(ncov, 10) / log(2)
                if not array_2[0] == "X" or not array_2[0] == "Y":
                    cnvlog = cnvlog - 1
                    covlog = covlog - 1

                OUTFILE.write("WGS\t%s\t%s\t%s\t%s" % (array_2[0], array_2[1], cnvlog, covlog))


call("gunzip /tmp/canvas/outdir/CNV.vcf.gz", shell=True)
call("mv /tmp/canvas/outdir/CNV.vcf %s" % vcf_out, shell=True)
call("mv /tmp/canvas/outdir/CNV_log2.cn %s" % cnv_copynumber, shell=True)
call("mv /tmp/canvas/outdir/CNV.CoverageAndVariantFrequency.txt %s" % cnv_text, shell=True)
call("rm -rf /tmp/canvas", shell=True)
