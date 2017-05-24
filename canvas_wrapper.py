#!/usr/bin/env python

import os
from subprocess import call
from sys import argv
from math import log
import socket

bam_file = argv[1]
mode = argv[2]
# vcf_in_1 = argv[3]
# vcf_in_2 = argv[4]
vcf_out = argv[3]
cnv_text = argv[4]
cnv_copynumber_obs = argv[5]
cnv_copynumber_call = argv[6]
uname = argv[7]

igv_data_folder = "/medstore/IGV_Folders/igv/data/%s" % uname


def igv_modification(user, infile):
    with open("/medstore/IGV_Folders/igv/users/%s_igv.xml" % user, "r+") as userfile:
        lines_of_file = userfile.readlines()
        bam = os.path.basename(infile)
        lines_of_file.insert(-2,
                             '\t\t<Resource name="%s" path="http://medstore.sahlgrenska.gu.se:8008/data/%s/%s" />\n' % (
                                 bam, user, bam))
        userfile.seek(0)
        userfile.truncate()
        userfile.writelines(lines_of_file)

        # Old solution (not the best)
        # bam = os.path.basename(infile)
        # newfile = []
        # for line in userfile.readlines()[:-2]:
        #     # if "<Resource name=" in line:
        #     newfile.append(line.rstrip("\n"))
        # newfile.append('\t\t<Resource name="%s" path="http://medstore.sahlgrenska.gu.se:8008/data/%s/%s" />' % (bam,
        #                                                                                                        user,
        #                                                                                                        bam))
        # newfile.append("\t</Category>")
        # newfile.append("</Global>")
        # for j in newfile:
        #     userfile.write(j)


error_file = open("/tmp/canvaserror.log", 'w+')
error_file.write(str(socket.gethostname()))

for i in argv:
    print(i)
    error_file.write(i + "\n")

if "Somatic" in mode:
    mode = "Somatic-WGS"
else:
    mode = "Germline-WGS"

call("mkdir /tmp/canvas_dir", shell=True)
call("mkdir /tmp/canvas_dir/bam", shell=True)
call("mkdir /tmp/canvas_dir/outdir", shell=True)
call("hostname")

indexed = 0

array = bam_file.split("/")
filename = array[-1]
bam_path = "/tmp/canvas_dir/bam/%s" % filename

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
    bam_path = "/tmp/canvas_dir/bam/"
    indexed = 1
    print("bam_file: " + str(bam_file))
    error_file.write("Bamfile after looking in text: %s\n" % bam_file)

error_file.write("Bampath after looking in text: %s" % bam_path)

call(["cp", str(bam_file), "-t", str(bam_path)])
call(["cp", str(bam_file) + ".bai", "-t", str(bam_path)])
call("cp -r /medstore/External_References/Canvas_CLC_HG19_Dataset /tmp/canvas_dir/", shell=True)
call(
    "cp -r /medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta* /tmp/canvas_dir/Canvas_CLC_HG19_Dataset",
    shell=True)

if not indexed:
    call("module load samtools/1.3.1", shell=True)
    call("/medstore/IGV_Folders/samtools index /tmp/canvas_dir/bam/%s" % filename, shell=True)

call(["/usr/bin/mono", "/apps/CLC_ExternalApps/canvas/1.11.0/Canvas.exe", str(mode), "-b",
      "/tmp/canvas_dir/bam/" + str(filename),
      "--b-allele-vcf=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/dbsnp_common_all_20160601.vcf",
      "-o", "/tmp/canvas_dir/outdir", "--reference=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/kmer.fa",
      "-g", "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/", "-f", "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/filter13.bed",
      "-n", "WGS", "--custom-parameters=CanvasBin,-p"])

with open("/tmp/canvas_dir/outdir/CNV.CoverageAndVariantFrequency.txt", "r") as INFILE:
    with open("/tmp/canvas_dir/outdir/CNV_observed.seg", "w+") as OUTFILE:
        OUTFILE.write("#track graphType=points maxHeightPixels=300:300:300 color=0,220,0 altColor=220,0,0\n")
        OUTFILE.write("Sample\tChromosome\tStart\tEnd\tCNV_Observed\n")
        with open("/tmp/canvas_dir/outdir/CNV_called.seg", "w+") as OUTFILE2:
            OUTFILE2.write("#track graphType=points maxHeightPixels=300:300:300 color=0,220,0 altColor=220,0,0\n")
            OUTFILE2.write("Sample\tChromosome\tStart\tEnd\tCNV_Called\n")
            for line in INFILE:
                array_2 = line.split("\t")
                print("array_2: ", array_2)
                length = len(array_2)
                cnv = 0
                ncov = 0

                try:
                    cnv = float(array_2[3])
                    ncov = float(array_2[6])
                    print("cnv: ", cnv)
                    print("ncov: ", ncov)
                except (IndexError, ValueError) as e:
                    print(e)
                    continue

                if ncov > 0 and cnv > 0:
                    print("PASSED")
                    cnvlog = log(cnv, 2)
                    covlog = log(ncov, 2)
                    if not array_2[0] == "X" or not array_2[0] == "Y":
                        cnvlog -= 1
                        covlog -= 1

                    OUTFILE.write("Observed_CNVs\t%s\t%s\t%s\t%s\n" % (array_2[0], array_2[1], array_2[2], covlog))
                    OUTFILE2.write("Called_CNVs\t%s\t%s\t%s\t%s\n" % (array_2[0], array_2[1], array_2[2], cnvlog))

call("gunzip /tmp/canvas_dir/outdir/CNV.vcf.gz", shell=True)
call("mv /tmp/canvas_dir/outdir/CNV.vcf %s" % vcf_out, shell=True)
call("cp /tmp/canvas_dir/outdir/CNV_observed.seg %s" % igv_data_folder, shell=True)
call("mv /tmp/canvas_dir/outdir/CNV_observed.seg %s" % cnv_copynumber_obs, shell=True)
call("cp /tmp/canvas_dir/outdir/CNV_called.seg %s" % igv_data_folder, shell=True)
call("mv /tmp/canvas_dir/outdir/CNV_called.seg %s" % cnv_copynumber_call, shell=True)
call("mv /tmp/canvas_dir/outdir/CNV.CoverageAndVariantFrequency.txt %s" % cnv_text, shell=True)

igv_modification(uname, igv_data_folder + "/CNV_observed.seg")
igv_modification(uname, igv_data_folder + "/CNV_called.seg")

call("rm -rf /tmp/canvas_dir", shell=True)
