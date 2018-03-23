#!/usr/bin/env python

import os
from subprocess import call
from sys import argv
from math import log
import socket
import datetime
import argparse
<<<<<<< HEAD
import glob
import shutil
=======
>>>>>>> 80c8c242ab5a15af3e8012b3d55375fbc93e8d23

timestring = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")
# bam_file = argv[1]
# mode = argv[2]
# # vcf_in_1 = argv[3]
# # vcf_in_2 = argv[4]
# vcf_out = argv[3]
# cnv_text = argv[4]
# cnv_copynumber_obs = argv[5]
# cnv_copynumber_call = argv[6]
# uname = argv[7]
# custom_uname = argv[8]
# manifest = argv[9]
# normal_bam = argv[10]
#
# print argv


parser = argparse.ArgumentParser()

parser.add_argument("-b", "--bam_file", nargs="?", type=str, action='store', help='Full path to your bam input file')
parser.add_argument("-m", "--mode", nargs="?", type=str, action='store', help='Specify which mode canvas will run in')
parser.add_argument("-v", "--vcf_out", nargs="?", action='store', type=str, help='Specify the path where'
                                                                                 ' the vcf file will be sent')
parser.add_argument("-t", "--cnv_text", nargs="?", action='store', type=str, help='Specify the path where'
                                                                                  ' the CNV textfile file will be sent')
parser.add_argument("-o", "--cnv_obs", nargs="?", action='store', type=str, help='Specify the path where the'
                                                                                 ' CNV_observed.seg file will be sent')
parser.add_argument("-c", "--cnv_call", nargs="?", action='store', type=str, help='Specify the path where the'
                                                                                  ' CNV_called.seg file will be sent')
parser.add_argument("-u", "--uname", nargs="?", action='store', type=str, help='Default username taken from CLC')
parser.add_argument("-n", "--custom_uname", nargs="?", action='store', type=str, help='Selected IGV username.'
                                                                                      ' Overrides default')
parser.add_argument("-a", "--manifest", nargs="?", action='store', type=str, help='Specify the path to the exome'
                                                                                  ' manifest file')
parser.add_argument("-r", "--normal_bam", nargs="?", action='store', type=str, help='Full path to bam_normal input file')
<<<<<<< HEAD
parser.add_argument("-s", "--sex", nargs="?", action='store', type=str, help='Specify sex of input sample')
=======
>>>>>>> 80c8c242ab5a15af3e8012b3d55375fbc93e8d23

args = parser.parse_args()

bam_file = str(args.bam_file)
print(bam_file)
mode = str(args.mode)
print(mode)
vcf_out = str(args.vcf_out)
print(vcf_out)
cnv_text = str(args.cnv_text)
print(cnv_text)
cnv_copynumber_obs = str(args.cnv_obs)
print(cnv_copynumber_obs)
cnv_copynumber_call = str(args.cnv_call)
print(cnv_copynumber_call)
uname = str(args.uname)
print(uname)
custom_uname = str(args.custom_uname)
print(custom_uname)
manifest = str(args.manifest)
if manifest:
    print(manifest)
normal_bam = str(args.normal_bam)
if normal_bam:
    print(normal_bam)
<<<<<<< HEAD
sex = str(args.sex)
print(sex)
=======
>>>>>>> 80c8c242ab5a15af3e8012b3d55375fbc93e8d23
print args

igv_data_folder = "/medstore/IGV_Folders/igv/data/%s" % uname


def igv_modification(user, infile):
    """
    Adds entries for seg files to user IGV xml list.
    """
    with open("/medstore/IGV_Folders/igv/users/%s_igv.xml" % user, "r+") as userfile:
        lines_of_file = userfile.readlines()
        bam = os.path.basename(infile)
        lines_of_file.insert(-2,
                             '\t\t<Resource name="%s" path="http://medstore.sahlgrenska.gu.se:8008/data/%s/%s" />\n' % (
                                 bam, user, bam))
        userfile.seek(0)
        userfile.truncate()
        userfile.writelines(lines_of_file)

error_file = open("/tmp/canvaserror.log", 'w+')
error_file.write(str(socket.gethostname()))

# for i in argv:
#    print(i)
#    error_file.write(i + "\n")

# if "Somatic-WGS" in mode:
#     mode = "Somatic-WGS"
# elif "Enrichment" in mode:
#     mode = "Somatic-Enrichment"
# elif "enrichment" in mode:
#     mode = "Tumor-normal-enrichment"
# else:
#     mode = "Germline-WGS"

shutil.rmtree("/tmp/canvas_dir")
call("mkdir /tmp/canvas_dir", shell=True)
call("mkdir /tmp/canvas_dir/bam", shell=True)
call("mkdir /tmp/canvas_dir/outdir", shell=True)
call("hostname")

indexed_bam = 0
indexed_normal = 0

array = bam_file.split("/")
bam_filename = array[-1]
bam_path = "/tmp/canvas_dir/bam/%s" % bam_filename

normal_array = normal_bam.split("/")
normal_filename = normal_array[-1]
normal_path = "/tmp/canvas_dir/bam/%s" % normal_filename

error_file.write("Bampath before looking in text: %s" % bam_path)
error_file.write("\n")

if str(bam_file).endswith(".txt"):
    print("Bam in textfile")
    bam_temp = bam_file.split("/")
    bam_file = "/" + str("/".join(bam_temp[1:]))
    bam_text_file = open(bam_file, "r")
    bam_file = bam_text_file.readline().rstrip()
    bam_text_file.close()
    array = bam_file.split("/")
    bam_filename = array[-1]
    #    bam_path = "/tmp/canvas/bam/%s" % filename.rstrip()
    bam_path = "/tmp/canvas_dir/bam/"
    indexed_bam = 1
    print("bam_file: " + str(bam_file))
    error_file.write("Bamfile after looking in text: %s\n" % bam_file)

if str(normal_bam).endswith(".txt"):
    print("Normal in textfile")
    normal_temp = normal_bam.split("/")
    normal_bam = "/" + str("/".join(normal_temp[1:]))
    normal_text_file = open(normal_bam, "r")
    normal_bam = normal_text_file.readline().rstrip()
    normal_text_file.close()
    array_n = normal_bam.split("/")
    normal_filename = array_n[-1]
    #    bam_path = "/tmp/canvas/bam/%s" % filename.rstrip()
    bam_path = "/tmp/canvas_dir/bam/"
    indexed_normal = 1
    print("normal_bam_file: " + str(bam_file))
    error_file.write("Normal bamfile after looking in text: %s\n" % normal_bam)

error_file.write("Bampath after looking in text: %s" % bam_path)

# Copy the bams either from CLC tmp or from IGV folder
call(["cp", str(bam_file), "-t", "/tmp/canvas_dir/bam"])
if normal_bam:
    call(["cp", str(normal_bam), "-t", "/tmp/canvas_dir/bam"])

# call(["cp", str(manifest), "-t", str(manifest_path)])

# Will only be run on text bams
if indexed_bam:
    print("bam.bai being copied from %s to %s" % (bam_file, bam_path))
    call(["cp", str(bam_file) + ".bai", "-t", str(bam_path)])
if indexed_normal and normal_bam:
    print("normal.bai being copied from %s to %s" % (bam_file, bam_path))
    call(["cp", str(normal_bam) + ".bai", "-t", str(bam_path)])


call("cp -r /medstore/External_References/Canvas_CLC_HG19_Dataset /tmp/canvas_dir/", shell=True)
call(
    "cp /medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta /tmp/canvas_dir/Canvas_CLC_HG19_Dataset",
    shell=True)
call("cp /medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta.fai /tmp/canvas_dir/Canvas_CLC_HG19_Dataset", shell=True)

# Will only be run on non-text bam files.
if not indexed_bam:
    # call("module load samtools/1.3.1", shell=True)
    call("/medstore/IGV_Folders/samtools index /tmp/canvas_dir/bam/%s" % bam_filename, shell=True)
<<<<<<< HEAD

if not indexed_normal:
    # call("module load samtools/1.3.1", shell=True)
    call("/medstore/IGV_Folders/samtools index /tmp/canvas_dir/bam/%s" % normal_filename, shell=True)

=======

if not indexed_normal:
    # call("module load samtools/1.3.1", shell=True)
    call("/medstore/IGV_Folders/samtools index /tmp/canvas_dir/bam/%s" % normal_filename, shell=True)

>>>>>>> 80c8c242ab5a15af3e8012b3d55375fbc93e8d23
if "Somatic-WGS" in mode:
    mode = "Somatic-WGS"
    print("Somatic-WGS selected")
    call(["/usr/bin/mono", "/apps/CLC_ExternalApps/canvas/1.11.0/Canvas.exe", str(mode), "-b",
          "/tmp/canvas_dir/bam/" + str(bam_filename),
          "--b-allele-vcf=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/dbsnp_common_all_20160601.vcf",
          "--exclude-non-het-b-allele-sites",
          "-o", "/tmp/canvas_dir/outdir", "--reference=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/kmer.fa",
          "-g", "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/", "-f", "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/filter13.bed",
          "-n", "WGS", "--custom-parameters=CanvasBin,-p"])
elif "Enrichment" in mode:
    mode = "Somatic-Enrichment"
    print("Somatic-Enrichment selected")
    call(["cp", str(manifest), "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/manifest.txt"])
<<<<<<< HEAD
    if sex == "male":
        print("male")
        pattern = "/medstore/CLC_Import_Export/Alvar_Almstedt/canvas_related/control_samples/binned/male/*/*"
        paths = glob.glob(pattern)
    else:
        print("female")
=======
>>>>>>> 80c8c242ab5a15af3e8012b3d55375fbc93e8d23
    call(["/usr/bin/mono", "/apps/CLC_ExternalApps/canvas/1.11.0/Canvas.exe", str(mode), "-b",
          "/tmp/canvas_dir/bam/" + str(bam_filename),
          "--b-allele-vcf=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/dbsnp_common_all_20160601.vcf",
          "--exclude-non-het-b-allele-sites",
          "-o", "/tmp/canvas_dir/outdir", "--reference=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/kmer.fa",
          "--manifest", "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/manifest.txt",
          "-g", "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/", "-f",
          "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/filter13.bed",
          "-n", "WGS", "--custom-parameters=CanvasBin,-p"])
elif "enrichment" in mode:
    mode = "Tumor-normal-enrichment"
    print("Tumor-normal-enrichment selected")
    call(["cp", normal_bam, "-t", "/tmp/canvas_dir/bam"])
    call(["cp", str(manifest), "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/manifest.txt"])
    call(["/usr/bin/mono", "/apps/CLC_ExternalApps/canvas/1.11.0/Canvas.exe", str(mode), "-b",
          "/tmp/canvas_dir/bam/" + str(bam_filename),
          "--normal-bam", str(normal_path),
          "--b-allele-vcf=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/dbsnp_common_all_20160601.vcf",
          "--exclude-non-het-b-allele-sites",
          "-o", "/tmp/canvas_dir/outdir", "--reference=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/kmer.fa",
          "--manifest", "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/manifest.txt",
          "-g", "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/", "-f",
          "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/filter13.bed",
          "-n", "WGS", "--custom-parameters=CanvasBin,-p"])
else:
    mode = "Germline-WGS"
    print("Germline-WGS selected")
    call(["/usr/bin/mono", "/apps/CLC_ExternalApps/canvas/1.11.0/Canvas.exe", str(mode), "-b",
          "/tmp/canvas_dir/bam/" + str(bam_filename),
          "--b-allele-vcf=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/dbsnp_common_all_20160601.vcf",
          "--exclude-non-het-b-allele-sites",
          "-o", "/tmp/canvas_dir/outdir", "--reference=/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/kmer.fa",
          "-g", "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/", "-f",
          "/tmp/canvas_dir/Canvas_CLC_HG19_Dataset/filter13.bed",
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
                    print("Exception: ", e)
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
call("cp /tmp/canvas_dir/outdir/CNV_observed.seg %s/CNV_observed_%s.seg" % (igv_data_folder, timestring), shell=True)
call("mv /tmp/canvas_dir/outdir/CNV_observed.seg %s" % cnv_copynumber_obs, shell=True)
call("cp /tmp/canvas_dir/outdir/CNV_called.seg %s/CNV_called_%s.seg" % (igv_data_folder, timestring), shell=True)
call("mv /tmp/canvas_dir/outdir/CNV_called.seg %s" % cnv_copynumber_call, shell=True)
call("mv /tmp/canvas_dir/outdir/CNV.CoverageAndVariantFrequency.txt %s" % cnv_text, shell=True)

if os.path.isfile("/medstore/IGV_Folders/igv/users/{}_igv.xml".format(custom_uname)):
    igv_modification(custom_uname, igv_data_folder + "/CNV_observed_{}.seg".format(timestring))
    igv_modification(custom_uname, igv_data_folder + "/CNV_called_{}.seg".format(timestring))
else:
    print("{} is not a valid user. IGV destination")
    igv_modification(uname, igv_data_folder + "/CNV_observed_{}.seg".format(timestring))
    igv_modification(uname, igv_data_folder + "/CNV_called_{}.seg".format(timestring))

# call("rm -rf /tmp/canvas_dir", shell=True)
