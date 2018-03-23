#!/usr/bin/env python

import os
from subprocess import call
from sys import argv
from math import log
import socket

with open("/tmp/canvas/outdir/CNV.CoverageAndVariantFrequency.txt", "r") as INFILE:
    with open("/tmp/canvas/outdir/CNV_log2.cn", "w+") as OUTFILE:
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
                        print(str(array_2[3]) + " Cannot be converted to an into or float")
                        print(str(array_2[6]) + " Cannot be converted to an in or a float")
                        continue
                    except IndexError:
                        continue
            if ncov > 0 and cnv > 0:
                cnvlog = log(cnv, 10) / log(2)
                covlog = log(ncov, 10) / log(2)
                if not array_2[0] == "X" or not array_2[0] == "Y":
                    cnvlog = cnvlog - 1
                    covlog = covlog - 1

                OUTFILE.write("WGS\t%s\t%s\t%s\t%s\n" % (array_2[0], array_2[1], cnvlog, covlog))
