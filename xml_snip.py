#!/usr/bin/env python

from sys import argv
import os


def igv_modification(user, infile):
    with open("/Users/alvaralmstedt/PycharmProjects/CLC_External_Apps/test_folder/alvar.xml", "r+") as userfile:

        lines_of_file = userfile.readlines()
        bam = os.path.basename(infile)
        lines_of_file.insert(-2, '\t\t<Resource name="%s" path="http://medstore.sahlgrenska.gu.se:8008/data/%s/%s" />\n' % (bam, user, bam))
        userfile.seek(0)
        userfile.truncate()
        userfile.writelines(lines_of_file)

        # newfile = []
        # for line in userfile.readlines()[:-2]:
            # if "<Resource name=" in line:
            # newfile.append(line.rstrip("\n"))
        # newfile.append('\t\t<Resource name="%s" path="http://medstore.sahlgrenska.gu.se:8008/data/%s/%s" />' % (bam,
        #                                                                                                        user,
        #                                                                                                        bam))
        # newfile.append("\t</Category>")
        # newfile.append("</Global>")
        # for j in newfile:
            # userfile.write(j)

# infile = argv[1]
infile = "/Users/alvaralmstedt/PycharmProjects/CLC_External_Apps/test_folder/alvar.xml"
# user = argv[2]
user = "alvar.almstedt"

igv_modification(user, infile)