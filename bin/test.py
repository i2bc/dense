#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 08:44:36 2023

@author: paul.roginski
"""

import sys

#2.7 sec with Mmus.gff
#with open(sys.argv[1], 'r') as gff:
#    for line in gff:
#        if not line.startswith("#"):
#            if line.split('\t')[2] == "CDS":
#                print(line.strip())
   
#30 sec!             
from pybedtools.bedtool import BedTool,Interval
#HOME = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/" 
#gff_filename = HOME+"Eukaryotes/MODELS/SCER/GENOMES/Scer_NCBI.gff"
#gff = BedTool(gff_filename)
#gff = BedTool(sys.argv[1])
#gff_CDS=BedTool([interval for interval in gff if interval.fields[2] == "CDS"])
#print(len(gff_CDS))
#flanked=gff.flank(g=sys.argv[2], b=75)
#print(flanked)

#1.5 sec or less
#import subprocess
##HOME = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/" 
##gff_path = HOME+"Eukaryotes/MODELS/SCER/GENOMES/Scer_NCBI.gff"
#gff_path = sys.argv[1]
#gff_CDS=subprocess.check_output(["grep", "CDS", gff_path])
#gff_CDS = gff_CDS.decode("ascii",errors="ignore")
##for line in gff_CDS.splitlines():
##    print(line)
##    break
#test = gff_CDS.splitlines()


import subprocess
flank=subprocess.check_output(["bedtools", "flank", "-i", sys.argv[1], "-g", sys.argv[2], "-s", '-b', str(75)])
flank = flank.decode("ascii",errors="ignore")
for i in flank.splitlines():
    print(i)