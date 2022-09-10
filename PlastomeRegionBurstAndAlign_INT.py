#!/usr/bin/env python3
'''Extract And Align Introns Across Multiple Plastomes'''
__version__ = 'm.gruenstaeudl@fu-berlin.de|2022-09-09T16:22:47 CEST'

#-----------------------------------------------------------------#
## IMPORTS
import argparse
import Bio
from Bio import SeqIO  # line is necessary for similar reason stated in: https://www.biostars.org/p/13099/
from Bio.Align import Applications  # line is necessary for similar reason stated in: https://www.biostars.org/p/13099/
import collections
import copy
import os
import pdb
import sys

#-----------------------------------------------------------------#
# DEBUGGING HELP
import pdb
#pdb.set_trace()

#-----------------------------------------------------------------#
# FUNCTIONS
def extract_INT(rec, feature, gene_name, offset):
            try:
                    feature.location = Bio.SeqFeature.FeatureLocation(feature.location.parts[offset].end, feature.location.parts[offset+1].start)
            except:
                    feature.location = Bio.SeqFeature.FeatureLocation(feature.location.parts[offset+1].start, feature.location.parts[offset].end)
            try:
                seq_name = gene_name + '_' + rec.name
                seq_obj = feature.extract(rec).seq   # Here the actual extraction is conducted
                seq_rec = Bio.SeqRecord.SeqRecord(seq_obj, id=seq_name, name='', description='')
                return(seq_rec,gene_name)
            except:
                raise Exception("  ERROR: Cannot conduct intron extraction for %s" % feature.qualifiers['gene'])
