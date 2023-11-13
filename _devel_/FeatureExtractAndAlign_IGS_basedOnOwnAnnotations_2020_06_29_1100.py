#!/usr/bin/env python2.7

''' Extract And Align Features Regions Across Multiple Plastomes '''

#####################
# IMPORT OPERATIONS #
#####################

import Bio
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqRecord import SeqRecord

from Bio.Align.Applications import MafftCommandline

from Bio.Seq import Seq
from Bio import AlignIO  # For function 'AlignIO.convert'
from StringIO import StringIO  # For file handles
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Nexus import Nexus  # For functions 'Nexus.combine' and 'Nexus.Nexus'

import collections  # For 'collections.OrderedDict()'
#import subprocess
import os
import re
import sys

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2020 Michael Gruenstaeudl'
__info__ = 'Extract And Align Features Regions Across Multiple Plastomes'
__version__ = '2020.06.29.1100'

#############
# DEBUGGING #
#############

import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

###########
# CLASSES #
###########

#############
# FUNCTIONS #
#############

########################################################################
# The concept of extracting introns based on the location of the flanking 
# exons was abandoned because the occurrence of incorrect extractions was 
# too high. For example, there are cases (e.g., the trnK intron), where 
# the situtation is far from clear. Hence, it is easier to extract based 
# on annotations only.

'''
def extract_collect_features(masterdict_nucl, inFn, fmt='genbank'):
    rec = SeqIO.read(inFn, fmt)
    
    # EXTRACT ALL GENES FROM RECORD (I.E., CDS, TRNA, RRNA)
    all_introns = [feature for feature in rec.features if feature.type=='intron']  # Resulting list contains adjacent introns in order of appearance on genome

    # LOOP THROUGH GENES
    for idx in range(0, len(all_introns)-1):
        current_feature = all_introns[idx]
        adjacent_feature = all_introns[idx+1]

        # Define names of Introns
        if 'gene' in current_feature.qualifiers and 'gene' in adjacent_feature.qualifiers:
            current_feature_name = current_feature.qualifiers['gene'][0]
            current_feature_name_SAFE = current_feature_name.replace('-','_')
            current_feature_name_SAFE = re.sub(r'[^\w]', '', current_feature_name_SAFE)
            adjacent_feature_name = adjacent_feature.qualifiers['gene'][0]
            adjacent_feature_name_SAFE = adjacent_feature_name.replace('-','_')
            adjacent_feature_name_SAFE = re.sub(r'[^\w]', '', adjacent_feature_name_SAFE)
            feature_name = current_feature_name_SAFE + '_' + adjacent_feature_name_SAFE
            inv_feature_name = adjacent_feature_name_SAFE + '_' + current_feature_name_SAFE

            # Exclude all genes with compound locations (as it only messes things up)
            if type(current_feature.location) is not Bio.SeqFeature.CompoundLocation and \
            type(adjacent_feature.location) is not Bio.SeqFeature.CompoundLocation:
                
                # Make Introns SeqFeature
                start_pos = SeqFeature.ExactPosition(current_feature.location.end+1)
                end_pos = SeqFeature.ExactPosition(adjacent_feature.location.start)
                try:
                    exact_location = SeqFeature.FeatureLocation(start_pos, end_pos)
                except Exception:
                    continue

                # Make Introns SeqRecord
                seq_obj = exact_location.extract(rec).seq
                seq_name = feature_name + '_' + rec.name
                seq_rec = SeqRecord(seq_obj, id=seq_name, name='', description='')
                
                # IF LENGTH OF Introns > 100 BP, ATTACH SEQRECORD TO GROWING DICTIONARY 
                if len(seq_obj) > 100:
                    if feature_name in masterdict_nucl.keys() or inv_feature_name in masterdict_nucl.keys():
                        if feature_name in masterdict_nucl.keys():
                            tmp = masterdict_nucl[feature_name]
                            tmp.append(seq_rec)
                            masterdict_nucl[feature_name] = tmp
                        if inv_feature_name in masterdict_nucl.keys():
                            pass  # Don't count Introns in the IRs twice
                    else:
                        masterdict_nucl[feature_name] = [seq_rec]
'''
########################################################################


def extract_collect_features(masterdict_nucl, inFn, feattype, fileext):
    rec = SeqIO.read(inFn, fileext)
    IRexclusion = [] # For excluding sequences that occur twice in a record due to IR duplication
    for feature in rec.features:
        if feature.type == feattype:
            if 'gene' in feature.qualifiers:
                feature_name = feature.qualifiers['gene'][0]
            elif 'note' in feature.qualifiers:
                feature_name = feature.qualifiers['note'][0]
            elif 'standard_name' in feature.qualifiers:
                feature_name = feature.qualifiers['standard_name'][0]
            else:
                sys.exit("  ERROR: Can find neither `gene` nor `note` nor `standard_name` in some or all %s qualifiers of record `%s`." % (feattype, rec.name))                
            seq_name = feature_name + '_' + rec.name
            # Nucleotide sequences
            seq_obj = feature.extract(rec).seq
            seq_rec = SeqRecord(seq_obj, id=seq_name, name='', description='')
            if feature_name in masterdict_nucl.keys():
                if feature_name not in IRexclusion: # For excluding sequences that occur twice in a record due to IR duplication
                    tmp = masterdict_nucl[feature_name]
                    tmp.append(seq_rec)
                    masterdict_nucl[feature_name] = tmp
                    IRexclusion.append(feature_name) # For excluding sequences that occur twice in a record due to IR duplication
            else:
                masterdict_nucl[feature_name] = [seq_rec]
                IRexclusion.append(feature_name) # For excluding sequences that occur twice in a record due to IR duplication


def minimal_n_taxa(my_dict, files):
    supposed_taxa = [f.rstrip(".gb") for f in files]
    for feature_name,v in my_dict.iteritems():
        actual_taxa = ["_".join(e.id.split("_")[2:]) for e in v]
        if len(actual_taxa) < len(supposed_taxa):
            missing = list(set(supposed_taxa).difference(set(actual_taxa)))
            print("For feature `%s`, the following taxa are missing: %s" % (feature_name, missing))
            del my_dict[feature_name]
    return my_dict


def main(inDir, feattype, fileext):

    # MAKE OUTPUT FOLDER
    outDir = os.path.join(inDir, 'output')
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    # EXTRACT AND COLLECT FEATURES FROM RECORDS
    files = [f for f in os.listdir(inDir) if f.endswith(fileext)]
    masterdict_nucl = collections.OrderedDict()
    for f in files:
        extract_collect_features(masterdict_nucl, os.path.join(inDir, f), feattype, fileext)

    # ONLY CONSIDER SUCH FEATURES THAT EXIST IN ALL TAXA
    #pdb.set_trace()
    masterdict_nucl = minimal_n_taxa(masterdict_nucl, files)
    

    # WRITE FEATURES TO FILE AND ALIGN
    print("WRITING FEATURES TO FILE AND ALIGNING THEM")
    if masterdict_nucl.items():
        for feature_name,v in masterdict_nucl.iteritems():
            outFn_unalign_nucl = os.path.join(outDir, 'nucl_'+feature_name+'.unalign.fasta')
            outFn_aligned_nucl = os.path.join(outDir, 'nucl_'+feature_name+'.aligned.fasta')
            # WRITE UNALIGNED PROTEIN SEQUENCES
            print("Writing and aligning `%s` ..." % (feature_name))
            with open(outFn_unalign_nucl, 'w') as hndl:
                SeqIO.write(v, hndl, 'fasta')
            # ALIGN SEQUENCES
            #import subprocess
            #subprocess.call(['mafft', '--auto', outFn_unalign_nucl, '>', outFn_aligned_nucl])
            mafft_cline = MafftCommandline(input=outFn_unalign_nucl)
            stdout, stderr = mafft_cline()
            with open(outFn_aligned_nucl, 'w') as hndl:
                hndl.write(stdout)
    if not masterdict_nucl.items():
        sys.exit('  ERROR: No items in nucleotide masterdictionary.')

    # IMPORT FEATURE ALIGNMENTS AND CONVERT (FASTA->NEXUS)
    print("CONVERTING FEATURE-WISE ALIGNMENTS (FASTA->NEXUS)")
    alignm_L = []
    for feature_name in masterdict_nucl.keys():
        aligned_nucl_fasta = os.path.join(outDir, 'nucl_'+feature_name+'.aligned.fasta')
        aligned_nucl_nexus = os.path.join(outDir, 'nucl_'+feature_name+'.aligned.nexus')
        # Convert from fasta to nexus
        print("Reading and converting alignment of `%s` ..." % (feature_name))
        try:
            alignm_fasta = AlignIO.read(aligned_nucl_fasta, 'fasta', alphabet=Alphabet.generic_dna)
            hndl = StringIO()
            AlignIO.write(alignm_fasta, hndl, 'nexus')
            nexus_string = hndl.getvalue()
            
            # Stripping the gene name from the sequence name:
            if "\n'trn" in nexus_string:
                nexus_string = nexus_string.replace("\n'"+feature_name+"_", "\n'combined_")
            else:
                nexus_string = nexus_string.replace('\n'+feature_name+'_', '\ncombined_')
            alignm_nexus = Nexus.Nexus(nexus_string)
            alignm_L.append((feature_name, alignm_nexus)) # Function 'Nexus.combine' needs a tuple.
        except:
            print '  ERROR: Cannot process alignment of', feature_name

    print("CONCATENATING ALL FEATURE-WISE ALIGNMENTS")
    # Combine the NEXUS alignments (in no particular order)
    n_aligned_features = len(alignm_L)
    alignm_combined = Nexus.combine(alignm_L) # Function 'Nexus.combine' needs a tuple.
    outFn_nucl_combined_fasta = os.path.join(outDir, 'nucl_'+str(n_aligned_features)+'combined.aligned.fasta')
    outFn_nucl_combined_nexus = os.path.join(outDir, 'nucl_'+str(n_aligned_features)+'combined.aligned.nexus')
    alignm_combined.write_nexus_data(filename=open(outFn_nucl_combined_nexus, 'w'))

    # NOTE: The "charpartition"-line of the file outFn_nucl_combined_nexus does not have single quotes around marker names that contain a dash. In such cases, it is manually necessary to add these single quotes before the below command (i.e., AlignIO.convert) is executed!

    #pdb.set_trace()
    AlignIO.convert(outFn_nucl_combined_nexus, 'nexus', outFn_nucl_combined_fasta, 'fasta')

    
############
# ARGPARSE #
############

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    
    # Required
    parser.add_argument('-i',
                        '--inp',
                        help='relative path to input directory; contains GenBank files; Example: /path_to_input/',
                        default='/home/username/Desktop/',
                        required=True)
                        
    parser.add_argument('-t',
                        '--feattype',
                        help='feature type to be extracted and aligned; choices: `intron` and `IGS`',
                        choices=['intron', 'IGS'],
                        default='IGS',
                        required=True)

    # Optional
    parser.add_argument('-f',
                        '--fileext',
                        help='File extension of input files', 
                        default='gb',
                        required=False)

    parser.add_argument('--version', 
                        help='Print version information and exit',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

########
# MAIN #
########
    main(args.inp, args.feattype, args.fileext)
