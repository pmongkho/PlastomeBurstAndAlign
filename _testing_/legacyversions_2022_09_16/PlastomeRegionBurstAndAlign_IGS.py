#!/usr/bin/env python3
'''Extract And Align Intergenic Spacer Regions Across Multiple Plastomes'''
__version__ = 'm.gruenstaeudl@fu-berlin.de|2022-09-16T17:22:33 CEST'

#-----------------------------------------------------------------#
## IMPORTS
import argparse
import Bio
from Bio import SeqIO  # line is necessary for similar reason stated in: https://www.biostars.org/p/13099/
from Bio.Align import Applications
import collections
import copy
import io
import os
import re
import sys

#-----------------------------------------------------------------#
# DEBUGGING HELP
import pdb
#pdb.set_trace()

#-----------------------------------------------------------------#
# FUNCTIONS

# Loop through all records
# Save all item of single record to same dictionary

# For each key in dictionary
# Align records in list

def extract_collect_IGS(fname, masterdict_nucl, inFn, len_cutoff):
    rec = SeqIO.read(inFn, 'genbank')

    # EXTRACT ALL GENES FROM RECORD (I.E., CDS, TRNA, RRNA)
    # Resulting list contains adjacent features in order of appearance on genome
    # Note 1: No need to include "if feature.type=='tRNA'", because all tRNAs are also annotated as genes
    # Note 2: The statement "if feature.qualifier['gene'][0] not 'matK'" is necessary, as matK is located inside trnK
    all_genes = [feature for feature in rec.features if (feature.type=='gene' and feature.qualifiers['gene'][0]!='matK')]

    # LOOP THROUGH GENES
    for count, idx in enumerate(range(0, len(all_genes)-1), 1):
        cur_feat = all_genes[idx]
        cur_feat_name = cur_feat.qualifiers['gene'][0]
        adj_feat = all_genes[idx+1]
        adj_feat_name = adj_feat.qualifiers['gene'][0]
        #print('  %s %s Analyzing genes `%s` and `%s`' % (fname, count, cur_feat_name, adj_feat_name))

        # Define names of IGS
        if 'gene' in cur_feat.qualifiers and 'gene' in adj_feat.qualifiers:
            cur_feat_name = cur_feat_name
            cur_feat_name_SAFE = cur_feat_name.replace('-','_')
            cur_feat_name_SAFE = re.sub(r'[^\w]', '', cur_feat_name_SAFE)
            adj_feat_name = adj_feat_name
            adj_feat_name_SAFE = adj_feat_name.replace('-','_')
            adj_feat_name_SAFE = re.sub(r'[^\w]', '', adj_feat_name_SAFE)
            IGS_name = cur_feat_name_SAFE + '_' + adj_feat_name_SAFE
            inv_IGS_name = adj_feat_name_SAFE + '_' + cur_feat_name_SAFE

            # Exclude all genes with compound locations (as it only messes things up)
            if type(cur_feat.location) is not Bio.SeqFeature.CompoundLocation and \
            type(adj_feat.location) is not Bio.SeqFeature.CompoundLocation:

                # Make IGS SeqFeature
                start_pos = Bio.SeqFeature.ExactPosition(cur_feat.location.end) #+1)  Note: It's unclear if a +1 is needed here.
                end_pos = Bio.SeqFeature.ExactPosition(adj_feat.location.start)
                if int(start_pos) >= int(end_pos):
                    continue    # If there is no IGS, then simply skip this gene pair
                else:
                    try:
                        exact_location = Bio.SeqFeature.FeatureLocation(start_pos, end_pos)
                    except:
                        print('  WARNING: %s Exception occurred for IGS between `%s` (start pos: %s) and `%s` (end pos:%s ) . Skipping this IGS ...' % (fname, cur_feat_name, start_pos, adj_feat_name, end_pos))
                        continue

                # Make IGS SeqRecord
                seq_obj = exact_location.extract(rec).seq
                seq_name = IGS_name + '_' + rec.name
                seq_rec = Bio.SeqRecord.SeqRecord(seq_obj, id=seq_name, name='', description='')

                # Length cutoff value
                if len(seq_obj) >= len_cutoff:
                    pass
                else:
                    continue

                # ATTACH SEQRECORD TO GROWING DICTIONARY
                if IGS_name in masterdict_nucl.keys() or inv_IGS_name in masterdict_nucl.keys():
                    if IGS_name in masterdict_nucl.keys():
                        tmp = masterdict_nucl[IGS_name]
                        tmp.append(seq_rec)
                        masterdict_nucl[IGS_name] = tmp
                        #print('  %s `%s` processed successfully' % (fname, IGS_name))
                    if inv_IGS_name in masterdict_nucl.keys():
                        pass  # Don't count IGS in the IRs twice
                else:
                    masterdict_nucl[IGS_name] = [seq_rec]

            # Handle genes with compound locations
            else:
                print('  WARNING: %s The IGS between `%s` and `%s` is currently not handled and has to be manually extracted.' % (fname, cur_feat_name, adj_feat_name))
                continue


#def minimal_n_taxa(my_dict, minimal_n):
#    for k,v in my_dict.items():
#        idtags = []
#        if len(v) < minimal_n:
#            del my_dict[k]
#    return my_dict


def main(args):

    # UNPACKING INPUT PARAMETERS
    inDir = args.inpd
    if not os.path.exists(inDir):
        raise Exception("  ERROR: Input directory `%s` does not exist." % inDir)
    outDir = args.outd
    if not os.path.exists(outDir):
        raise Exception("  ERROR: Output directory `%s` does not exist." % outDir)
    fileext = args.fileext
    exclude_list = args.excllist
    len_cutoff = args.lencutoff

    # SET UP EMPTY ORDERED DICTIONARIES
    masterdict_nucl = collections.OrderedDict()

    # EXTRACT AND COLLECT IGS FROM RECORDS
    files = [f for f in os.listdir(inDir) if f.endswith(fileext)]
    for f in files:
        print('\nProcessing file `%s`' % (f))
        extract_collect_IGS(f, masterdict_nucl, os.path.join(inDir, f), len_cutoff)

    # REMOVE ORFs (OPTIONAL!)
    list_of_orfs = [orf for orf in masterdict_nucl.keys() if "orf" in orf]
    for orf in list_of_orfs:
        del masterdict_nucl[orf]

    # ONLY CONSIDER SUCH IGS THAT EXIST IN ALL TAXA
    #masterdict_nucl = minimal_n_taxa(masterdict_nucl, len(files))

    # REMOVE SPECIFIC GENES (OPTIONAL!)
    if exclude_list:
        for excluded in exclude_list:
            if excluded in masterdict_nucl:
                del masterdict_nucl[excluded]
            else:
                raise Exception("  ERROR: Region to be excluded (%s) not found in infile." % excluded)

    # CONVERT FASTA ALIGNMENT TO NEXUS ALIGNMENT AND APPEND FOR CONCATENATION
    alignm_L = []
    for k in masterdict_prot.keys():
        aligned_nucl_fasta = os.path.join(outDir, 'nucl_'+k+'.aligned.fasta')
        aligned_nucl_nexus = os.path.join(outDir, 'nucl_'+k+'.aligned.nexus')

        # CONVERT FASTA ALIGNMENT TO NEXUS ALIGNMENT
        try:
            Bio.AlignIO.convert(aligned_nucl_fasta, 'fasta', aligned_nucl_nexus, 'nexus', molecule_type='DNA')
        except:
            raise Exception("  ERROR: Cannot convert alignment of `%s` from FASTA to NEXUS" % k)

        # IMPORT NEXUS AND APPEND TO LIST FOR CONCATENATION
        try:
            alignm_nexus = Bio.AlignIO.read(aligned_nucl_nexus, 'nexus')
            hndl = io.StringIO()
            Bio.AlignIO.write(alignm_nexus, hndl, 'nexus')
            nexus_string = hndl.getvalue()
            nexus_string = nexus_string.replace('\n'+k+'_', '\ncombined_')  # IMPORTANT: Stripping the gene name from the sequence name
            alignm_nexus = Bio.Nexus.Nexus.Nexus(nexus_string)
            alignm_L.append((k, alignm_nexus)) # Function 'Bio.Nexus.Nexus.combine' needs a tuple.
        except:
            raise Exception("  ERROR: Cannot add alignment of `%s` to concatenation" % k)

    # COMBINE NEXUS ALIGNMENTS (IN NO PARTICULAR ORDER)
    alignm_combined = Bio.Nexus.Nexus.combine(alignm_L) # Function 'Bio.Nexus.Nexus.combine' needs a tuple.
    outFn_nucl_combined_fasta = os.path.join(outDir, 'nucl_'+str(len(alignm_L))+'combined.aligned.fasta')
    outFn_nucl_combined_nexus = os.path.join(outDir, 'nucl_'+str(len(alignm_L))+'combined.aligned.nexus')
    alignm_combined.write_nexus_data(filename=open(outFn_nucl_combined_nexus, 'w'))
    Bio.AlignIO.convert(outFn_nucl_combined_nexus, 'nexus', outFn_nucl_combined_fasta, 'fasta')

#-----------------------------------------------------------------#
# MAIN
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Author|Version: '+__version__)
    # Required
    parser.add_argument("--inpd", "-i", type=str, required=True,
                        help="path to input directory (which contains the GenBank files)",
                        default="./input")
    # Optional
    parser.add_argument("--outd", "-o", type=str, required=False,
                        help="(Optional) Path to output directory",
                        default="./output")
    parser.add_argument("--fileext", "-f", type=str, required=False,
                        help="(Optional) File extension of input files",
                        default=".gb")
    parser.add_argument("--excllist", "-e", type=list, required=False,
                        default=['rpl12_rpl2'],
                        #default=[
                        #'trnI_CAU_ycf2',
                        #'rpl23_trnI_CAU',
                        #'trnM_CAU_ycf2',
                        #'rpl23_trnM_CAU',
                        #'ndhI_ndhH',
                        #'clpP_rpl20',
                        #'rpl2_trnI_CAU',
                        #'ycf4_cemA',
                        #'ycf3_trnS_GGA',
                        #'trnV_GAC_trnI_GAU',
                        #'trnI_CAU_trnI_CAU',
                        #'trnH_GUG_trnK_UUU',
                        #'trnH_GUG_rpl2',
                        #'psaI_ycf4',
                        #'psaA_ycf3',
                        #'trnM_CAU_rps14',
                        #'trnl_CAU_ycf2',
                        #'trnI_CAU_trnL_CAA',
                        #'trnI_CAU_rpl12',
                        #'trnG_UCC_trnM_CAU',
                        #'trnG_UCC_atpA',
                        #'rrn16_rps7',
                        #'rps19_trnI_CAU',
                        #'rps19_rpl23',
                        #'rps15_ndhF',
                        #'rpl23_trnl_CAU',
                        #'rpl12_rpl2',
                        #'psbZ_trnG_UCC',
                        #'psbT_pbf1',
                        #'psbH_rpoA',
                        #'psbH_petBgene',
                        #'petBgene_petD',
                        #'pbf1_psbH',
                        #'matK_rps16']
                        help="(Optional) List of regions to be excluded")
    parser.add_argument("--lencutoff", "-l", type=int, required=False,
                        help="(Optional) Sequence length below which regions will not be extracted",
                        default=1)
    parser.add_argument("--verbose", "-v", action="version", version="%(prog)s "+__version__,
                        help="(Optional) Enable verbose logging", default=True)
    args = parser.parse_args()
    main(args)

#-----------------------------------------------------------------#
#EOF
#-----------------------------------------------------------------#
