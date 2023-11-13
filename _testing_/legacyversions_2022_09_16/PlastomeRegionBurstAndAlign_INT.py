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

def remove_duplicates(my_dict):
    for k,v in my_dict.items():
        idtags = []
        for counter, seqrec in enumerate(v):
            if seqrec.id in idtags:
                #v.pop(counter)
                del v[counter]
            else:
                idtags.append(seqrec.id)
        my_dict[k] = v
    #return my_dict

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

    # SET UP EMPTY ORDERED DICTIONARIES
    masterdict_intron1 = collections.OrderedDict()
    masterdict_intron2 = collections.OrderedDict()

    # EXTRACT AND COLLECT INT FROM RECORDS
    files = [f for f in os.listdir(inDir) if f.endswith(fileext)]
    for f in files:
        rec = Bio.SeqIO.read(os.path.join(inDir, f), 'genbank')
        for feature in rec.features:
            if feature.type == 'CDS':
                try:
                    gene_name_base = feature.qualifiers['gene'][0]
                except:
                    #pdb.set_trace()
                    print("WARNING: Unable to extract gene name for CDS starting at `%s` of `%s`. Skipping feature ..." % (feature.location.start, rec.id))
                    continue

                ## Let's limit the search to CDS containing introns
                if len(feature.location.parts)==2:
                    try:
                        gene_name = gene_name_base + "_intron1"
                        seq_rec,gene_name = extract_INT(rec, feature, gene_name, 0)
                        if gene_name not in masterdict_intron1.keys():
                            masterdict_intron1[gene_name] = [seq_rec]
                        else:
                            masterdict_intron1[gene_name].append(seq_rec)
                    except:
                        print("An error for `%s` occurred." % (gene_name))
                        pass

                if len(feature.location.parts)==3:
                    copy_feature = copy.deepcopy(feature)  ## Important b/c feature is overwritten in extract_INT()
                    try:
                        gene_name = gene_name_base + "_intron1"
                        seq_rec,gene_name = extract_INT(rec, feature, gene_name, 0)
                        if gene_name not in masterdict_intron1.keys():
                            masterdict_intron1[gene_name] = [seq_rec]
                        else:
                            masterdict_intron1[gene_name].append(seq_rec)
                    except:
                        raise Exception("  ERROR: An error occurred for `%s`" % (gene_name))
                        #pass
                    feature = copy_feature
                    try:
                        gene_name = gene_name_base + "_intron2"
                        seq_rec,gene_name = extract_INT(rec, feature, gene_name, 1)
                        if gene_name not in masterdict_intron2.keys():
                            masterdict_intron2[gene_name] = [seq_rec]
                        else:
                            masterdict_intron2[gene_name].append(seq_rec)
                    except:
                            print("An error for `%s` occurred." % (gene_name))
                            pass

    for introndata in [masterdict_intron1, masterdict_intron2]:
        if len(introndata) > 0:

            # REMOVE ALL DUPLICATE ENTRIES
            # Note: Not sure why I have to run this removal twice, but not all
            #       duplicates are removed first time around.
            remove_duplicates(introndata)
            remove_duplicates(introndata)

            # REMOVE ORFs
            list_of_orfs = [orf for orf in introndata.keys() if "orf" in orf]
            for orf in list_of_orfs:
                del introndata[orf]

            # REMOVE GENES OF EXCLUDE-LIST
            if exclude_list:
                #pdb.set_trace()
                to_be_excluded = [i+"_intron1" for i in exclude_list]+[i+"_intron2" for i in exclude_list]
                for excluded in to_be_excluded:
                    if excluded in introndata:
                        del introndata[excluded]
                        print('  NOTE: Gene `%s` of exclude-list detected here; thus, was excluded from output.' % (excluded))

            # ALIGN AND WRITE TO FILE
            if introndata.items():
                for k,v in introndata.items():
                    outFn_unalign_nucl = os.path.join(outDir, ''+k+'.unalign.fasta')
                    outFn_aligned_nucl = os.path.join(outDir, ''+k+'.aligned.fasta')
                    with open(outFn_unalign_nucl, 'w') as hndl:
                        Bio.SeqIO.write(v, hndl, 'fasta')
                    # ALIGN SEQUENCES
                    #import subprocess
                    #subprocess.call(['mafft', '--auto', outFn_unalign_nucl, '>', outFn_aligned_nucl])
                    mafft_cline = Bio.Align.Applications.MafftCommandline(input=outFn_unalign_nucl, adjustdirection=True)
                    stdout, stderr = mafft_cline()
                    with open(outFn_aligned_nucl, 'w') as hndl:
                        hndl.write(stdout)
            if not introndata.items():
                raise Exception("  ERROR: No items in nucleotide masterdictionary for %s." % (introndata))

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
                        default=['rps12'],
                        help="(Optional) List of regions to be excluded")
    parser.add_argument("--verbose", "-v", action="version", version="%(prog)s "+__version__,
                        help="(Optional) Enable verbose logging", default=True)
    args = parser.parse_args()
    main(args)

#-----------------------------------------------------------------#
#EOF
#-----------------------------------------------------------------#
