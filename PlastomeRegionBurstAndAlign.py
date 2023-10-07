#!/usr/bin/env python3
'''Extract and align coding and non-coding regions across multiple plastid genomes'''
__version__ = 'm.gruenstaeudl@fu-berlin.de|2022-09-19T11:07:41 CEST'
#------------------------------------------------------------------------------#
## IMPORTS
import argparse
import Bio
from Bio import SeqIO  # line is necessary for similar reason stated in: https://www.biostars.org/p/13099/
from Bio.Align import Applications  # line is necessary for similar reason stated in: https://www.biostars.org/p/13099/
import coloredlogs
import collections
import copy
import io
import logging
import os
import re
import subprocess
import sys
#------------------------------------------------------------------------------#
## DEBUGGING HELP
import ipdb
#ipdb.set_trace()
#-----------------------------------------------------------------#
# CLASSES AND FUNCTIONS

class ExtractAndCollect:

    def __init__(self, mainD_nucl):
        self.mainD_nucl = mainD_nucl

    def do_CDS(self, mainD_prot, rec, len_cutoff, log):
        for feature in rec.features:
            if feature.type == 'CDS':
                if 'gene' in feature.qualifiers:
                    gene_name = feature.qualifiers['gene'][0]
                    seq_name = gene_name + '_' + rec.name
                    # Nucleotide sequences
                    
                    ## TO DO ##
                    # If warning 'BiopythonWarning: Partial codon, len(sequence) not a multiple of three.' occurs in line below:
                    # Add function to check if len(sequence) not a multiple of three, as error message for NC_013553.gb suggests
                    
                    seq_obj = feature.extract(rec).seq
                    seq_rec = Bio.SeqRecord.SeqRecord(seq_obj, id=seq_name, name='', description='')
                    if gene_name in self.mainD_nucl.keys():
                        tmp = self.mainD_nucl[gene_name]
                        tmp.append(seq_rec)
                        self.mainD_nucl[gene_name] = tmp
                    else:
                        self.mainD_nucl[gene_name] = [seq_rec]
                    # Protein sequences
                    seq_obj = feature.extract(rec).seq.translate(table=11)#, cds=True)
                    seq_rec = Bio.SeqRecord.SeqRecord(seq_obj, id=seq_name, name='', description='')
                    # Length cutoff value
                    if len(seq_obj) >= len_cutoff:
                        pass
                    else:
                        continue
                    # Foo bar baz
                    if gene_name in mainD_prot.keys():
                        tmp = mainD_prot[gene_name]
                        tmp.append(seq_rec)
                        mainD_prot[gene_name] = tmp
                    else:
                        mainD_prot[gene_name] = [seq_rec]

    def do_IGS(self, rec, fname, len_cutoff, log):
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
                            log.warning('\t%s: Exception occurred for IGS between `%s` (start pos: %s) and `%s` (end pos:%s ). Skipping this IGS ...' % (fname, cur_feat_name, start_pos, adj_feat_name, end_pos))
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
                    if IGS_name in self.mainD_nucl.keys() or inv_IGS_name in self.mainD_nucl.keys():
                        if IGS_name in self.mainD_nucl.keys():
                            tmp = self.mainD_nucl[IGS_name]
                            tmp.append(seq_rec)
                            self.mainD_nucl[IGS_name] = tmp
                        if inv_IGS_name in self.mainD_nucl.keys():
                            pass  # Don't count IGS in the IRs twice
                    else:
                        self.mainD_nucl[IGS_name] = [seq_rec]
                # Handle genes with compound locations
                else:
                    log.warning("\t%s: the IGS between `%s` and `%s` is currently not handled and has to be extracted manually. Skipping this IGS ..." % (fname, cur_feat_name, adj_feat_name))
                    continue

    def do_INT(self, mainD_intron2, rec, log):
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'tRNA':
                try:
                    gene_name_base = feature.qualifiers['gene'][0]
                except:
                    log.warning("\tUnable to extract gene name for CDS starting at `%s` of `%s`. Skipping feature ..." % (feature.location.start, rec.id))
                    continue
                # Limiting the search to CDS containing introns
                if len(feature.location.parts)==2:
                    try:
                        gene_name = gene_name_base + "_intron1"
                        seq_rec,gene_name = extract_INT_internal(rec, feature, gene_name, 0, log)
                        if gene_name not in self.mainD_nucl.keys():
                            self.mainD_nucl[gene_name] = [seq_rec]
                        else:
                            self.mainD_nucl[gene_name].append(seq_rec)
                    except:
                        log.warning("\tAn error for `%s` occurred" % (gene_name))
                        pass
                if len(feature.location.parts)==3:
                    copy_feature = copy.deepcopy(feature)  ## Important b/c feature is overwritten in extract_INT_internal()
                    try:
                        gene_name = gene_name_base + "_intron1"
                        seq_rec,gene_name = extract_INT_internal(rec, feature, gene_name, 0, log)
                        if gene_name not in self.mainD_nucl.keys():
                            self.mainD_nucl[gene_name] = [seq_rec]
                        else:
                            self.mainD_nucl[gene_name].append(seq_rec)
                    except:
                        log.critical("\tAn error occurred for `%s`" % (gene_name))
                        raise Exception()
                        #pass
                    feature = copy_feature
                    try:
                        gene_name = gene_name_base + "_intron2"
                        seq_rec,gene_name = extract_INT_internal(rec, feature, gene_name, 1, log)
                        if gene_name not in mainD_intron2.keys():
                            mainD_intron2[gene_name] = [seq_rec]
                        else:
                            mainD_intron2[gene_name].append(seq_rec)
                    except:
                        log.warning("\tAn issue occurred for `%s`" % (gene_name))
                        pass

#-----------------------------------------------------------------#

def conduct_backtranslation(mainD_prot, outDir, log):
        # Step 1. Write unaligned protein sequences to file
        if mainD_prot.items():
            for k,v in mainD_prot.items():
                outFn_unalign_prot = os.path.join(outDir, 'prot_'+k+'.unalign.fasta')
                outFn_aligned_prot = os.path.join(outDir, 'prot_'+k+'.aligned.fasta')
                with open(outFn_unalign_prot, 'w') as hndl:
                    Bio.SeqIO.write(v, hndl, 'fasta')
                # Step 2. Align sequences
                #import subprocess
                #subprocess.call(['mafft', '--auto', outFn_unalign_prot, '>', outFn_aligned_prot])
                mafft_cline = Bio.Align.Applications.MafftCommandline(input=outFn_unalign_prot)
                stdout, stderr = mafft_cline()
                # Step 3. Write aligned protein sequences to file
                with open(outFn_aligned_prot, 'w') as hndl:
                    hndl.write(stdout)
        else:
            log.critical("\tNo items in protein main dictionary")
            raise Exception()

        # Step 4. Check if back-translation script existant
        path_to_this_script = os.path.dirname(os.path.realpath(__file__))
        path_to_back_transl_helper = path_to_this_script + '/align_back_trans.py'
        if not os.path.isfile(path_to_back_transl_helper):
            log.critical("\tCannot find `align_back_trans.py` not alongside this script")
            raise Exception()

        # Step 5. Conduct actual back-translation via Python script by Peter Cook
        # https://github.com/peterjc/pico_galaxy/tree/master/tools/align_back_trans
        for k,v in mainD_prot.items():
            outFn_unalign_nucl = os.path.join(outDir, 'nucl_'+k+'.unalign.fasta')
            outFn_aligned_nucl = os.path.join(outDir, 'nucl_'+k+'.aligned.fasta')
            outFn_aligned_prot = os.path.join(outDir, 'prot_'+k+'.aligned.fasta')
            cmd = ['python3', path_to_back_transl_helper, 'fasta', outFn_aligned_prot, outFn_unalign_nucl, outFn_aligned_nucl, '11']
            try:
                
                ## TO DO ##
                # For some reason, stderr is not properly saved to log_msg in cases when the back-translation python package does not work; then log_msg is empty
                
                log_msg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except:
                cmd_prt = ' '.join(cmd)
                log_prt = str(log_msg, 'utf-8')
                log.critical("\tCannot conduct back-translation of `%s`. Command used: %s. Error message: %s" % (k, cmd_prt, log_prt))
                #raise Exception()

#------------------------------------------------------------------------------#

def extract_INT_internal(rec, feature, gene_name, offset, log):
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
                log.critical("\tCannot conduct intron extraction for %s" % feature.qualifiers['gene'])
                raise Exception()

#------------------------------------------------------------------------------#
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

#------------------------------------------------------------------------------#
## MAIN
#------------------------------------------------------------------------------#
def main(args):

    #--------------------------------------------------------------------------#
    # UNPACKING INPUT PARAMETERS
    inDir = args.inpd
    if not os.path.exists(inDir):
        log.critical("\tInput directory `%s` does not exist." % inDir)
        raise Exception()
    outDir = args.outd
    if not os.path.exists(outDir):
        log.critical("\tOutput directory `%s` does not exist." % outDir)
        raise Exception()
    fileext = args.fileext
    exclude_list = args.excllist
    len_cutoff = args.lencutoff
    tax_cutoff = args.taxcutoff
    select_mode = args.selectmode.lower()
    verbose = args.verbose
    #--------------------------------------------------------------------------#
    # SET UP LOGGER
    log = logging.getLogger(__name__)
    if verbose:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.DEBUG, logger=log)
    else:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO, logger=log)
    #--------------------------------------------------------------------------#
    # SET UP EMPTY ORDERED DICTIONARIES
    mainD_nucl = collections.OrderedDict()
    if select_mode == 'cds':
        mainD_prot = collections.OrderedDict()
    if select_mode == 'int':
        mainD_intron2 = collections.OrderedDict()
    #--------------------------------------------------------------------------#
    # EXTRACT AND COLLECT ANNOTATIONS FROM RECORDS
    action = "extracting annotations from genome records"
    log.info("%s" % action)
    ###
    files = [f for f in os.listdir(inDir) if f.endswith(fileext)]
    for f in files:
        log.info('reading file `%s`' % (f))
        rec = Bio.SeqIO.read(os.path.join(inDir, f), 'genbank')

        ## TO DO ##   
        # If warning 'BiopythonWarning: Partial codon, len(sequence) not a multiple of three.' occurs in line above:
        # Is there a way to suppress the warning in line above but activate a flag which would allow us to solve it in individual function
        
        if select_mode == 'cds':
            ExtractAndCollect(mainD_nucl).do_CDS(mainD_prot, rec, len_cutoff, log)
        if select_mode == 'igs':
            ExtractAndCollect(mainD_nucl).do_IGS(rec, f, len_cutoff, log)
        if select_mode == 'int':
            ExtractAndCollect(mainD_nucl).do_INT(mainD_intron2, rec, log)
            mainD_nucl.update(mainD_intron2)
        # Check if main dictionary empty
        if not mainD_nucl.items():
            log.critical("\tNo items in main dictionary" % outDir)
            raise Exception()
    #--------------------------------------------------------------------------#
    # REMOVE DUPLICATE ANNOTATIONS
        # Note: Not sure why I have to run this removal twice, but not all
        #       duplicates are removed first time around.
    action = "removing duplicate annotations"
    log.info("%s" % action)
    ###
    remove_duplicates(mainD_nucl)
    remove_duplicates(mainD_nucl)
    if select_mode == 'cds':
        remove_duplicates(mainD_prot)
        remove_duplicates(mainD_prot)
    #--------------------------------------------------------------------------#
    # REMOVE ANNOTATIONS THAT OCCUR IN FEWER THAN X TAXA
    action = ("removing annotations that occur in fewer than %s taxa" % tax_cutoff)
    log.info("%s" % action)
    ###
    for k,v in mainD_nucl.items():
        if len(v) < tax_cutoff:
            del mainD_nucl[k]
            if select_mode == 'cds':
                del mainD_prot[k]
    #--------------------------------------------------------------------------#
    # REMOVE ORFs
    action = "removing ORFs"
    log.info("%s" % action)
    ###
    list_of_orfs = [orf for orf in mainD_nucl.keys() if "orf" in orf]
    for orf in list_of_orfs:
        del mainD_nucl[orf]
        if select_mode == 'cds':
            del mainD_prot[orf]
    #--------------------------------------------------------------------------#
    # REMOVE USER-DEFINED GENES
    action = "removing user-defined genes"
    log.info("%s" % action)
    ###
    if exclude_list:
        if select_mode == 'int':
            to_be_excluded = [i+"_intron1" for i in exclude_list]+[i+"_intron2" for i in exclude_list]
            exclude_list = to_be_excluded
        for excluded in exclude_list:
            if excluded in mainD_nucl:
                del mainD_nucl[excluded]
                if select_mode == 'cds':
                    del mainD_prot[excluded]
            else:
                log.warning("\tRegion `%s` to be excluded but cannot be found in infile." % excluded)
                pass
    #--------------------------------------------------------------------------#
    # ALIGNING SEQUENCES BASED ON NUCLEOTIDES
    action = "aligning sequences based on nucleotides"
    log.info("%s" % action)
    ###
    if mainD_nucl.items():
        for k,v in mainD_nucl.items():
            outFn_unalign_nucl = os.path.join(outDir, 'nucl_'+k+'.unalign.fasta')
            outFn_aligned_nucl = os.path.join(outDir, 'nucl_'+k+'.aligned.fasta')
            with open(outFn_unalign_nucl, 'w') as hndl:
                Bio.SeqIO.write(v, hndl, 'fasta')
            #import subprocess
            #subprocess.call(['mafft', '--auto', outFn_unalign_nucl, '>', outFn_aligned_nucl])
            mafft_cline = Bio.Align.Applications.MafftCommandline(input=outFn_unalign_nucl, adjustdirection=True)
            stdout, stderr = mafft_cline()
            with open(outFn_aligned_nucl, 'w') as hndl:
                hndl.write(stdout)
    else:
        log.critical("\tNo items in nucleotide main dictionary to process")
        raise Exception()
    #--------------------------------------------------------------------------#
    # CONDUCT BACK-TRANSLATION
    if select_mode == 'cds':
        action = "conducting back-translation"
        log.info("%s" % action)
        ###
        try:
            conduct_backtranslation(mainD_prot, outDir, log)
        except:
            log.critical("\tBack-translation of that item not produced ...")
    #--------------------------------------------------------------------------#
    # CONVERT FASTA ALIGNMENT TO NEXUS ALIGNMENT AND APPEND TO LIST
    action = "converting FASTA alignment to NEXUS alignment"
    log.info("%s" % action)
    ###
    alignm_L = []
    for k in mainD_nucl.keys():
        aligned_nucl_fasta = os.path.join(outDir, 'nucl_'+k+'.aligned.fasta')
        aligned_nucl_nexus = os.path.join(outDir, 'nucl_'+k+'.aligned.nexus')
        # Convert FASTA alignment to NEXUS alignment
        try:
            Bio.AlignIO.convert(aligned_nucl_fasta, 'fasta', aligned_nucl_nexus, 'nexus', molecule_type='DNA')
        except:
            log.critical("\tCannot convert alignment of `%s` from FASTA to NEXUS" % k)
            raise Exception()
        # Import NEXUS and append to list for concatenation
        try:
            alignm_nexus = Bio.AlignIO.read(aligned_nucl_nexus, 'nexus')
            hndl = io.StringIO()
            Bio.AlignIO.write(alignm_nexus, hndl, 'nexus')
            nexus_string = hndl.getvalue()
            nexus_string = nexus_string.replace('\n'+k+'_', '\ncombined_')  # IMPORTANT: Stripping the gene name from the sequence name
            alignm_nexus = Bio.Nexus.Nexus.Nexus(nexus_string)
            alignm_L.append((k, alignm_nexus)) # Function 'Bio.Nexus.Nexus.combine' needs a tuple.
        except:
            log.critical("\tCannot add alignment of `%s` to concatenation" % k)
            raise Exception()
    #--------------------------------------------------------------------------#
    # COMBINE NEXUS ALIGNMENTS (IN NO PARTICULAR ORDER)
    action = "combining NEXUS alignments (in no particular order)"
    log.info("%s" % action)
    ###
    alignm_combined = Bio.Nexus.Nexus.combine(alignm_L) # Function 'Bio.Nexus.Nexus.combine' needs a tuple.
    outFn_nucl_combined_fasta = os.path.join(outDir, 'nucl_'+str(len(alignm_L))+'combined.aligned.fasta')
    outFn_nucl_combined_nexus = os.path.join(outDir, 'nucl_'+str(len(alignm_L))+'combined.aligned.nexus')
    alignm_combined.write_nexus_data(filename=open(outFn_nucl_combined_nexus, 'w'))
    Bio.AlignIO.convert(outFn_nucl_combined_nexus, 'nexus', outFn_nucl_combined_fasta, 'fasta')
    #--------------------------------------------------------------------------#
    # CLOSING
    action = "end of script\n"
    log.info("%s" % action)
    quit()

#------------------------------------------------------------------------------#
# ARGPARSE
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
    parser.add_argument("--selectmode", "-s", type=str, required=False,
                        help="(Optional) Type of regions to be extracted (i.e. `cds`, `int`, or `igs`)",
                        default="cds")
    parser.add_argument("--fileext", "-f", type=str, required=False,
                        help="(Optional) File extension of input files",
                        default=".gb")
    parser.add_argument("--excllist", "-e", type=list, required=False,
                        default=['rps12'],
                        help="(Optional) List of genes to be excluded")
    parser.add_argument("--lencutoff", "-l", type=int, required=False,
                        help="(Optional) Sequence length below which regions will not be extracted",
                        default=1)
    parser.add_argument("--taxcutoff", "-t", type=int, required=False,
                        help="(Optional) Minimum number of taxa in which a region must be present to be extracted",
                        default=1)
    parser.add_argument("--verbose", "-v", action="version", version="%(prog)s "+__version__,
                        help="(Optional) Enable verbose logging", default=True)
    args = parser.parse_args()
    main(args)
#------------------------------------------------------------------------------#
#EOF
#------------------------------------------------------------------------------#
