#!/usr/bin/env python3
'''Extract and align coding and non-coding regions across multiple plastid genomes'''
__version__ = 'm_gruenstaeudl@fhsu.edu|Mon 09 Oct 2023 08:18:54 PM CDT'

# ------------------------------------------------------------------------------#
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
# ------------------------------------------------------------------------------#
## DEBUGGING HELP
import ipdb


# ipdb.set_trace()
# -----------------------------------------------------------------#
# CLASSES AND FUNCTIONS

class ExtractAndCollect:

    def __init__(self, main_d_nucl):
        self.main_d_nucl = main_d_nucl

    def do_CDS(self, main_d_prot, rec, len_cutoff, log):
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
                    if gene_name in self.main_d_nucl.keys():
                        tmp = self.main_d_nucl[gene_name]
                        tmp.append(seq_rec)
                        self.main_d_nucl[gene_name] = tmp
                    else:
                        self.main_d_nucl[gene_name] = [seq_rec]
                    # Protein sequences
                    seq_obj = feature.extract(rec).seq.translate(table=11)  # , cds=True)
                    seq_rec = Bio.SeqRecord.SeqRecord(seq_obj, id=seq_name, name='', description='')
                    # Length cutoff value
                    if len(seq_obj) >= len_cutoff:
                        pass
                    else:
                        continue
                    # Foo bar baz
                    if gene_name in main_d_prot.keys():
                        tmp = main_d_prot[gene_name]
                        tmp.append(seq_rec)
                        main_d_prot[gene_name] = tmp
                    else:
                        main_d_prot[gene_name] = [seq_rec]

    def do_IGS(self, rec, fname, len_cutoff, log):
        # EXTRACT ALL GENES FROM RECORD (I.E., CDS, TRNA, RRNA)
        # Resulting list contains adjacent features in order of appearance on genome
        # Note: No need to include "if feature.type=='tRNA'", because all tRNAs are also annotated as genes
        all_genes = [feature for feature in rec.features if (
                feature.type == 'gene' and
                'gene' in feature.qualifiers
        )]
        # Note: The statement "if feature.qualifier['gene'][0] not 'matK'" is necessary, as matK is located inside trnK
        ## TO DO ##
        # Isn't there a better way to handle matK?!
        all_genes_minus_matK = [feature for feature in all_genes if feature.qualifiers['gene'][0] != 'matK']
        # LOOP THROUGH GENES
        for count, idx in enumerate(range(0, len(all_genes_minus_matK) - 1), 1):
            cur_feat = all_genes_minus_matK[idx]
            cur_feat_name = cur_feat.qualifiers['gene'][0]
            adj_feat = all_genes_minus_matK[idx + 1]
            adj_feat_name = adj_feat.qualifiers['gene'][0]
            # Define names of IGS
            if 'gene' in cur_feat.qualifiers and 'gene' in adj_feat.qualifiers:
                cur_feat_name = cur_feat_name
                cur_feat_name_SAFE = cur_feat_name.replace('-', '_')
                cur_feat_name_SAFE = re.sub(r'[^\w]', '', cur_feat_name_SAFE)
                adj_feat_name = adj_feat_name
                adj_feat_name_SAFE = adj_feat_name.replace('-', '_')
                adj_feat_name_SAFE = re.sub(r'[^\w]', '', adj_feat_name_SAFE)
                IGS_name = cur_feat_name_SAFE + '_' + adj_feat_name_SAFE
                inv_IGS_name = adj_feat_name_SAFE + '_' + cur_feat_name_SAFE
                # Exclude all genes with compound locations (as it only messes things up)
                if type(cur_feat.location) is not Bio.SeqFeature.CompoundLocation and \
                        type(adj_feat.location) is not Bio.SeqFeature.CompoundLocation:
                    # Make IGS SeqFeature
                    start_pos = Bio.SeqFeature.ExactPosition(
                        cur_feat.location.end)  # +1)  Note: It's unclear if a +1 is needed here.
                    end_pos = Bio.SeqFeature.ExactPosition(adj_feat.location.start)
                    if int(start_pos) >= int(end_pos):
                        continue  # If there is no IGS, then simply skip this gene pair
                    else:
                        try:
                            exact_location = Bio.SeqFeature.FeatureLocation(start_pos, end_pos)
                        except:
                            log.warning(
                                '\t%s: Exception occurred for IGS between `%s` (start pos: %s) and `%s` (end pos:%s ). Skipping this IGS ...' % (
                                fname, cur_feat_name, start_pos, adj_feat_name, end_pos))
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
                    if IGS_name in self.main_d_nucl.keys() or inv_IGS_name in self.main_d_nucl.keys():
                        if IGS_name in self.main_d_nucl.keys():
                            tmp = self.main_d_nucl[IGS_name]
                            tmp.append(seq_rec)
                            self.main_d_nucl[IGS_name] = tmp
                        if inv_IGS_name in self.main_d_nucl.keys():
                            pass  # Don't count IGS in the IRs twice
                    else:
                        self.main_d_nucl[IGS_name] = [seq_rec]
                # Handle genes with compound locations
                else:
                    log.warning(
                        "%s: the IGS between `%s` and `%s` is currently not handled and has to be extracted manually. Skipping this IGS ..." % (
                        fname, cur_feat_name, adj_feat_name))
                    continue

    def do_INT(self, main_d_intron2, rec, log):
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'tRNA':
                try:
                    gene_name_base = feature.qualifiers['gene'][0]
                    gene_name_base_SAFE = gene_name_base.replace('-', '_')
                    gene_name_base_SAFE = re.sub(r'[^\w]', '', gene_name_base_SAFE)
                except:
                    log.warning("Unable to extract gene name for CDS starting at `%s` of `%s`. Skipping feature ..." % (
                    feature.location.start, rec.id))
                    continue
                # Limiting the search to CDS containing introns
                if len(feature.location.parts) == 2:
                    try:
                        gene_name = gene_name_base_SAFE + "_intron1"
                        seq_rec, gene_name = extract_INT_internal(rec, feature, gene_name, 0, log)
                        if gene_name not in self.main_d_nucl.keys():
                            self.main_d_nucl[gene_name] = [seq_rec]
                        else:
                            self.main_d_nucl[gene_name].append(seq_rec)
                    except:
                        log.warning("An error for `%s` occurred" % (gene_name))
                        pass
                if len(feature.location.parts) == 3:
                    copy_feature = copy.deepcopy(
                        feature)  ## Important b/c feature is overwritten in extract_INT_internal()
                    try:
                        gene_name = gene_name_base_SAFE + "_intron1"
                        seq_rec, gene_name = extract_INT_internal(rec, feature, gene_name, 0, log)
                        if gene_name not in self.main_d_nucl.keys():
                            self.main_d_nucl[gene_name] = [seq_rec]
                        else:
                            self.main_d_nucl[gene_name].append(seq_rec)
                    except:
                        log.critical("An error occurred for `%s`" % (gene_name))
                        raise Exception()
                        # pass
                    feature = copy_feature
                    try:
                        gene_name = gene_name_base_SAFE + "_intron2"
                        seq_rec, gene_name = extract_INT_internal(rec, feature, gene_name, 1, log)
                        if gene_name not in main_d_intron2.keys():
                            main_d_intron2[gene_name] = [seq_rec]
                        else:
                            main_d_intron2[gene_name].append(seq_rec)
                    except:
                        log.warning("An issue occurred for `%s`" % (gene_name))
                        pass


# -----------------------------------------------------------------#

def proteinalign_and_backtranslate(main_d_prot, out_dir, log):
    # Step 1. Write unaligned protein sequences to file
    for k, v in main_d_prot.items():
        out_fn_unalign_prot = os.path.join(out_dir, f'prot_{k}.unalign.fasta')
        out_fn_aligned_prot = os.path.join(out_dir, f'prot_{k}.aligned.fasta')

        with open(out_fn_unalign_prot, 'w') as hndl:
            SeqIO.write(v, hndl, 'fasta')

        # Step 2. Align sequences with MAFFT
        # import subprocess
        # subprocess.call(['mafft', '--auto', out_fn_unalign_prot, '>', out_fn_aligned_prot])
        mafft_align(out_fn_unalign_prot, out_fn_aligned_prot)

    # Step 4. Check if back-translation script exists
    path_to_back_transl_helper = os.path.join(os.path.dirname(__file__), 'align_back_trans.py')
    if not os.path.isfile(path_to_back_transl_helper):
        log.critical("Unable to find `align_back_trans.py` alongside this script")
        raise Exception()

    # Step 5. Conduct actual back-translation
    for k, v in main_d_prot.items():
        out_fn_unalign_nucl = os.path.join(out_dir, f'nucl_{k}.unalign.fasta')
        out_fn_aligned_nucl = os.path.join(out_dir, f'nucl_{k}.aligned.fasta')
        out_fn_aligned_prot = os.path.join(out_dir, f'prot_{k}.aligned.fasta')

        ## TO DO - BUG! ##
        # For some reason, the path_to_back_transl_helper spits out only nexus files
        cmd = ['python3', path_to_back_transl_helper, 'fasta', out_fn_aligned_prot, out_fn_unalign_nucl,
               out_fn_aligned_nucl, '11']

        try:
            ## TO DO - BUG! ##
            # For some reason, stderr is not properly saved to log_msg in cases when the back-translation python package does not work; then log_msg is empty

            log_msg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            cmd_prt = ' '.join(cmd)
            ## TO DO - BUG! ##
            # The following lines do not work because log_msg is empty, as indicated above.

            # log_prt = str(log_msg, 'utf-8')
            # log.warning("Unable to conduct back-translation of `%s`. Command used: %s. Error message: %s" % (k, cmd_prt, log_prt))

            log.warning(
                f"Unable to conduct back-translation of `{k}`. Command used: {cmd_prt}. Error message: {e.output.decode('utf-8')}")


def mafft_align(input_file, output_file):
    # Perform sequence alignment using MAFFT
    mafft_cline = Applications.MafftCommandline(input=input_file)
    stdout, stderr = mafft_cline()

    with open(output_file, 'w') as hndl:
        hndl.write(stdout)


# ------------------------------------------------------------------------------#

def extract_INT_internal(rec, feature, gene_name, offset, log):
    try:
        feature.location = Bio.SeqFeature.FeatureLocation(feature.location.parts[offset].end,
                                                          feature.location.parts[offset + 1].start)
    except:
        feature.location = Bio.SeqFeature.FeatureLocation(feature.location.parts[offset + 1].start,
                                                          feature.location.parts[offset].end)
    try:
        seq_name = gene_name + '_' + rec.name
        seq_obj = feature.extract(rec).seq  # Here the actual extraction is conducted
        seq_rec = Bio.SeqRecord.SeqRecord(seq_obj, id=seq_name, name='', description='')
        return (seq_rec, gene_name)
    except:
        log.critical("Unable to conduct intron extraction for %s" % feature.qualifiers['gene'])
        raise Exception()


# ------------------------------------------------------------------------------#
def remove_duplicates(my_dict):
    for k, v in my_dict.items():
        idtags = []
        for counter, seqrec in enumerate(v):
            if seqrec.id in idtags:
                # v.pop(counter)
                del v[counter]
            else:
                idtags.append(seqrec.id)
        my_dict[k] = v
    # return my_dict


# ------------------------------------------------------------------------------#
## MAIN
# ------------------------------------------------------------------------------#
def main(args):
    # --------------------------------------------------------------------------#
    # UNPACKING INPUT PARAMETERS
    in_dir = args.inpd
    if not os.path.exists(in_dir):
        logging.critical("Input directory `%s` does not exist." % in_dir)
        raise Exception()
    out_dir = args.outd
    if not os.path.exists(out_dir):
        logging.critical("Output directory `%s` does not exist." % out_dir)
        raise Exception()
    fileext = args.fileext
    exclude_list = args.excllist
    len_cutoff = args.lencutoff
    tax_cutoff = args.taxcutoff
    select_mode = args.selectmode.lower()
    verbose = args.verbose
    # --------------------------------------------------------------------------#
    # SET UP LOGGER
    log = logging.getLogger(__name__)
    if verbose:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.DEBUG, logger=log)
    else:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO, logger=log)
    # --------------------------------------------------------------------------#
    # SET UP EMPTY ORDERED DICTIONARIES
    main_d_nucl = collections.OrderedDict()
    if select_mode == 'cds':
        main_d_prot = collections.OrderedDict()
    if select_mode == 'int':
        main_d_intron2 = collections.OrderedDict()
    # --------------------------------------------------------------------------#
    # EXTRACT AND COLLECT ANNOTATIONS FROM RECORDS
    action = "extracting annotations from genome records"
    log.info("%s" % action)
    ###
    files = [f for f in os.listdir(in_dir) if f.endswith(fileext)]
    for f in files:
        log.info('reading file `%s`' % (f))
        rec = Bio.SeqIO.read(os.path.join(in_dir, f), 'genbank')

        ## TO DO ##
        # Integrate taxcutoff in ExtractAndCollect(main_d_nucl).do_CDS
        #						 ExtractAndCollect(main_d_nucl).do_IGS
        #						 ExtractAndCollect(main_d_nucl).do_INT

        ## TO DO ##   
        # If warning 'BiopythonWarning: Partial codon, len(sequence) not a multiple of three.' occurs in line above:
        # Is there a way to suppress the warning in line above but activate a flag which would allow us to solve it in individual function

        if select_mode == 'cds':
            ExtractAndCollect(main_d_nucl).do_CDS(main_d_prot, rec, len_cutoff, log)
        if select_mode == 'igs':
            ExtractAndCollect(main_d_nucl).do_IGS(rec, f, len_cutoff, log)
        if select_mode == 'int':
            ExtractAndCollect(main_d_nucl).do_INT(main_d_intron2, rec, log)
            main_d_nucl.update(main_d_intron2)
        # Check if main dictionary empty
        if not main_d_nucl.items():
            log.critical("No items in main dictionary" % out_dir)
            raise Exception()
    # --------------------------------------------------------------------------#
    # REMOVE DUPLICATE ANNOTATIONS
    # Note: Not sure why I have to run this removal twice, but not all
    #       duplicates are removed first time around.
    action = "removing duplicate annotations"
    log.info("%s" % action)
    ###
    remove_duplicates(main_d_nucl)
    remove_duplicates(main_d_nucl)
    if select_mode == 'cds':
        remove_duplicates(main_d_prot)
        remove_duplicates(main_d_prot)
    # --------------------------------------------------------------------------#
    # REMOVE ANNOTATIONS THAT OCCUR IN FEWER THAN X TAXA
    action = ("removing annotations that occur in fewer than %s taxa" % tax_cutoff)
    log.info("%s" % action)
    ###
    for k, v in main_d_nucl.items():
        if len(v) < tax_cutoff:
            del main_d_nucl[k]
            if select_mode == 'cds':
                del main_d_prot[k]
    # --------------------------------------------------------------------------#
    # REMOVE ORFs
    action = "removing ORFs"
    log.info("%s" % action)
    ###
    list_of_orfs = [orf for orf in main_d_nucl.keys() if "orf" in orf]
    for orf in list_of_orfs:
        del main_d_nucl[orf]
        if select_mode == 'cds':
            del main_d_prot[orf]
    # --------------------------------------------------------------------------#
    # REMOVE USER-DEFINED GENES
    action = "removing user-defined genes"
    log.info("%s" % action)
    ###
    if exclude_list:
        if select_mode == 'int':
            to_be_excluded = [i + "_intron1" for i in exclude_list] + [i + "_intron2" for i in exclude_list]
            exclude_list = to_be_excluded
        for excluded in exclude_list:
            if excluded in main_d_nucl:
                del main_d_nucl[excluded]
                if select_mode == 'cds':
                    del main_d_prot[excluded]
            else:
                log.warning("Region `%s` to be excluded but unable to be found in infile." % excluded)
                pass
    # --------------------------------------------------------------------------#
    # EXTRACTING AND COMBINING NUCLEOTIDE SEQUENCES GENEWISE
    action = "extracting and combining nucleotide sequences genewise"
    log.info("%s" % action)
    ###
    if main_d_nucl.items():
        for k, v in main_d_nucl.items():
            out_fn_unalign_nucl = os.path.join(out_dir, 'nucl_' + k + '.unalign.fasta')
            with open(out_fn_unalign_nucl, 'w') as hndl:
                Bio.SeqIO.write(v, hndl, 'fasta')

    # --------------------------------------------------------------------------#
    # MULTIPLE SEQUENCE ALIGNMENT BASED ON NUCLEOTIDE SEQUENCE DATA
    if not select_mode == 'cds':
        action = "conducting multiple sequence alignment based on nucleotide sequence data"
        log.info("%s" % action)
        ###
        if main_d_nucl.items():
            for k, v in main_d_nucl.items():
                out_fn_unalign_nucl = os.path.join(out_dir, 'nucl_' + k + '.unalign.fasta')
                out_fn_aligned_nucl = os.path.join(out_dir, 'nucl_' + k + '.aligned.fasta')
                # import subprocess
                # subprocess.call(['mafft', '--auto', out_fn_unalign_nucl, '>', out_fn_aligned_nucl])

                ## TO DO ##
                # Automatically determine number of threads available #
                # Have the number of threads saved as num_threads
                num_threads = 1

                ## TO DO ##
                # Let user choose if alignment conducted with MAFFT, MUSCLE, CLUSTAL, etc.; use a new argparse argument and if statements in the ine below

                alignm_cmdlexec = Bio.Align.Applications.MafftCommandline(input=out_fn_unalign_nucl,
                                                                          adjustdirection=True, thread=num_threads)
                stdout, stderr = alignm_cmdlexec()
                with open(out_fn_aligned_nucl, 'w') as hndl:
                    hndl.write(stdout)
        else:
            log.critical("No items in nucleotide main dictionary to process")
            raise Exception()
    # --------------------------------------------------------------------------#
    # CONDUCT PROTEIN ALIGNMENT AND BACK-TRANSLATION TO NUCLEOTIDES
    if select_mode == 'cds':
        action = "conducting multiple sequence alignment based on protein sequence data, followed by a back-translation to nucleotides"
        log.info("%s" % action)
        ###
        try:
            proteinalign_and_backtranslate(main_d_prot, out_dir, log)
        except:
            raise Exception()
    # --------------------------------------------------------------------------#
    # CONVERT FASTA ALIGNMENT TO NEXUS ALIGNMENT AND APPEND TO LIST
    action = "converting FASTA alignment to NEXUS alignment"
    log.info("%s" % action)
    ###
    successList = []
    for k in main_d_nucl.keys():
        aligned_nucl_fasta = os.path.join(out_dir, 'nucl_' + k + '.aligned.fasta')
        aligned_nucl_nexus = os.path.join(out_dir, 'nucl_' + k + '.aligned.nexus')
        # Convert FASTA alignment to NEXUS alignment
        try:
            Bio.AlignIO.convert(aligned_nucl_fasta, 'fasta', aligned_nucl_nexus, 'nexus', molecule_type='DNA')
        except:
            log.warning("Unable to convert alignment of `%s` from FASTA to NEXUS" % k)
            # raise Exception()
            continue  # skip to next k in loop, so that k is not included in successList
        # Import NEXUS and append to list for concatenation
        try:
            # ipdb.set_trace()
            alignm_nexus = Bio.AlignIO.read(aligned_nucl_nexus, 'nexus')
            hndl = io.StringIO()
            Bio.AlignIO.write(alignm_nexus, hndl, 'nexus')
            nexus_string = hndl.getvalue()
            nexus_string = nexus_string.replace('\n' + k + '_',
                                                '\nconcatenated_')  # IMPORTANT: Stripping the gene name from the sequence name
            alignm_nexus = Bio.Nexus.Nexus.Nexus(nexus_string)
            successList.append((k, alignm_nexus))  # Function 'Bio.Nexus.Nexus.combine' needs a tuple.
        except:
            log.warning("Unable to add alignment of `%s` to concatenation" % k)
            # raise Exception()
            pass
    # --------------------------------------------------------------------------#
    # CONCATENATE ALIGNMENTS (IN NO PARTICULAR ORDER)
    action = "concatenate alignments (in no particular order)"
    log.info("%s" % action)
    ###
    # Define output names
    out_fn_nucl_concatenated_fasta = os.path.join(out_dir, 'nucl_' + str(len(successList)) + 'concatenated.aligned.fasta')
    out_fn_nucl_concatenated_nexus = os.path.join(out_dir, 'nucl_' + str(len(successList)) + 'concatenated.aligned.nexus')
    # Do concatenation
    try:
        alignm_concatenated = Bio.Nexus.Nexus.combine(successList)  # Function 'Bio.Nexus.Nexus.combine' needs a tuple.
    except:
        log.critical("Unable to concatenate alignments")
        raise Exception()
    # Write concatenated alignment to file in NEXUS format
    alignm_concatenated.write_nexus_data(filename=open(out_fn_nucl_concatenated_nexus, 'w'))
    # Convert concatenated alignment in NEXUS format to FASTA format
    Bio.AlignIO.convert(out_fn_nucl_concatenated_nexus, 'nexus', out_fn_nucl_concatenated_fasta, 'fasta')
    # --------------------------------------------------------------------------#
    # CLOSING
    action = "end of script\n"
    log.info("%s" % action)
    quit()


# ------------------------------------------------------------------------------#
# ARGPARSE
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Author|Version: ' + __version__)
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
    parser.add_argument("--verbose", "-v", action="version", version="%(prog)s " + __version__,
                        help="(Optional) Enable verbose logging", default=True)
    args = parser.parse_args()
    main(args)
# ------------------------------------------------------------------------------#
# EOF
# ------------------------------------------------------------------------------#
