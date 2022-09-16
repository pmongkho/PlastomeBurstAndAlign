def main(args):

    # ONLY CONSIDER SUCH IGS THAT EXIST IN ALL TAXA


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
