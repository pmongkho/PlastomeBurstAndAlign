def main(args):

            # REMOVE GENES OF EXCLUDE-LIST
            if exclude_list:
                #pdb.set_trace()
                to_be_excluded = [i+"_intron1" for i in exclude_list]+[i+"_intron2" for i in exclude_list]
                for excluded in to_be_excluded:
                    if excluded in masterdict_nucl:
                        del masterdict_nucl[excluded]
                        print('  NOTE: Gene `%s` of exclude-list detected here; thus, was excluded from output.' % (excluded))

            # ALIGN AND WRITE TO FILE
            if masterdict_nucl.items():
                for k,v in masterdict_nucl.items():
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
