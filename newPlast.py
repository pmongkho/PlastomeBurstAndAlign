def process_protein_alignment(k, v, out_dir, path_to_back_transl_helper, num_threads):
    # Define input and output names
    out_fn_unalign_prot = os.path.join(out_dir, f"prot_{k}.unalign.fasta")
    out_fn_aligned_prot = os.path.join(out_dir, f"prot_{k}.aligned.fasta")
    out_fn_unalign_nucl = os.path.join(out_dir, f"nucl_{k}.unalign.fasta")
    out_fn_aligned_nucl = os.path.join(out_dir, f"nucl_{k}.aligned.fasta")

    # Step 1. Write unaligned protein sequences to file
    with open(out_fn_unalign_prot, "w") as hndl:
        SeqIO.write(v, hndl, "fasta")

    # Step 2. Align matrices based on their PROTEIN sequences via third-party alignment tool
    mafft_align(out_fn_unalign_prot, out_fn_aligned_prot, num_threads)

    # Step 4. Conduct actual back-translation from PROTEINS TO NUCLEOTIDES
    # Note: For some reason, the path_to_back_transl_helper spits only works if FASTA files are specified, not if NEXUS files are specified
    cmd = [
        "python3",
        path_to_back_transl_helper,
        "fasta",
        out_fn_aligned_prot,
        out_fn_unalign_nucl,
        out_fn_aligned_nucl,
        "11",
    ]
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        cmd_prt = " ".join(cmd)
        log.warning(
            f"Unable to conduct back-translation of `{k}`. Command used: {cmd_prt}. Error message received: {e.output.decode('utf-8').strip()}."
        )


def conduct_protein_alignment_and_back_translation(main_dict_prot, out_dir):
    """
    Iterates over all unaligned PROTEIN matrices, aligns them as proteins via third-party software, and back-translates each alignment to NUCLEOTIDES
    Input:  dictionary of sorted PROTEIN sequences of all regions
    Output: aligned nucleotide matrices (present as files in NEXUS format)
    """
    action = "conducting multiple sequence alignment based on protein sequence data, followed by back-translation to nucleotides"

    log.info(action)

    # Step X. Determine number of CPU core available
    ## TO DO ##
    # Automatically determine number of threads available #
    # Have the number of threads saved as num_threads
    num_threads = os.cpu_count()  # Automatically determine number of threads available

    # Step 3. Check if back-translation script exists
    path_to_back_transl_helper = os.path.join(
        os.path.dirname(__file__), "align_back_trans.py"
    )
    if not os.path.isfile(path_to_back_transl_helper):
        log.critical("Unable to find `align_back_trans.py` alongside this script")
        raise Exception("Back-translation helper script not found")

    # Use ThreadPoolExecutor to parallelize the alignment and back-translation tasks
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        future_to_protein = {
            executor.submit(
                process_protein_alignment,
                k,
                v,
                out_dir,
                path_to_back_transl_helper,
                num_threads,
            ): k
            for k, v in main_dict_prot.items()
        }

        for future in as_completed(future_to_protein):
            k = future_to_protein[future]
            try:
                result = future.result()  # If needed, you can handle results here
                log.info(result)
            except Exception as exc:
                log.error("%r generated an exception: %s" % (k, exc))
