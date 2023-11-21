#!/usr/bin/env python
"""Back-translate a protein alignment to nucleotides
This tool is a short Python script (using Biopython library functions) to
load a protein alignment, and matching nucleotide FASTA file of unaligned
sequences, in order to produce a codon aware nucleotide alignment - which
can be viewed as a back translation.
"""

import sys
from Bio.Seq import Seq

# from Bio.Alphabet import generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Data.CodonTable import ambiguous_generic_by_id


def check_trans(identifier, nuc, prot, table):
    """Returns nucleotide sequence if works (can remove trailing stop)"""
    if len(nuc) % 3:
        sys.exit(
            f"Nucleotide sequence for {identifier} is length {len(nuc)} (not a multiple of three)"
        )

    p = str(prot).upper().replace("*", "X")
    t = str(nuc.translate(table)).upper().replace("*", "X")
    if len(t) == len(p) + 1:
        if str(nuc)[-3:].upper() in ambiguous_generic_by_id[table].stop_codons:
            # Allow this...
            t = t[:-1]
            nuc = nuc[:-3]  # edit return value
    if len(t) != len(p):
        err = (
            f"Inconsistent lengths for {identifier}, ungapped protein {len(p)}, "
            f"tripled {len(p) * 3} vs ungapped nucleotide {len(nuc)}."
        )
        if t.endswith(p):
            err += f"\nThere are {len(t) - len(p)} extra nucleotides at the start."
        elif t.startswith(p):
            err += f"\nThere are {len(t) - len(p)} extra nucleotides at the end."
        elif p in t:
            err += "\nHowever, protein sequence found within translated nucleotides."
        elif p[1:] in t:
            err += "\nHowever, ignoring first amino acid, protein sequence found within translated nucleotides."
        sys.exit(err)

    if t == p:
        return nuc
    elif p.startswith("M") and "M" + t[1:] == p:
        if str(nuc[0:3]).upper() in ambiguous_generic_by_id[table].start_codons:
            return nuc
        else:
            sys.exit(
                f"Translation check failed for {identifier}\n"
                f"Would match if {nuc[0:3].upper()} was a start codon (check correct table used)"
            )

    else:
        m = "".join("." if x == y else "!" for (x, y) in zip(p, t))
        if len(prot) < 70:
            sys.stderr.write(f"Protein:     {p}\n")
            sys.stderr.write(f"             {m}\n")
            sys.stderr.write(f"Translation: {t}\n")
        else:
            for offset in range(0, len(p), 60):
                sys.stderr.write(f"Protein:     {p[offset:offset + 60]}\n")
                sys.stderr.write(f"             {m[offset:offset + 60]}\n")
                sys.stderr.write(f"Translation: {t[offset:offset + 60]}\n\n")
        sys.exit(f"Translation check failed for {identifier}\n")


def sequence_back_translate(
    aligned_protein_record, unaligned_nucleotide_record, gap, table=0
):
    if not gap or len(gap) != 1:
        raise ValueError("Please supply a single gap character")

    ######
    # Modification on 09-Sep-2022 by Michael Gruenstaeudl
    # alpha = unaligned_nucleotide_record.seq.alphabet
    # if hasattr(alpha, "gap_char"):
    #    gap_codon = alpha.gap_char * 3
    #    assert len(gap_codon) == 3
    # else:
    #    from Bio.Alphabet import Gapped
    #    alpha = Gapped(alpha, gap)
    #    gap_codon = gap * 3
    ######

    gap_codon = "-" * 3

    ungapped_protein = aligned_protein_record.seq.ungap(gap)
    ungapped_nucleotide = unaligned_nucleotide_record.seq
    if table:
        ungapped_nucleotide = check_trans(
            aligned_protein_record.id, ungapped_nucleotide, ungapped_protein, table
        )
    elif len(ungapped_protein) * 3 != len(ungapped_nucleotide):
        sys.exit(
            f"Inconsistent lengths for {aligned_protein_record.id}, ungapped protein {len(ungapped_protein)}, "
            f"tripled {len(ungapped_protein) * 3} vs ungapped nucleotide {len(ungapped_nucleotide)}"
        )

    seq = []
    nuc = str(ungapped_nucleotide)
    for amino_acid in aligned_protein_record.seq:
        if amino_acid == gap:
            seq.append(gap_codon)
        else:
            seq.append(nuc[:3])
            nuc = nuc[3:]
    assert (
        not nuc
    ), f"Nucleotide sequence for {unaligned_nucleotide_record.id} longer than protein {aligned_protein_record.id}"

    aligned_nuc = unaligned_nucleotide_record[:]  # copy for most annotation
    aligned_nuc.letter_annotation = {}  # clear this
    aligned_nuc.seq = Seq(
        "".join(seq)
    )  # , alpha)  # Modification on 09-Sep-2022 by Michael Gruenstaeudl
    assert len(aligned_protein_record.seq) * 3 == len(aligned_nuc)
    return aligned_nuc


def alignment_back_translate(
    protein_alignment, nucleotide_records, key_function=None, gap=None, table=0
):
    """Thread nucleotide sequences onto a protein alignment."""
    if key_function is None:
        key_function = lambda x: x
    if gap is None:
        gap = "-"

    aligned = []
    for protein in protein_alignment:
        try:
            nucleotide = nucleotide_records[key_function(protein.id)]
        except KeyError:
            raise ValueError(
                f"Could not find nucleotide sequence for protein {protein.id}"
            )
        aligned.append(sequence_back_translate(protein, nucleotide, gap, table))
    return MultipleSeqAlignment(aligned)


if len(sys.argv) < 4:
    sys.exit(0)

align_format, prot_align_file, nuc_fasta_file = sys.argv[1:4]
nuc_align_file = sys.stdout
table = 0

if len(sys.argv) >= 5:
    nuc_align_file = sys.argv[4]

if len(sys.argv) == 6:
    table = int(sys.argv[5])

try:
    table = int(table)
except ValueError:
    sys.exit(f"Bad table argument {table}")

prot_align = AlignIO.read(
    prot_align_file, align_format
)  # , alphabet=generic_protein)  # Modification on 09-Sep-2022 by Michael Gruenstaeudl
nuc_dict = SeqIO.index(nuc_fasta_file, "fasta")
nuc_align = alignment_back_translate(prot_align, nuc_dict, gap="-", table=table)
AlignIO.write(nuc_align, nuc_align_file, align_format)
