# flake8: noqa
from rich import print
import random


DNA_CODON_TABLE = {
        # T                    # C                    # A                    # G
    "TTT": ["Phe", "F"],   "TCT": ["Ser", "S"],   "TAT": ["Tyr", "Y"],   "TGT": ["Cys", "C"],  # T
    "TTC": ["Phe", "F"],   "TCC": ["Ser", "S"],   "TAC": ["Tyr", "Y"],   "TGC": ["Cys", "C"],
    "TTA": ["Leu", "L"],   "TCA": ["Ser", "S"],   "TAA": ["STP", "*"],   "TGA": ["STP", "*"],
    "TTG": ["Leu", "L"],   "TCG": ["Ser", "S"],   "TAG": ["STP", "*"],   "TGG": ["Trp", "W"],
    # --------------------------------------------------------------------------------------------
    "CTT": ["Leu", "L"],   "CCT": ["Pro", "P"],   "CAT": ["His", "H"],   "CGT": ["Arg", "R"],  # C
    "CTC": ["Leu", "L"],   "CCC": ["Pro", "P"],   "CAC": ["His", "H"],   "CGC": ["Arg", "R"],
    "CTA": ["Leu", "L"],   "CCA": ["Pro", "P"],   "CAA": ["Gln", "Q"],   "CGA": ["Arg", "R"],
    "CTG": ["Leu", "L"],   "CCG": ["Pro", "P"],   "CAG": ["Gln", "Q"],   "CGG": ["Arg", "R"],
    # --------------------------------------------------------------------------------------------
    "ATT": ["Ile", "I"],   "ACT": ["Thr", "T"],   "AAT": ["Asn", "N"],   "AGT": ["Ser", "S"],  # A
    "ATC": ["Ile", "I"],   "ACC": ["Thr", "T"],   "AAC": ["Asn", "N"],   "AGC": ["Ser", "S"],
    "ATA": ["Ile", "I"],   "ACA": ["Thr", "T"],   "AAA": ["Lys", "K"],   "AGA": ["Arg", "R"],
    "ATG": ["Met", "M"],   "ACG": ["Thr", "T"],   "AAG": ["Lys", "K"],   "AGG": ["Arg", "R"],
    # --------------------------------------------------------------------------------------------
    "GTT": ["Val", "V"],   "GCT": ["Ala", "A"],   "GAT": ["Asp", "D"],   "GGT": ["Gly", "G"],  # G
    "GTC": ["Val", "V"],   "GCC": ["Ala", "A"],   "GAC": ["Asp", "D"],   "GGC": ["Gly", "G"],
    "GTA": ["Val", "V"],   "GCA": ["Ala", "A"],   "GAA": ["Glu", "E"],   "GGA": ["Gly", "G"],
    "GTG": ["Val", "V"],   "GCG": ["Ala", "A"],   "GAG": ["Glu", "E"],   "GGG": ["Gly", "G"],
}

AMINO_ACID_NAMES = {
    "F": ["phenylalanine", "phenylalanyl"],
    "L": ["leucine", "leucyl"],
    "I": ["isoleucine", "isoleucyl"],
    "M": ["methionine", "methionyl"],
    "V": ["valine", "valyl"],
    "S": ["serine", "seryl"],
    "P": ["proline", "prolyl"],
    "T": ["threonine", "threonyl"],
    "A": ["alanine", "alanyl"],
    "Y": ["tyrosine", "tyrosyl"],
    "H": ["histidine", "histidyl"],
    "Q": ["glutamine", "glutamyl"],
    "N": ["aspargine", "asparaginyl"],
    "K": ["lysine", "lysyl"],
    "D": ["aspartic acid", "aspartyl"],
    "E": ["glutamic acid", "glutamyl"],
    "C": ["cysteine", "cysteinyl"],
    "W": ["tytrophan", "tryptophyl"],
    "R": ["arginine", "arginyl"],
    "G": ["glycine", "glycyl"]
}

DNA_NUCLEOTIDE_PAIRS = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}


def print_format_sequence_str(s):
    """
    Color code the sequence by nucleotide
    """
    print(str := (s \
        .replace("A", "[yellow]A[/]") \
        .replace("T", "[green]T[/]") \
        .replace("U", "[salmon1]U[/]") \
        .replace("G", "[blue]G[/]") \
        .replace("C", "[magenta]C[/]")))


def print_format_amino_chain_str(a):
    """
    Color code the chain by amino acid
    """
    print(str := (a \
        .replace("F", "[cyan]F[/]") \
        .replace("L", "[yellow]L[/]") \
        .replace("I", "[gold1]I[/]") \
        .replace("M", "[green]M[/]") \
        .replace("V", "[green1]V[/]") \
        .replace("S", "[aquamarine1]S[/]") \
        .replace("P", "[magenta]P[/]") \
        .replace("T", "[light_coral]T[/]") \
        .replace("A", "[sky_blue2]A[/]") \
        .replace("Y", "[purple]Y[/]") \
        .replace("H", "[aquamarine3]H[/]") \
        .replace("Q", "[hot_pink3]Q[/]") \
        .replace("N", "[salmon1]N[/]") \
        .replace("K", "[orchid]K[/]") \
        .replace("D", "[red]D[/]") \
        .replace("E", "[dark_orange]E[/]") \
        .replace("C", "[spring_green2]C[/]") \
        .replace("W", "[sea_green1]W[/]") \
        .replace("R", "[plum4]R[/]") \
        .replace("G", "[chartreuse4]G[/]")))


def dna_to_rna(s) -> str:
    """
    Converts a DNA nucleotide sequence to an RNA nucleotide sequence.
    """
    return s \
        .replace("T", "U")


def rna_to_dna(s) -> str:
    """
    Converts a RNA nucleotide sequence to an DNA nucleotide sequence.
    """
    return s \
        .replace("U", "T")


def codons(s) -> str:
    """
    Returns a nucleotide sequence separated by codon.
    """
    s = list(s)
    i = 3
    while i < len(s):
        s.insert(i, " ")
        i += 4
    return "".join(s)


def dna_to_amino(s, letterOutput=True) -> str:
    """
    Converts a DNA nucleotide sequence to an amino acid chain.

    Essentially, protein biosynthesis without a visit to the Golgi apparatus.
    """
    global DNA_CODON_TABLE

    amino_acids = ""

    for codon in codons(s).split(" "):
        amino_acid = DNA_CODON_TABLE.get(codon)[letterOutput]

        if amino_acid == "*" or amino_acid == "STP":
            break
        amino_acids += amino_acid

        if letterOutput is False:
            amino_acids += " "
    return amino_acids


def rna_to_amino(s, letterOutput=True) -> str:
    """
    Converts an RNA nucleotide sequence to an amino acid chain.

    The subprocess of translation in protein biosynthesis.
    """
    sRna = rna_to_dna(s)
    return dna_to_amino(sRna, letterOutput=letterOutput)


def dna_inverse_seq(s) -> str:
    """
    Returns antinucleotide sequence for a DNA nucleotide sequence.

    `A->T; G->C; T->A; C->G`
    """
    global DNA_NUCLEOTIDE_PAIRS
    inverse = ""
    for nucleotide in s:
        inverse += DNA_NUCLEOTIDE_PAIRS.get(nucleotide)
    return inverse


def rna_inverse_seq(s) -> str:
    """
    Returns antinucleotide sequence for an RNA nucleotide sequence.

    `A->U; G->C; U->A; C->G`
    """
    sRna = rna_to_dna(s)
    return dna_to_rna(dna_inverse_seq(sRna))


def chemical_name_amino_chain(a) -> str:
    """
    Returns the chemical name of an amino acid sequence with letters.
    """
    global AMINO_ACID_NAMES
    name = ""
    i = 0
    for amino in a:
        if i != len(a) - 1:
            name += AMINO_ACID_NAMES[amino][1]
        else:
            name += AMINO_ACID_NAMES[amino][0]
        i += 1
    return name


def amino_to_dna(a) -> str:
    """
    Converts a letter-based amino acid chain to a DNA nucleotide sequence.

    Results may vary each time the function is run, owing to the fact that
    multiple codons may correspond to a single amino acid.
    """
    global DNA_CODON_TABLE
    s = ""
    for amino in a:
        codon_set = []
        for codon, amino_letter in DNA_CODON_TABLE.items():
            if amino_letter[1] == amino:
                codon_set.append(codon)
        s += codon_set[random.randint(0, len(codon_set) - 1)]
    return s


def amino_to_rna(a) -> str:
    """
    Converts a letter-based amino acid chain to an RNA nucleotide sequence.

    Results may vary each time the function is run, owing to the fact that
    multiple codons may correspond to a single amino acid.
    """
    return dna_to_rna(amino_to_dna(a))


def chemical_name_to_dna(c) -> str:
    """
    Converts the chemical name of an amino acid chain into a
    sequence of DNA nucleotides.
    """
    global AMINO_ACID_NAMES
    c = c.replace("yl", "yl-")
    amino_letter_string = ""
    amino_chem_names = c.split("-")
    for amino_chem_name in amino_chem_names:
        for letter_code, amino_full_name in AMINO_ACID_NAMES.items():
            if amino_full_name[0] == amino_chem_name or amino_full_name[1] == amino_chem_name:
                amino_letter_string += letter_code
    return amino_to_dna(amino_letter_string)


def chemical_name_to_rna(c) -> str:
    """
    Converts the chemical name of an amino acid chain into a
    sequence of RNA nucleotides.
    """
    return dna_to_rna(chemical_name_to_dna(c))
