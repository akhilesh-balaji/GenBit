# GenBit

Operate on DNA/RNA genetic sequences using a computer, the same way protein biosynthesis would happen within a cell. For instance, a given DNA sequence can be converted into an RNA sequence, and then an amino acid sequence. The full name of the protein then generated can also be found. You can also find the complentary genetic sequence, of the non-coding strand involved.

## Coming Soon

- The ability to randomly mutate sequences
- Randomly generate sequences
- Finding the non-chemical name of the protein, and its functionality
- Building a system of proteins that can work together (protein simulations)

## Usage

Install the package with `pip install dnaseq-genbit`. Now, it can be imported with `import genbit`, and you now have access to all the methods in GenBit. Sequences are represented by strings, and it is these strings that are operated on.

### Example

```python
import genbit

aa_seq = 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRXXVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRXXVHFTAEEKAAVTSLWSKMNVEEAGGEALGRLLVVYPWTQRFFDSFGNLSSPSAILGNPKVKAHGKKVLTSFGDAIKNMDNLKPAFAKLSELHCDKLHVDPENFKLLGNVMVIILATHFGKEFTPEVQAAWQKLVSAVAIALAHKYXXVHFTAEEKAAVTSLWSKMNVEEAGGEALGRLLVVYPWTQRFFDSFGNLSSPSAILGNPKVKAHGKKVLTSFGDAIKNMDNLKPAFAKLSELHCDKLHVDPENFKLLGNVMVIILATHFGKEFTPEVQAAWQKLVSAVAIALAHKYXX'  # Amino Acid Sequence

rna_seq = genbit.amino_to_rna(aa_seq)  # RNA Sequence
words_amino_seq = genbit.chemical_name_amino_chain(aa_seq)  # Chemical Name of Amino Acid Sequence
dna_seq = genbit.chemical_name_to_dna(words_amino_seq)  # DNA Sequence
codon_dna_seq = genbit.codons(dna_seq)  # Codon Sequence
```
