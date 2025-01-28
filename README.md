# Staphylococus epidermidis AMR Tolerance Study

This repository contains code and data from my contribution to a study on antimicrobial resistance (AMR) tolerance in *Staphylococcus epidermidis*. For full methods, refer to the original manuscript or contact the corresponding author (carolin.kobras@path.ox.ac.uk).

## Usage `find_amino_acid_differences.py` Script

This script identifies amino acid differences and insertions/deletions between sequences in a protein multiple sequence alignment with respect to a specified reference. This was used on core genes in the study to identiify variants associated with resistant phenotypes.

```bash
python3 find_amino_acid_differences.py [-h] fasta_file ref_header output_file

positional arguments:
  fasta_file   Input multiple sequence alignment protein FASTA file. Please note that it is expected
               that the header will contain the strain name appended with an equals sign
               e.g. >AnythingYouWant=MyStrainName. The multiple alignment is expected to be in the
               same orientation for each strain starting at the N-terminus (i.e. with methionine as
               the leftmost amino acid)
  ref_header   Name of the reference strain appended to respective header e.g. MyStrainName
  output_file  Output TSV file to store the results

options:
  -h, --help   show this help message and exit
```

## Data

Data (supplementary material) referred to in the methods section have been included in the 'data' folder of this repository.

These include:
  - S1: spreadsheet containing the source for the gene queries used in the study
  - S2: the protein fasta file countanining these queries
  - S3: the gene coordinates in the original strain assemblies
  - S4: the gene presence-absence matrix, with putitive frameshifts, premature stop codons, and low identity
  - S5: the variants for each core AMR-associated gene
