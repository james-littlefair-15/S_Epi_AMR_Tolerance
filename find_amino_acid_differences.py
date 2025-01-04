from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import argparse

def load_msa(fasta_file):
    """Load a multiple sequence alignment from a FASTA file."""
    alignment = list(SeqIO.parse(fasta_file, "fasta"))
    return MultipleSeqAlignment(alignment)

def set_reference(alignment, ref_header):
    """Set the reference sequence based on the given header name."""
    for record in alignment:
        header_parts = record.id.split('=')
        strain_name = header_parts[1]
        if strain_name == ref_header:
            return record
    raise ValueError(f"Reference sequence with header '{ref_header}' not found.")

def get_strain_name(header):
    """Extract strain name from the header."""
    return header.split('=')[1]

def is_multiple_of_three(seq):
    """Check if the sequence (excluding gaps) is a multiple of 3."""
    non_gap_seq = seq.replace('-', '')
    return len(non_gap_seq) % 3 == 0

def translate_codon(codon):
    """Translate a codon into its corresponding amino acid."""
    if codon == '---':
        return '-'
    return str(Seq(codon).translate())

def find_amino_acid_differences(reference, alignment, frameshifted_strains):
    """Find amino acid differences between each strain and the reference, excluding frameshifted strains."""
    differences = {}
    stop_codon_changes = {}
    ambiguous_changes = {}  # New column for ambiguous changes
    ref_seq = str(reference.seq)

    for record in alignment:
        strain_name = get_strain_name(record.id)

        seq = str(record.seq)
        
        if strain_name in frameshifted_strains:
            continue
        
        diff = []
        stop_changes = []
        ambiguous = []  # List to store ambiguous changes
        ref_index = 0
        strain_index = 0
        codon_index = 0

        while ref_index < len(ref_seq):
            ref_codon = ""
            ref_positions = []  # Keep track of reference positions used for this codon

            # Build reference codon dynamically (skipping gaps)
            while len(ref_codon) < 3 and ref_index < len(ref_seq):
                if ref_seq[ref_index] != '-':
                    ref_codon += ref_seq[ref_index]
                    ref_positions.append(ref_index)
                ref_index += 1

            # If codon is incomplete, stop processing
            if len(ref_codon) < 3:
                break

            # Use equivalent alignment positions from the strain
            strain_codon = "".join(seq[pos] for pos in ref_positions)

            # Only consider substitutions if the strain codon matches the reference codon in non-gap length
            if strain_codon.count('-') == 0 and len(strain_codon) == len(ref_codon):
                ref_aa = translate_codon(ref_codon)
                strain_aa = translate_codon(strain_codon)

                if ref_aa != strain_aa:
                    if strain_aa == '*':
                        stop_changes.append(f"{ref_aa}{codon_index + 1}*")
                    elif strain_aa == 'X':
                        ambiguous.append(f"{ref_aa}{codon_index + 1}X")
                    else:
                        diff.append(f"{ref_aa}{codon_index + 1}{strain_aa}")

            codon_index += 1

        differences[strain_name] = diff
        stop_codon_changes[strain_name] = stop_changes
        ambiguous_changes[strain_name] = ambiguous

    return differences, stop_codon_changes, ambiguous_changes


def find_insertions_deletions(reference, alignment):
    """Find insertions and deletions in nucleotide space relative to the reference, ignoring gaps in the reference index."""
    indels = {}
    ref_seq = str(reference.seq)
    
    for record in alignment:
        seq = str(record.seq)
        indel_events = []
        ref_index = 0
        strain_index = 0
        ref_pos = 1  # 1-based position in the reference, ignoring gaps
        
        while ref_index < len(ref_seq) and strain_index < len(seq):
            if ref_seq[ref_index] == '-' and seq[strain_index] != '-':
                start = ref_pos
                insertion = []
                while ref_index < len(ref_seq) and ref_seq[ref_index] == '-':
                    insertion.append(seq[strain_index])
                    ref_index += 1
                    strain_index += 1
                end = ref_pos
                indel_events.append(f"{start}_{end}ins{''.join(insertion)}")
            elif seq[strain_index] == '-' and ref_seq[ref_index] != '-':
                start = ref_pos
                deletion = []
                while strain_index < len(seq) and seq[strain_index] == '-':
                    deletion.append(ref_seq[ref_index])
                    ref_index += 1
                    strain_index += 1
                end = ref_pos + len(deletion) - 1
                indel_events.append(f"{start}_{end}del{''.join(deletion)}")
            else:
                if ref_seq[ref_index] != '-':
                    ref_pos += 1
                ref_index += 1
                strain_index += 1
        
        strain_name = get_strain_name(record.id)
        indels[strain_name] = indel_events
    
    return indels

def main(fasta_file, ref_header, output_file):
    alignment = load_msa(fasta_file)
    reference = set_reference(alignment, ref_header)
    
    # Check if the reference strain is frameshifted
    if not is_multiple_of_three(str(reference.seq)):
        print("Error: reference strain appears to be frameshifted. Please select a non-frameshifted reference for amino acid substitutions to be accurately determined")
        exit(1)

    frameshifted_strains = []
    
    # Check for frameshifts in other strains
    for record in alignment:
        strain_name = get_strain_name(record.id)
        if not is_multiple_of_three(str(record.seq)):
            frameshifted_strains.append(strain_name)

    differences, stop_codon_changes, ambiguous_changes = find_amino_acid_differences(reference, alignment, frameshifted_strains)
    indels = find_insertions_deletions(reference, alignment)
    
    # Collect all strain names and sort them
    strain_names = sorted([get_strain_name(record.id) for record in alignment])
    
    with open(output_file, 'w') as f:
        f.write("Strain\tSubstitutions (aa)\tAmbiguous Substitutions (aa)\tPremature Stops (aa)\tIndels (nuc)\tFrameshifted\n")
        for strain_name in strain_names:
            diff_str = ','.join(differences.get(strain_name, []))
            ambiguous_str = ','.join(ambiguous_changes.get(strain_name, []))
            stop_changes_str = ','.join(stop_codon_changes.get(strain_name, []))
            indel_str = ','.join(indels.get(strain_name, []))
            frameshifted = "Yes" if strain_name in frameshifted_strains else ""
            f.write(f"{strain_name}\t{diff_str}\t{ambiguous_str}\t{stop_changes_str}\t{indel_str}\t{frameshifted}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find amino acid differences and insertions/deletions between sequences in a multiple sequence alignment and a reference sequence")
    parser.add_argument("fasta_file", help="Input multiple sequence alignment FASTA file. Please note that it is expected that the header will contain the strain name appended with an equals sign e.g. >AnythingYouWant=MyStrainName")
    parser.add_argument("ref_header", help="Name of the reference Strain appended to respective header e.g. MyStrainName")
    parser.add_argument("output_file", help="Output TSV file to store the results")
    
    args = parser.parse_args()
    
    main(args.fasta_file, args.ref_header, args.output_file)

