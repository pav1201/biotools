# import libraries
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os.path
from itertools import product
from Bio import SeqIO
from pybedtools import BedTool
import sys

# Parse command line arguments
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
                        description="Script to generate bed file with given sequences coordinates\nthat intersect with your target regions")
parser.add_argument("sequence", type=str, help="Your sequence")
parser.add_argument("bed", type=str, help="Path to yor bed file with regions of interest")
parser.add_argument("fasta", type=str, help="Path to reference fasta file")
parser.add_argument("-t", "--input_type", type=str, help="Type of input file", default="bed")


args = vars(parser.parse_args())


def generate_combinations(input_string):
    # define iupac nucleotide code
    iupac_nucleotides = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'U',
    'R': ['A', 'G'],
    'Y': ['C', 'T'] ,
    'S': ['G' , 'C'],
    'W': ['A' , 'T'],
    'K': ['G' , 'T'],
    'M': ['A' , 'C'],
    'B': ['C' , 'G' , 'T'],
    'D': ['A', 'G' , 'T'],
    'H': ['A' , 'C' , 'T'],
    'V': ['A', 'C' , 'G' ],
    'N': ['A', 'C', 'G', 'T']
    }
    # generate all possible sequences
    options = [iupac_nucleotides.get(char, ['Not IUPAC nucleotide']) for char in input_string]
    # check if all nucleotides are coded with UIPAC
    if ['Not IUPAC nucleotide'] in options:
        print("Input contains wrong characters")
        return
    # if everithing OK, return sqequences
    return [''.join(combination) for combination in product(*options)]

# Get all possible strings
resulting_strings = generate_combinations(args['sequence'].upper())
print("List of input sequences")
print(resulting_strings)

# define function to generate bed with query sequence
# that intersects with your bed regions
def generate_bed(ref_fasta, sequence):
    with open(f"{sequence}.bed", "w") as bed_file:
        for record in SeqIO.parse(ref_fasta, "fasta"):
            seq = str(record.seq)
            start = 0
            while True:
                start = seq.find(sequence, start)
                if start == -1:
                    break
                end = start + len(sequence)
                bed_file.write(f"{record.id}\t{start}\t{end}\n")
                start += 1  # Move to the next position
    return

# generate bed file for each sequence
for seq in resulting_strings:
    print(f"Sequence to generate bed: {seq}")
    generate_bed(args['fasta'], seq)

# intersect sequence beds
if len(resulting_strings) > 1:
    seq_bed = BedTool(f"{resulting_strings[0]}.bed")
    print(seq_bed)
    for bed in resulting_strings[1:]:
        bed_b = BedTool(f"{bed}.bed")
        print(bed_b)
        seq_bed = seq_bed.cat(bed_b, postmerge=False)
        os.remove(f"{bed}.bed")
elif len(resulting_strings) == 1:
    seq_bed = BedTool(f"{resulting_strings[0]}.bed")
else:
    print("Something wrong")

print("Intersecting reference and sequence beds")
bed_regs_of_interest = BedTool(args['bed'])
result = seq_bed.intersect(bed_regs_of_interest)

file_prefix = os.path.basename(args['bed']).replace(f".{args['input_type']}", '')
result.saveas(f"{file_prefix}_{args['sequence'].upper()}.bed")

print("Your bed is ready")
print(f"Head of final bed:")
print(f"{result.head()}")

os.remove(f"{resulting_strings[0]}.bed")
