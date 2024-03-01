from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def compute_distance_matrix(sequences):
    # Compute a distance matrix based on pairwise sequence comparisons
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(sequences)
    return distance_matrix

def build_guide_tree(distance_matrix):
    # Build a guide tree using the distance matrix (e.g., using neighbor-joining algorithm)
    constructor = DistanceTreeConstructor()
    guide_tree = constructor.upgma(distance_matrix)
    return guide_tree

def progressive_alignment(sequences, guide_tree):
    # Perform progressive alignment using the guide tree
    aligned_sequences = []
    for record in SeqIO.parse(sequences, 'fasta'):
        aligned_sequences.append(record.seq)
    align = MultipleSeqAlignment(aligned_sequences)
    return align

def refine_alignment(alignment):
    # Optionally, refine the alignment through iterative refinement steps
    muscle_cline = MuscleCommandline(input="input.fasta", out="output.fasta")
    muscle_cline()
    alignment = AlignIO.read("output.fasta", "fasta")
    return alignment

def muscle_alignment(input_file):
    # Perform MUSCLE alignment
    distance_matrix = compute_distance_matrix(input_file)
    guide_tree = build_guide_tree(distance_matrix)
    alignment = progressive_alignment(input_file, guide_tree)
    alignment = refine_alignment(alignment)
    return alignment

# Example usage:
input_file = "sequences.fasta"
alignment = muscle_alignment(input_file)
print(alignment)
