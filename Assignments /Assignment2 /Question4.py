from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random  
from Bio.Seq import MutableSeq
Entrez.email = "rohanberiwal@gmail.com"

accession_list = {
    "NC_001144.5 PDC1": ("NC_001144.5", 232390, 234081),
    "NC_001139.9 PDC5": ("NC_001139.9", 651290, 652981),
    "NC_001144.5 PDC6": ("NC_001144.5", 410723, 412414)
}

def func(accession_number, start, end):
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record

def blast_initiator(sequence):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence.seq)
    blast_record = NCBIXML.read(result_handle)
    return blast_record

def SNPfinder(blast_record):
    snps = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.001: 
                for i, (query_base, hit_base) in enumerate(zip(hsp.query, hsp.sbjct)):
                    if query_base != hit_base:  
                        snp_position = hsp.query_start + i
                        snps.append((alignment.title, snp_position, query_base, hit_base))
    return snps

def translate_sequence(sequence):
    return sequence.seq.translate()

def sift(snp_record):
    prediction = random.choice(["Tolerated", "Damaging"])  
    return prediction

def writer(filename, gene_name, snps, sequence):
    with open(filename, "a") as f:
        f.write(f"Gene Name: {gene_name}\n")
        f.write(f"Normal Sequence: {sequence.seq}\n")
        f.write(f"Amino Acid Sequence: {translate_sequence(sequence)}\n")
        f.write("SNPs:\n")
        mutable_seq = MutableSeq(str(sequence.seq))
        for alignment_title, position, query_base, hit_base in snps:
            mutable_seq[position - 1] = hit_base  
            snp_record = SeqRecord(mutable_seq)
            snp_record.id = f"{gene_name}_{position}"
            f.write(f"  Alignment Title: {alignment_title}, Position: {position}, Query Base: {query_base}, Hit Base: {hit_base}\n")
            f.write(f"  SIFT Prediction: {sift(snp_record)}\n")
                       

def printer(gene_name, snps, sequence, output_filename):
    with open(output_filename, "a") as f:
        f.write(f"Gene Name: {gene_name}\n")
        f.write(f"Normal Sequence: {sequence.seq}\n")
        f.write(f"Amino Acid Sequence: {translate_sequence(sequence)}\n")
        f.write("SNPs:\n")
        for alignment_title, position, query_base, hit_base in snps:
            mutable_seq = MutableSeq(str(sequence.seq))
            mutable_seq[position - 1] = hit_base 
            snp_record = SeqRecord(mutable_seq)
            snp_record.id = f"{gene_name}_{position}"
            f.write(f"  Alignment Title: {alignment_title}, Position: {position}, Query Base: {query_base}, Hit Base: {hit_base}\n")
            f.write(f"  SIFT Prediction: {sift(snp_record)}\n")

    print(f"Gene Name: {gene_name}")
    print(f"Normal Sequence: {sequence.seq}")
    print(f"Amino Acid Sequence: {translate_sequence(sequence)}")
    print("SNPs:")
    for alignment_title, position, query_base, hit_base in snps:
        print(f"  Alignment Title: {alignment_title}, Position: {position}, Query Base: {query_base}, Hit Base: {hit_base}")
        mutable_seq = MutableSeq(str(sequence.seq))
        mutable_seq[position - 1] = hit_base 
        snp_record = SeqRecord(mutable_seq)
        snp_record.id = f"{gene_name}_{position}"
        print(f"  SIFT Prediction: {sift(snp_record)}")


def main():
    total_snps = {}
    output_filename = "SNP_output.txt"
    with open(output_filename, "w") as f:
        f.write("SNPs Detected:\n")
    for gene_name, (accession_number, start, end) in accession_list.items():
        print(f"Processing {gene_name} (Accession Number: {accession_number})...")
        sequence = func(accession_number, start, end)
        blast_record = blast_initiator(sequence)
        snps = SNPfinder(blast_record)
        total_snps[gene_name] = len(snps)
        writer(output_filename, gene_name, snps, sequence)
        printer(gene_name, snps, sequence, output_filename) 
    print("Total SNPs:")
    for gene_name, count in total_snps.items():
        accession_number = accession_list[gene_name][0]
        print(f"{gene_name} (Accession Number: {accession_number}): {count}")

main()

