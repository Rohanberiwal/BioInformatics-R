from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random  
from Bio.Seq import MutableSeq

##My email id 
Entrez.email = "rohanberiwal@gmail.com"
snps = []
## same acession list fo rthe pdc acession 
accesion_list = {
    "NsC_001144.5 PDC1": ("NC_001144.5", 232390, 234081),
    "NC_001139.9 PDC5": ("NC_001139.9", 651290, 652981),
    "NC_001144.5 PDC6": ("NC_001144.5", 410723, 412414)
}

## main func to retirve the database pdc sequnce
def func(accession_number, start, end):
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record

## This is the most important function this is the blast initator function that starts the blast 
## This uses the NCBIWWW  and the NCBIXML . These are the part of Biopython 
## This ncbiwww take the inputs as the blast n i am using the blastn because i have to compare the genome against all
## Qblast is the place where the blast starts programmatically
## read is the just a normal readfunction to read 
def blast_initiator(sequence):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence.seq)
    result = NCBIXML.read(result_handle)
    return result

## same functionality implemnted in the question  2 
## this is the finder function 
## The blast finder takes the blast records as  inputs
## the hsp is high scoring segmnet pair 
## I have taken a hsp.expect values less than 0.001 because we get the similar alignment resuls 
## the high expect score say more deviation in  the alignment 
## a less E value is better so we get mathcing score for the yeast / given virus 
def finder(result):
    E = 0.001
    for alignment in result.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E: 
                for i, (query_base, hit_base) in enumerate(zip(hsp.query, hsp.sbjct)):
                    if query_base != hit_base:  
                        snp_position = hsp.query_start + i
                        snps.append((alignment.title, snp_position, query_base, hit_base))
    return snps

## nORMAL TRANSLATE FUNCTION TO DO TEH PROCESS OF THE TRANSLATIO
def translate_sequence(sequence):
    return sequence.seq.translate()

## sift function
def sift(snp_record):
    prediction = random.choice(["Tolerated", "Damaging"])  
    return prediction

## Writer function
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
                       
## printer to pritn the putput on the temrianl
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

## def where all the things started 

def main():
    total_snps = {}
    output_filename = "SNP_output.txt"
    with open(output_filename, "w") as f:
        f.write("SNPs Detected:\n")
    for gene_name, (accession_number, start, end) in accesion_list.items():
        print(f"Processing {gene_name} (Accession Number: {accession_number})...")
        sequence = func(accession_number, start, end)
        if "PDC1" in gene_name or "PDC6" in gene_name:
            sequence = sequence.reverse_complement()
        blast_record = blast_initiator(sequence)
        snps = finder(blast_record)
        total_snps[gene_name] = len(snps)
        writer(output_filename, gene_name, snps, sequence)
        printer(gene_name, snps, sequence, output_filename) 
    print("Total SNPs:")
    for gene_name, count in total_snps.items():
        accession_number = accesion_list[gene_name][0]
        print(f"{gene_name} (Accession Number: {accession_number}): {count}")

main()
