from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
## my email id 

Entrez.email = "rohanberiwal@gmail.com"
## acession ist wiht all the pdc gene acession numebr and the pdc gene indexed where i can get the sequences 
snps = []
accession_list = {
    "NC_001144.5 PDC1": ("NC_001144.5", 232389, 234081),
    "NC_001139.9 PDC6": ("NC_001139.9", 651289, 652981),
    "NC_001144.5 PDC5": ("NC_001144.5", 410722, 412414)
}
## A dictionary for the appedning of the Snps counts for teh pdc 1 , 5,6 
total_snps = {}
## I am writing the SNp detals for the pdc into a file named SNP_output 
## here i AM appeding teh acession numebr , snp genome that has the pdc and poition of the snp 
output_filename = "SNP_output.txt"

## function to fetch the pdc gene from the NCBI database 
def func(accession_number, start, end):
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
    record = handle.read()
    handle.close()
    return record

## modified cleaner function to clean teh file 
## it gebnerates the output otherwise 
def cleaner(filename):
    try:
        with open(filename, "w") as f:
            f.truncate(0)
    except Exception as e:
        print("Error")
        print(e)

## a printer function to print on  the termianl window / termianl
def printer(accession_number, gene_name, sequence):
    print("\n")
    output = f"Accession Number: {accession_number}\nGene Name: {gene_name}\nSequence: {sequence}\n"
    print(output)
    print("\n")

## This is the most important function this is the blast initator function that starts the blast 
## This uses the NCBIWWW  and the NCBIXML . These are the part of Biopython 
## This ncbiwww take the inputs as the blast n i am using the blastn because i have to compare the genome against all
## Qblast is the place where the blast starts programmatically
## read is the just a normal readfunction to read 
def blast_initiator(sequence):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    result = NCBIXML.read(result_handle)
    return result

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

## Just a vasic function to just reverse complment 
def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())

## a normal writer fucntion to write the fucntion output in the file  / termianl 
def writer(filename, gene_name, snps, sequence):
    with open(filename, "a") as f:
        f.write(f"Gene Name: {gene_name}\n")
        f.write(f"Normal Sequence: {sequence}\n")
        f.write("SNPs:\n")
        for alignment_title, position, query_base, hit_base in snps:
            f.write(f"  Alignment Title: {alignment_title}, Position: {position}, Query Base: {query_base}, Hit Base: {hit_base}\n")
            f.write(f"  SNP-containing Sequence: {sequence[:position-1]}{hit_base}{sequence[position:]}\n")


## main function where everything  starts 
def main():
    with open(output_filename, "w") as f:
        f.write("SNPs Detected:\n")
    for gene_name, (accession_number, start, end) in accession_list.items():
        ## debugging statements 
        print(f"Processing {gene_name} (Accession Number: {accession_number})...")
        sequence = func(accession_number, start, end)
        if "PDC1" in gene_name or "PDC6" in gene_name:
            sequence = reverse_complement(sequence)
        blast_record = blast_initiator(sequence)
        snps = finder(blast_record)
        total_snps[gene_name] = len(snps)
        writer(output_filename, gene_name, snps, sequence)
    print("Total SNPs:")
    for gene_name, count in total_snps.items():
        accession_number = accession_list[gene_name][0]
        print(f"{gene_name} (Accession Number: {accession_number}): {count}")

## normal routine main call 
main()
