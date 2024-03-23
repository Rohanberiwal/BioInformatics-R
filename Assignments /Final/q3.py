from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
import requests 
## My mail for the Nbci server Request 
Entrez.email = "rohanberiwal@gmail.com"
## dictionary for the saving  of the total Snp counts 
total_snps = {}
## dicationry to store the translated sequence 
## what are these trasnlted sequnce ?
## This disctionary stores the normal and the SNP contiaing gene transltion
## what is translation  ? conversion fo the mrna to protiens 
translated_sequences = {}  
## Just a file for the saving of the ouput protiens 
output_filename = "SNP_output.txt"
## A query to retrieve the sequence from the database
query = '"Saccharomyces cerevisiae"[Organism] AND pyruvate decarboxylase[PDC]'
## This is a Url that uses the EXpasy tool kit 
## I have used this website to get teh infrotiuon about  the -> https://web.expasy.org/protparam/protparam-doc.html
## I used the documation for writing the Expasy code 
## please refer to the readme for more

url = "https://web.expasy.org/cgi-bin/translate/dna2aa.cgi"
snps = []
accession_list = {
    "NC_001144.5 PDC1": ("NC_001144.5", 232389, 234081),
    "NC_001139.9 PDC6": ("NC_001139.9", 651289, 652981),
    "NC_001144.5 PDC5": ("NC_001144.5", 410722, 412414)
}


## function to fetch the pdc gene from the NCBI database 
def func(accession_number, start, end):
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
    record = handle.read()
    handle.close()
    return record

## This is the most important function this is the blast initator function that starts the blast 
## This uses the NCBIWWW  and the NCBIXML . These are the part of Biopython 
## This ncbiwww take the inputs as the blast n i am using the blastn because i have to compare the genome against all
## Qblast is the place where the blast starts programmatically
## read is the just a normal readfunction to read 

def blast_initiator(sequence):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence , entrez_query = "Saccharomyces cerevisiae[Organism]")
    blast_record = NCBIXML.read(result_handle)
    return blast_record

## same functionality implemnted in the question  2 
## this is the finder function 
## The blast finder takes the blast records as  inputs
## the hsp is high scoring segmnet pair 
## I have taken a hsp.expect values less than 0.001 because we get the similar alignment resuls 
## the high expect score say more deviation in  the alignment 
## a less E value is better so we get mathcing score for the yeast / given virus 
def finder(blast_record):
    E = 0.001
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.001: 
                for i, (query_base, hit_base) in enumerate(zip(hsp.query, hsp.sbjct)):
                    if query_base != hit_base:  
                        snp_position = hsp.query_start + i
                        snps.append((alignment.title, snp_position, query_base, hit_base))
    return snps

## revrse the pdc 1 and pdc 6 just mentioned on the ncbi 
def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())

## There is a dict_extra that send teh data to the sevrer where the EXpasy does the computaiton  
## the dict_exta has teh format defined i choose the fasta format and the sequnce for the translation 
## the translated-out is the rseult 
## finally return the output 

def translate(sequence):
        dict_extra = {
            "Original Snp": sequence,
            "format": "fasta"
        }
        results = requests.post(url, data=dict_extra)
        translated_out = results.text.split('<pre>')[-1].split('</pre>')[0].strip()
        return translated_out

## just a writer funcrtionb to write the output 
def writer(filename, gene_name, snps, sequence):
    with open(filename, "a") as f:
        f.write(f"Gene Name: {gene_name}\n")
        f.write(f"Normal Sequence: {sequence}\n")
        f.write("Translated Normal Sequence: " + translate(sequence) + "\n")
        f.write("SNPs:\n")
        for alignment_title, position, query_base, hit_base in snps:
            f.write(f"  Alignment Title: {alignment_title}, Position: {position}, Query Base: {query_base}, Hit Base: {hit_base}\n")
            modified_sequence = sequence[:position-1] + hit_base + sequence[position:]
            f.write(f"  SNP-containing Sequence: {modified_sequence}\n")
            f.write("  Translated SNP-containing Sequence: " + translate(modified_sequence) + "\n")

## main frunction where it all starts 
def main():
    with open(output_filename, "w") as f:
        f.write("SNPs Detected:\n")
    for gene_name, (accession_number, start, end) in accession_list.items():
        print(f"Processing {gene_name} (Accession Number: {accession_number})...")
        sequence = func(accession_number, start, end)
        if "PDC1" in gene_name or "PDC6" in gene_name:
            sequence = reverse_complement(sequence)
        blast_record = blast_initiator(sequence)
        snps = finder(blast_record)
        total_snps[gene_name] = len(snps)
        writer(output_filename, gene_name, snps, sequence)
        translated_sequences[gene_name] = translate(sequence)  
    print("Total SNPs:")
    for gene_name, count in total_snps.items():
        accession_number = accession_list[gene_name][0]
        print(f"{gene_name} (Accession Number: {accession_number}): {count}")
    print("\nTranslated Sequences:")
    for gene_name, sequence in translated_sequences.items():
        print(f"{gene_name}: {sequence}")

## main routine 
main()

