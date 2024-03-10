from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

Entrez.email = "rohanberiwal.email@example.com" 
query = '"Saccharomyces cerevisiae"[Organism] AND pyruvate decarboxylase[PDC]'
accession_list = []

def blast_search(sequence):
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        blast_records = NCBIXML.parse(result_handle)
        return blast_records
    except Exception as e:
        print("Error during BLAST search:")
        print(e)
        return None

def func(filename):
    try:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=100)  
        record = Entrez.read(handle)
        handle.close()

        if record["Count"] == "0":
            print("No sequences found for the query.")
            return None
        
        liststore = record["IdList"]
        for accession_number in liststore:
            accession_list.append(accession_number)
            handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
            pdc_seq = SeqIO.read(handle, "fasta")
            handle.close()
            result_handle = NCBIWWW.qblast("blastn", "nt", pdc_seq.seq)
            blast_records = NCBIXML.parse(result_handle)
            with open("output.txt", "a") as output_file:
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            output_file.write("****Alignment****\n")
                            output_file.write(f"Accession: {alignment.accession}\n")
                            output_file.write(f"Alignment length: {hsp.align_length}\n")
                            output_file.write(f"Identity: {hsp.identities}/{hsp.align_length} ({hsp.identities / hsp.align_length:.2%})\n")
                            output_file.write(f"Query start: {hsp.query_start}\n")
                            output_file.write(f"Query end: {hsp.query_end}\n")
                            output_file.write(f"Subject start: {hsp.sbjct_start}\n")
                            output_file.write(f"Subject end: {hsp.sbjct_end}\n")
                            output_file.write(f"E-value: {hsp.expect}\n")
                            output_file.write(f"Score: {hsp.score}\n")
                            output_file.write(f"Query: {hsp.query}\n")
                            output_file.write(f"Subject: {hsp.sbjct}\n\n")

        return len(accession_list)

    except Exception as e:
        print("Error:")
        print(e)
        return None

def main():
    filename = "Sequences.txt"
    accession_count = func(filename)
    if accession_count:
        print("Code compiled successfully. Please check the file 'output.txt' for alignment blast scores.")

main()




