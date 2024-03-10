from Bio import Entrez
from Bio import SeqIO
Entrez.email = "rohanberiwal.email@example.com" 
query = '"Saccharomyces cerevisiae"[Organism] AND pyruvate decarboxylase[PDC]'
accession_list = []

def cleaner(filename):
    with open(filename, "w") as f:
        f.truncate(0)

def printer(accession_number, sequence):
    print("\n")
    print(f"Accession Number: {str(accession_number)}")
    print(f"DNA Sequence ")
    print(sequence)
    print("\n")

def writer(accession_number, sequence, filename):
    with open(filename, "a") as f:
        f.write(f"Accession Number: {accession_number}\n")
        f.write("DNA Sequence :  \n")
        f.write(str(sequence) + "\n\n")

def func(filename):
    try:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=100)  
        record = Entrez.read(handle)
        handle.close()

        if record["Count"] == "0":
            liststore = []
        else:
            liststore = record["IdList"]
            
        for accession_number in liststore:
            accession_list.append(accession_number)
            handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
            pdc_seq = SeqIO.read(handle, "fasta")
            handle.close()
            printer(accession_number, pdc_seq.seq)
            writer(accession_number, pdc_seq.seq, filename)

        return len(accession_list)

    except Exception as e:
        print("Error:")
        print(e)
        return None

def main() :
    filename="Sequences.txt"
    cleaner(filename)
    accession_count = func(filename)
    print("Total count:", accession_count)

main()
