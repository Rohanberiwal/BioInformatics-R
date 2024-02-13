from Bio import Entrez
from Bio import SeqIO
from typing import List
Entrez.email = "rohanberiwal.com"
accession_list = [
    "NC_001133", "NC_001134", "NC_001135", "NC_001136", "NC_001137",
    "NC_001138", "NC_001139", "NC_001140", "NC_001141", "NC_001142",
    "NC_001143", "NC_001144", "NC_001145", "NC_001146", "NC_001147",
    "NC_001148", "NC_001224"
]
def cleaner(file_path):
    with open(file_path, "w") as f:
        f.truncate(0)

def writer(accession, result, output_file):
    with open(output_file, "a") as f:
        f.write(f"Name: {result.name}\n")
        f.write(f"ID: {result.id}\n")
        f.write(f"Description: {result.description}\n")
        f.write(f"Seq('{result.seq[:-1]}')")
        f.write("\n")
        f.write("\n")

def printer(accession, result):
    print(f"Genome chromosome {accession}:")
    print(f"Name: {result.name}")
    print(f"ID: {result.id}")
    print(f"Description: {result.description}")
    print(f"Seq('{result.seq[:50]}')")
    print("\n")
    
def count_summary(count):
    print("Genome for Baker yeast Done")
    print(f"Total genome sequences printed: {count}")

def main(output_file):
    count = 0
    cleaner(output_file)
    try:
        for accession in accession_list:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            result  = SeqIO.read(handle, "fasta")
            print("\n")
            print("Generating sequence for ",result.id)
            writer(accession, result , output_file)
            printer(accession, result)
            count += 1
            handle.close()
            
        count_summary(count)

    except Exception as exceptions:
        print("Not able to fetch results")
        print("Error -> "+exceptions)

output_file = "Saccharomyces.txt"
main(output_file)
