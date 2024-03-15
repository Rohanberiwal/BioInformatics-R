from Bio import Entrez

Entrez.email = "rohanberiwal@gmail.com"
query = '"Saccharomyces cerevisiae"[Organism] AND pyruvate decarboxylase[PDC]'
accession_list = {
    "NC_001144.5": (232390, 234081),
    "NC_001139.9": (651290, 652981),
    "NC_001144.5 PDC Gene 1": (410723, 412414)
}

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

def func(accession_number, start, end):
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
    record = handle.read()
    handle.close()
    return record

def main():
    output_filename = "output.txt"
    cleaner(output_filename)
    for accession_number, index_range in accession_list.items():
        start, end = index_range
        sequence = func(accession_number, start, end)
        writer(accession_number, sequence, output_filename)
        printer(accession_number, sequence)

main()

