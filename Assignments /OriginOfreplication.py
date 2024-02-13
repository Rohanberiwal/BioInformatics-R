from itertools import product
from Bio import Entrez
from Bio import SeqIO
import matplotlib.pyplot as plt

Entrez.email = "rohanberiwal.com"
###(T/A)TTTAT(A/G)TTT(T/A) this is the genome sequnce that i used as a autonomous replication sequnce ARS 
listDNA = [
        "TTTTATATTTT",
        "TTTTATATTTA",
        "TTTTATGTTTT",
        "TTTTATGTTTA",
        "ATTTATATTTT",
        "ATTTATGTTTT",
        "ATTTATGTTTA"
    ]

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
        
    
"""    
The algorithm iterates over the given DNA sequence in chunks of length equal to the length 
of each motif in the listDNA. For each chunk, it checks if the chunks has the listDNA element
If it is, it adds the start and end positions of the chunk to the listextra. In the end it return the psotion 
of the origin of the sequence .
The size of the chunks is found using the len(listDNA[0]) AND is equal to the len of the "TTTTATATTTT"
"""
def generator(sequence):
    listDNA = [
        "TTTTATATTTT",
        "TTTTATATTTA",
        "TTTTATGTTTT",
        "TTTTATGTTTA",
        "ATTTATATTTT",
        "ATTTATGTTTT",
        "ATTTATGTTTA"
    ]
    seqlen = len(listDNA[0])
    listextra = []
    for i in range(len(sequence) - seqlen + 1):
        nums = sequence[i:i + seqlen]
        if nums in listDNA:
            listextra.append((i, i + seqlen))
    return listextra


def printer_seq(accession, ori_sites):
    print("\n")
    print(f"Oigin of replication for chromosome {accession}:")
    counter = 1
    for start, end in ori_sites:
        print(counter,".)"f"INDEX (STARTED): {start} to  INDEX(ENDING): {end}")
        counter = counter + 1
    print("\n")

def func(total):
    print(f"Total ORI count across all chromosomes: {total}")

    
    """
    This function prints the ORI count across all chromosomes on  the console or the temrianl  .
    """
def printer(i, result):
    print(f"Genome chromosome {i}:")
    print(f"Name: {result.name}")
    print(f"ID: {result.id}")
    print(f"Description: {result.description}")
    print(f"Seq('{result.seq[:50]}')")
    print("\n")
    
def count_summary(count):
    print("Algo completed ")
    print(f"Total genome sequences printed: {count}")
    

##Typically this function is for the Plotting of the Graph of the ORI SEQUNCES 
##THIS PART OF THE CODE TAKE A TIME OF 30 SECOND - 1 minute  TO GIVE THE OUTPUT GRAPH ! 
##I HAVE MADE THIS GRAPH TO GET A BROADER AND A WIDER VIEW OF THE ORI - SEQUNCES INDEX IN Each of the chromosome

def GraphComparator(accession_list):
    plt.figure(figsize=(10, 6))
    for i in accession_list:
        handle = Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
        res = SeqIO.read(handle, "fasta")
        result =  generator(res.seq)
        x_axis = [site[0] for site in result]
        y_axis = [site[1] for site in result]
        plt.scatter(x_axis, y_axis, label=f'Chromosome {i}')
        handle.close()
    plt.title('ORI Sites for Each Chromosome')
    plt.xlabel('Start index ')
    plt.ylabel('End index ')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('ori_sites.png') 
    plt.show()
    

    """
    Main part of the code where the code start is the main() funtion 
    this prcoesses the ori sequnces using the generator function  
    and the writer function to write the output to a file
    and the count_summary function to print the total number of ORI sites across all chromosomes
    and the GraphComparator function to plot the ORI sites for each chromosome .
    the main() function takes the output file name as an argument .
    """
def main(output_file):
    total = 0
    try:
        for i in accession_list:
            handle = Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
            result = SeqIO.read(handle, "fasta")
            ori_sites = generator(result.seq)
            total = total +  len(ori_sites)
            #writer(i, result, output_file)
            printer_seq(i, ori_sites) 
            handle.close()
        func(total)
        GraphComparator(accession_list)  
        
    except Exception as exceptions:
        print("Not able to fetch results")
        print("Error -> "+str(exceptions))


output_file = "Saccharomyces.txt"
main(output_file)
