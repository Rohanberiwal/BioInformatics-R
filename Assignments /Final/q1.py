from Bio import Entrez
from Bio import SeqIO
Entrez.email = "rohanberiwal@gmail.com"
query = '"Saccharomyces cerevisiae"[Organism] AND pyruvate decarboxylase[PDC]'
output_filename = "output.txt"
complement =["PDC1", "PDC6"]
### The part just assigns a list named acession list and has the pdc gene acession number 
### Why does the left side of the index 1 less than the mentioend in the articles ? This is because it is wrtten a zero based
### indexing that is done in the code 

accession_list = {
        "PDC1": ("NC_001144.5", 232389, 234081),
        "PDC5": ("NC_001144.5", 410722, 412414),
        "PDC6": ("NC_001139.9", 651289, 652981)
    }

## function to clean the file befroe writng the output to it 
def cleaner(filename):
    try:
        with open(filename, "w") as f:
            f.truncate(0)
    except Exception as e:
        print("Error")
        print(e)
        
### This is the  function that bascailly fethces the pdc sequence from the NCBI website 
### Output it text file named as output.txt

def func(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return record
    except Exception as e:
        print(e)


## a printer function to print the result on the termianl alogn wiht the index and the acession number 
def printer(accession_number, gene_name, sequence):
    print("\n")
    output = f"Accession Number: {accession_number}\nGene Name: {gene_name}\nSequence: {sequence}\n"
    print(output)
    print("\n")


## A writer function to write  the result on the Output.txt file
def writer(output, filename):
    try:
        with open(filename, "a") as f:
            f.write(output)
    except Exception as e:
        print("Error")
        print(e)

## on the website it was writtne to take the complenet of the index of pdc 
### specially for that the computation function is made , it does the revser complment of the PDC 1 AND PDC 6 
## PDC 1 is namely nc_001144.5 (232389, 234081) and pdc 6 is  NC_001139.9( 651289, 652981) 
## Note that pdc 5 is left as it is 
## the index are written alongside so that we only stick to the pdc part of the gene 
   

def computation(sequence, start_pos, end_pos, gene_name):
    try:
        if gene_name in complement:
            indx = sequence[start_pos:end_pos].reverse_complement()
        else:
            indx = sequence[start_pos:end_pos]
        return indx
    except Exception as e:
        print(e)

### the basic  main function fromwhere the call to all teh function starts 

def main() :
    cleaner(output_filename)
    for name , (accession_number, start, end) in accession_list.items():
        sequence = func(accession_number)
        if sequence:
            indexed = computation(sequence, start, end, name)
            printer(accession_number, name, indexed.seq)
            output = f"Accession Number: {accession_number}\nGene Name: {name}\nSequence: {indexed.seq}\n\n"
            writer(output, output_filename)
  
### main function call           
main()
