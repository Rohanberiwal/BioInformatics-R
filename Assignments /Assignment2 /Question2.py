from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
## my email id 
Entrez.email = "rohanberiwal@gmail.com"
## acession ist wiht all the pdc gene acession numebr and the pdc gene indexed where i can get the sequences 
accession_list = {
    "NC_001144.5 PDC1": ("NC_001144.5", 232389, 234081),
    "NC_001139.9 PDC6": ("NC_001139.9", 651289, 652981),
    "NC_001144.5 PDC5": ("NC_001144.5", 410722, 412414)
}
liststorage  = ["PDC1" ,"PDC6"]
## I am writing the SNp detals for the pdc into a file named SNP_output 
## here i AM appeding teh acession numebr , snp genome that has the pdc and poition of the snp 
output_file = "SNP_output.txt"
snps = []
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
    for align in result.alignments:
        for hsp in align.hsps:
            if hsp.expect < E: 
                for i, (orignial_base, hit_base) in enumerate(zip(hsp.query, hsp.sbjct)):
                    if orignial_base != hit_base:  
                        snp_position = hsp.query_start + i
                        snps.append((align.title, snp_position, orignial_base, hit_base))
    return snps

## Just a vasic function to just reverse complment 
def reverse_complement(sequence_name):
    return str(Seq(sequence_name).reverse_complement())

## a normal writer fucntion to write the fucntion output in the file  / termianl 
## reasom for takign the pos - 1 since we are considerign the 1 based indexing and the last has o be avoided

def writer(filename, gname, snps, sequence_name):
    with open(filename, "a") as f:
        f.write(f"Gene Name: {gname}")
        f.write("\n")
        f.write(f"Normal Sequence: {sequence_name}")
        f.write("\n")
        f.write("SNPs:")
        f.write("\n")
        for align, pos, query_base, hits_base in snps:
            f.write(f"  Alignment Title: {align}, Position: {pos}, Query Base: {query_base}, Hit Base: {hits_base}")
            f.write("\n")
            f.write(f"  SNP-containing Sequence: {sequence_name[:pos-1]}{hits_base}{sequence_name[pos:]}")
            f.write("\n")

### just a function to print out which pdc is gettign executed 
def debugger(gname , accession_number):
    print(f"Processing {gname} (Accession Number: {accession_number})...")

## main function where everything  starts 
def main():
    ## A dictionary for the appedning of the Snps counts for teh pdc 1 , 5,6 
    total_snps = {}
    with open(output_file, "w") as f:
        f.write("SNPs Detected:\n")
        ##This below has gene name in the dic wiht the tuple as the value
    for gname, (accession_number, start, end) in accession_list.items():
        ## debugging statements 
        debugger(gname  , accession_number)
        sequence_name = func(accession_number, start, end)
        if "PDC1" in gname  or "PDC6" in gname:
            sequence_name = reverse_complement(sequence_name)
        res = blast_initiator(sequence_name)
        snps = finder(res)
        total_snps[gname] = len(snps)
        writer(output_file, gname, snps, sequence_name)

    differences = {}
    differences["pdc1"] = total_snps.get("NC_001144.5 PDC1", 0) 
    differences["pdc5"] = total_snps.get("NC_001144.5 PDC5", 0) - total_snps.get("NC_001139.9 PDC6", 0)
    differences["pdc6"] = total_snps.get("NC_001139.9 PDC6", 0) - total_snps.get("NC_001144.5 PDC1", 0)

    print("SNP counts :")
    for pair, difference in differences.items():
        print(f"{pair.upper()}: {difference}")
    print("\n")

    print("SNPs(Index in the file):")
    for gname, count in total_snps.items():
        accession_number = accession_list[gname][0]
        print(f"{gname} (Accession Number: {accession_number}): {count}")

## normal routine main call 
main()



