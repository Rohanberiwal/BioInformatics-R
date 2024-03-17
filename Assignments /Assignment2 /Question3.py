from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "rohanberiwal@gamil.com"
accession_list = {
    "NC_001144.5 PDC1": ("NC_001144.5", 232390, 234081),
    "NC_001139.9 PDC5": ("NC_001139.9", 651290, 652981),
    "NC_001144.5 PDC6": ("NC_001144.5", 410723, 412414)
}

def fetch_sequence(accession_number, start, end):
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
    sequence = handle.read()
    handle.close()
    return sequence

def blast_initiator(sequence):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_record = NCBIXML.read(result_handle)
    return blast_record

def find_snps(blast_record):
    snps = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.001: 
                for i, (query_base, hit_base) in enumerate(zip(hsp.query, hsp.sbjct)):
                    if query_base != hit_base:  
                        snp_position = hsp.query_start + i
                        snps.append((alignment.title, snp_position, query_base, hit_base))
    return snps

def write_to_file(filename, gene_name, snps):
    with open(filename, "a") as f:
        f.write(f"Gene Name: {gene_name}\n")
        for alignment_title, position, query_base, hit_base in snps:
            f.write(f"  Alignment Title: {alignment_title}, Position: {position}, Query Base: {query_base}, Hit Base: {hit_base}\n")

def print_information(gene_name, snps):
    print(f"Gene Name: {gene_name}")
    for alignment_title, position, query_base, hit_base in snps:
        print(f"  Alignment Title: {alignment_title}, Position: {position}, Query Base: {query_base}, Hit Base: {hit_base}")

def main():
    total_snps = {}
    output_filename = "SNP_output.txt"
    with open(output_filename, "w") as f:
        f.write("SNPs Detected:\n")
    for gene_name, (accession_number, start, end) in accession_list.items():
        print(f"Processing {gene_name} (Accession Number: {accession_number})...")
        sequence = fetch_sequence(accession_number, start, end)
        blast_record = blast_initiator(sequence)
        snps = find_snps(blast_record)
        total_snps[gene_name] = len(snps)
        write_to_file(output_filename, gene_name, snps)
        print_information(gene_name, snps)
    print("Total SNPs:")
    for gene_name, count in total_snps.items():
        accession_number = accession_list[gene_name][0]
        print(f"{gene_name} (Accession Number: {accession_number}): {count}")

if __name__ == "__main__":
    main()
