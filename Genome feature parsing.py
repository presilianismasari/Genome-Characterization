#Import package yang diperlukan
from Bio import SeqIO
import pandas as pd
from Bio.SeqUtils import gc_fraction
from collections import Counter
import re
import os

#Untuk menghitung jenis dan jumlah types yang ada di seq_record (genebank)
def feature_type_counter(seq_record):
    feature_types = {}
    for feature in seq_record.features:
        if feature.type in feature_types:
            feature_types[feature.type] += 1
        else:
            feature_types[feature.type] = 1
    for key, value in feature_types.items():
        print(f"{key}: {value}")
#Untuk buat file khusus trna SAJA
def count_tRNA_only(gbk_path):
    tRNA_list = []
    for record in SeqIO.parse(gbk_path, "genbank"):
        for feature in record.features:
            if feature.type == "tRNA":
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                tRNA_list.append(product)
    tRNA_counter = Counter(tRNA_list)
    return pd.DataFrame(tRNA_counter.items(), columns=["Jenis tRNA", "Jumlah"])
#Untuk buat file khusus trna rRNA
def count_rRNA_only(gbk_path):
    rRNA_list = []
    for record in SeqIO.parse(gbk_path, "genbank"):
        for feature in record.features:
            if feature.type == "rRNA":
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                rRNA_list.append(product)
    rRNA_counter = Counter(rRNA_list)
    return pd.DataFrame(rRNA_counter.items(), columns=["Jenis rRNA", "Jumlah"])

#Untuk membuat class "sequence" supaya memudahkan analisis
class sequence:
    def __init__(self,seq_id,seq,seq_product,seq_translation,cog,cog_category,start_location,end_location,strand_type,gene_name,gene_product,protein_sequence):
        self.id = seq_id
        self.seq = seq
        self.product =seq_product
        self.translation = seq_translation
        self.cog = cog
        self.cog_category = cog_category
        self.startloc = start_location
        self.endloc = end_location
        self.strand = strand_type
        self.gene = gene_name
        self.product = gene_product
        self.translation = protein_sequence

#Fungsi untuk ekstraksi informasi yang diperlukan di dalam tabel
#Ekstraksi informasi locus tag, produk, dan lokasi --> disimpan ke dalam sequence_list
def info_extraction(seq_record,feature_type):
    sequence_list = []
    for feature in seq_record.features:
        if feature.type == feature_type: 
            feature_id = feature.qualifiers['locus_tag'][0]
            feature_product = feature.qualifiers['product'][0]
            feature_sequence = feature.location.extract(seq_record).seq
            start_location = feature.location.start
            end_location = feature.location.end
            strand_type = feature.location.strand
           
            try:
             feature_translation = feature.qualifiers['translation'][0]
            except:
                feature_translation = "-"
            try:
                gene_name = feature.qualifiers['gene']
            except:
                gene_name = "-"
            try:
                gene_product = feature.qualifiers['product']
            except:
                gene_product = "-"
            try:
                protein_sequence = feature.qualifiers['translation']
            except:
                protein_sequence = "-"
            try:
                if feature.qualifiers['db_xref'][0].split(":")[0] == "COG":
                    feature_cog = feature.qualifiers['db_xref'][0].split(":")[1]
            except:
                feature_cog = "-"

            cog_category = "-"
            try:
                for note in feature.qualifiers.get('note', []):
                    if "COG:" in note and "Category:" in note:
                        match_category = re.search(r"\[Category:([A-Z]+)", note)
                        if match_category:
                            cog_category = match_category.group(1)
            except:
                pass

            sequence_list.append(sequence(feature_id,feature_sequence,feature_product,feature_translation,feature_cog,cog_category,start_location,end_location,strand_type,gene_name,gene_product,protein_sequence))
    return sequence_list

def append_information(sequence_list):
#Membuat list 
    locus_tag_list = []
    start_location_list = []
    end_location_list = []
    strand_type_list = []
    gene_name_list = []
    gene_product_list = []
    nucleotide_sequence_list = []
    protein_sequence_list = []
    cog_category_list = []

#Penambahan informasi ke dalam list
    i = 0
    for sequence in sequence_list:
        locus_tag = sequence_list[i].id
        start_location = sequence_list[i].startloc
        end_location = sequence_list[i].endloc
        strand_type = sequence_list[i].strand
        gene_name = sequence_list[i].gene
        gene_product = sequence_list[i].product
        nucleotide_sequence = sequence_list[i].seq
        protein_sequence = sequence_list[i].translation
        cog_category = sequence_list[i].cog_category
#print("====================")
#print(locus_tag)
#print(start_location)
#print(end_location)
#print(strand_type)
#print("====================")

        locus_tag_list.append(locus_tag)
        start_location_list.append(start_location)
        end_location_list.append(end_location)
        strand_type_list.append(strand_type)
        gene_name_list.append(gene_name)
        gene_product_list.append(gene_product)
        nucleotide_sequence_list.append(nucleotide_sequence)
        protein_sequence_list.append(protein_sequence)
        cog_category_list.append(cog_category) 
        i += 1
    return locus_tag_list, start_location_list, end_location_list, strand_type_list,gene_name_list,gene_product_list,nucleotide_sequence_list,protein_sequence_list,cog_category_list

##==Script untuk ambil sekuens 16S rRNA dari (.gbk)==
def extract_16s_rRNA_from_genbanks(gbk_file, output_fasta):
    count = 0
    if isinstance(gbk_file, str):
        gbk_file = [gbk_file]
    with open(output_fasta, "w") as fasta_out:
        for gb_file in gbk_file:
            if not os.path.exists(gb_file):
                print(f"File not found: {gb_file}")
                continue
            for record in SeqIO.parse(gb_file, "genbank"):
                for i, feature in enumerate(record.features):
                    if feature.type == "rRNA":
                        product = feature.qualifiers.get("product", [""])[0].lower()
                        if "16s" in product and "ribosomal rna" in product:
                            seq = feature.extract(record.seq)
                            start = int(feature.location.start) + 1
                            end = int(feature.location.end)
                            header = f"{record.id}_16S_{i+1}_{start}-{end}"
                            fasta_out.write(f">{header}\n{seq}\n")
                            count += 1
    print(f"Extracted {count} 16S rRNA sequences to '{output_fasta}'.")


#Membuka file genbank (.gbk)
gbk_file = "C:/Users/User/Documents/TA2/karakterisasi/all/Planococcus_Selat_Makassar.genbank"
seq_record_list = []
for seq_record in SeqIO.parse(gbk_file, "genbank"):
    seq_record_list.append(seq_record)
#Menghitung jumlah contig yang ada di dalam file genbank (.gbk)
number_of_contig_sequence = len(seq_record_list)
print("Jumlah contig:",number_of_contig_sequence)

#Membuat daftar jenis "type" apa saja yang berada di dalam genom beserta jumlahnya
for seq_record in seq_record_list:
    feature_type_counter(seq_record)
    print("-----")


#Ekstraksi informasi terkait "type" yang diinginkan
feature_type = 'CDS'
final_locus_tag_list = []
final_start_location_list = []
final_end_location_list = []
final_strand_type_list = []
final_gene_name_list = []
final_gene_product_list = []
final_nucleotide_sequence_list = []
final_protein_sequence_list = []
final_cog_category_list = []

for seq_record in seq_record_list:
    sequence_list = info_extraction(seq_record,feature_type)
    locus_tag_list, start_location_list, end_location_list, strand_type_list, gene_name_list, gene_product_list, nucleotide_sequence_list, protein_sequence_list, cog_category_list = append_information(sequence_list)
    final_locus_tag_list += locus_tag_list
    final_start_location_list += start_location_list
    final_end_location_list += end_location_list
    final_strand_type_list += strand_type_list
    final_gene_name_list += gene_name_list
    final_gene_product_list += gene_product_list
    final_nucleotide_sequence_list += nucleotide_sequence_list
    final_protein_sequence_list += protein_sequence_list
    final_cog_category_list += cog_category_list

#Membuat data frame hasil ekstraksi informasi "type" yang diinginkan
gbk_df = pd.DataFrame()
name_column = pd.Series(final_locus_tag_list, name="name")
start_column = pd.Series(final_start_location_list, name="start")
stop_column = pd.Series(final_end_location_list, name="stop")
strand_column = pd.Series(final_strand_type_list, name="strand")
gene_column = pd.Series(final_gene_name_list, name="gene")
product_column = pd.Series(final_gene_product_list, name="product")
sequence_column = pd.Series(final_nucleotide_sequence_list, name="sequence")
translation_column = pd.Series(final_protein_sequence_list, name="translation")
cog_category_column = pd.Series(final_cog_category_list, name="COG_category")
c1,c2,c3,c4,c5,c6,c7,c8,c9 = name_column.to_frame(),start_column.to_frame(),stop_column.to_frame(),strand_column.to_frame(),gene_column.to_frame(),product_column.to_frame(),sequence_column.to_frame(),translation_column.to_frame(),cog_category_column.to_frame()

gbk_df = pd.concat([c1,c2,c3,c4,c5,c6,c7,c8,c9], axis=1)
gbk_df['strand'] = gbk_df['strand'].replace(1, '+')
gbk_df['strand'] = gbk_df['strand'].replace(-1, '-')
gbk_df['start'] += 1
gbk_df

#Menyimpan hasil parsing ke dalam .csv
gbk_df.to_csv("Planococcus_Selat_Makassar_general.csv",index=False)
# Simpan data tRNA saja
tRNA_df = count_tRNA_only(gbk_file)
tRNA_df.to_csv("Planococcus_Selat_Makassar_tRNA.csv", index=False)
#Data rRNA saja
rRNA_df = count_rRNA_only(gbk_file)
rRNA_df.to_csv("Planococcus_Selat_Makassar_rRNA.csv", index=False)
#Ekstrak 16S rRNA
output_fasta = "Sekuens16S.fasta"
extract_16s_rRNA_from_genbanks(gbk_file, output_fasta)


#Menghitung panjang genom dan GC percentage
genome_sequence = ''
for i in range(len(seq_record_list)):
    genome_sequence+=seq_record_list[i].seq
genome_length = len(genome_sequence)
genome_GCpercent = gc_fraction(genome_sequence)
print("Ukuran genom:",genome_length)
print("GC-percentage:",genome_GCpercent)

# Load the original CSV file
df = pd.read_csv("Planococcus_Selat_Makassar_general.csv")

# Remove the "sequence" column if it exists
if "sequence" in df.columns:
    df.drop(columns=["sequence"], inplace=True)

if "translation" in df.columns:
    df.drop(columns=["translation"], inplace=True)

if "product" in df.columns:
    df.drop(columns=["product"], inplace=True)

# Rename "COG_category" column to "type" if it exists
if "COG_category" in df.columns:
    df.rename(columns={"COG_category": "type"}, inplace=True)

# Save to a new CSV file
output_path = "Planococcus_Selat_Makassar_proksee.csv"
df.to_csv(output_path, index=False)

print("File has been modified and saved successfully.")

