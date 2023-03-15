import os
import subprocess
import shutil
import pandas as pd
import itertools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

###############PART 1 Database establishment###############

print("Database establishing ...")
os.mkdir("DB")
for filename in os.listdir("database"):
    if filename.endswith(".faa"):
        db_name = os.path.splitext(filename)[0]
        db_path = os.path.join("DB", db_name)
        data_path = os.path.join("database",filename)
        cmd = ["./diamond", "makedb", "--in", data_path, "--db", db_path, "--quiet"]
        subprocess.run(cmd)
print("Database establishment finished")

###############PART 2 Query sequences process###############

print("Query sequences processing ...")
df = pd.read_excel("input.xlsx", header=0, usecols=[1, 2, 3], names=["ID", "Name", "Protein seq"], dtype={"Protein seq": str})
seq_records = []
for index, row in df.iterrows():
    gene_id = row["ID"]
    gene_name = row["Name"]
    protein_seq = Seq(row["Protein seq"].replace("\n",""))
    seq_record = SeqRecord(protein_seq, id=gene_id, name=gene_name, description="")
    seq_records.append(seq_record)
with open("protein_seqs.fasta", "w") as file:
    SeqIO.write(seq_records, file, "fasta")
print("Query sequences process finished")

###############PART 3 BLASTp alignment###############

print("BLASTp aligning")
if os.path.exists("blast_output"):
    shutil.rmtree("blast_output")
os.mkdir("blast_output")
for filename in os.listdir("DB"):
    out_file = os.path.join("blast_output", os.path.splitext(filename)[0])
    cmd = ["./diamond", "blastp", "-q", "protein_seqs.fasta", "--quiet", "-d", os.path.join("DB",filename), "-o", out_file, "-e", "1e-5", "-f", "6", "qseqid", "sseqid", "-k", "1000000"]
    subprocess.run(cmd)
print("BLASTp alignment finished")

###############PART 4 Result analysis###############
print("Result analysing")
def middle(size, start, end):
    if end > start:
        length = end - start
        middle = (end + start)/2
        return middle, length
    else:               
        left = size - start + 1
        right = end
        length = left + right
        if left > right:
            middle = start + length/2
            return middle, length
        else:
            middle = end - length/2
            return middle, length   
        
def distance(ls_middle, ls_length, size, topology):    
    if topology == "linear":
        size = 1e12 
    diss = []
    for i in range(len(ls_middle)-1):
        dis = ls_middle[i+1] - ls_middle[i] - (ls_length[i+1] + ls_length[i]) / 2
        diss.append(dis)
    final_dis = (size - ls_middle[len(ls_middle)-1] + ls_middle[0]) - (ls_length[len(ls_middle)-1] + ls_length[0])/ 2
    diss.append(final_dis)
    diss.remove(max(diss))
    distance = sum(diss) / len(diss)
    return round(distance, 2)    

output_columns = ["Species ID", "aligned genes", "gene number", "number score", "distance"]
output = pd.DataFrame(columns = output_columns)
Species_IDs = []
aligned_genes = []
gene_numbers = []
number_scores = []
distances = []
for filename in os.listdir("blast_output"):
    blast_file = os.path.join("blast_output",filename)
    if os.path.getsize(blast_file) > 0:
        Species_IDs.append(filename)
        df1 = pd.read_csv(blast_file, sep='\t', header=None, usecols=[0, 1], names=["query ID", "gene ID"], dtype={"query ID": str, "gene ID": str})
        gff_file = os.path.join("database",filename + ".gff")
        df2 = pd.read_csv(gff_file, sep='\t', header=None, comment="#", usecols=[0,3,4,6,8], names=["seqid", "start", "end", "strand", "attributes"], dtype={"strand": str, "attributes": str})
        seqids = []
        query_names = []
        starts = []
        ends = []
        middles = []
        lengths = []
        strands = []
        sizes = []
        topologys = []
        for i in range(len(df1)):
            query_id = df1.loc[i, "query ID"]
            gene_id = df1.loc[i, "gene ID"]
            query_row = df.loc[df["ID"] == query_id]
            query_name = query_row["Name"]
            df2_row = df2.loc[df2["attributes"].str.contains(gene_id)]
            seqid = df2_row["seqid"].values[0]  
            size = df2.loc[df2["seqid"] == seqid].iloc[0]["end"]
            info = df2.loc[df2["seqid"] == seqid].iloc[0]["attributes"] 
            attributes = df2.loc[df2["seqid"] == seqid].iloc[0]["attributes"]
            if "circular" in info:  
                topologys.append("circular")
            else:
                topologys.append("linear")
            seqids.append(seqid)
            query_names.append(query_name.values[0])
            starts.append(df2_row["start"].values[0])
            ends.append(df2_row["end"].values[0])
            middles.append(middle(size, starts[i], ends[i])[0])
            lengths.append(middle(size, starts[i], ends[i])[1])
            strands.append(df2_row["strand"].values[0])
            sizes.append(size)
        df1.insert(loc=0, column="seqid", value=seqids)
        df1.insert(loc=1, column="query name", value=query_names)
        df1.insert(loc=4, column="start", value=starts)
        df1.insert(loc=5, column="end", value=ends)
        df1.insert(loc=6, column="middle", value=middles)
        df1.insert(loc=7, column="length", value=lengths)
        df1.insert(loc=8, column="strand", value=strands)
        df1.insert(loc=9, column="size", value=sizes)
        df1.insert(loc=10, column="topology", value=topologys)
        ls1 = []
        ls2 = []
        ls3 = []   
        for i in range(len(df1)):
            if df1.loc[i,"query name"] not in ls1:
                ls1.append(df1.loc[i,"query name"])
                ls2.append([df1.loc[i,"gene ID"]])
            else:
                ls2[-1].append(df1.loc[i,"gene ID"])
        for i in range(len(ls1)):
            str1 = "|" + ls1[i] + "|" + " ".join(ls2[i])
            ls3.append(str1)
        gene = "\n".join(ls3)
        gene_number = len(ls1)
        number_score = round(gene_number / len(df) * 100, 2)
        aligned_genes.append(gene)
        gene_numbers.append(gene_number)
        number_scores.append(number_score)
        sub_distance = []
        for n in itertools.product(*ls2):
            sub_ids = list(n)
            sub_df1 = df1[df1["gene ID"].isin(sub_ids)][["seqid", "gene ID", "middle", "length", "size", "topology"]]
            sub_df1_dict = dict(tuple(sub_df1.groupby("seqid")))
            diss = []
            for key in sub_df1_dict.keys():
                sub_df1_sub = sub_df1_dict[key] 
                sub_df1_sub.drop_duplicates(inplace=True)
                sub_df1_sub = sub_df1_sub.sort_values(by="middle", ascending=True)
                if len(sub_df1_sub) > 1:
                    ls_middle = list(sub_df1_sub["middle"])
                    ls_length = list(sub_df1_sub["length"])
                    size = sub_df1_sub.iloc[0]["size"]
                    topology = sub_df1_sub.iloc[0]["topology"]
                    dis = distance(ls_middle, ls_length, size, topology)
                    diss.append(dis)
            N = len(sub_df1_dict)-1
            sub_distance.append(sum(diss) + N*6e6)
        if sum(sub_distance) != 0:
            distances.append(min(sub_distance))
        else:
            distances.append(1e12)
output["Species ID"] = Species_IDs
output["aligned genes"] = aligned_genes
output["gene number"] = gene_numbers
output["number score"] = number_scores
output["distance"] = distances
output1 = output.loc[output["distance"] == 1e12]
output2 = output.loc[output["distance"] != 1e12]
output2 = output2.sort_values(by="distance")
score1 = []
score2 = []
op1_rows = len(output1)
for i in range(op1_rows):
    score1.append(0)
output1.insert(loc=5, column="distance score", value=score1)
op2_rows = len(output2)
for i in range(op2_rows):
    percentage = i / op2_rows
    score = 100 - int(percentage*100)
    score2.append(score)
output2.insert(loc=5, column="distance score", value=score2)
output = pd.concat([output2, output1], axis=0, ignore_index=True)
x = len(df)
number_percentage = 1
distance_percentage = 1/x
total_scores = []
for i in range(len(output)):
    row = output.loc[i]
    total_score = round((row["number score"]*number_percentage + row["distance score"]*distance_percentage) / (1 + 1/x), 3)
    total_scores.append(total_score)
output.insert(loc=6, column="total score", value=total_scores)
output = output.sort_values(by="total score", ascending=False)
output= output.head(50)
output.to_excel('output.xlsx', index=False)
os.remove("protein_seqs.fasta")
shutil.rmtree("blast_output", ignore_errors=True)
shutil.rmtree("DB", ignore_errors=True)
print("Result analysis finished")
print("all works done")

