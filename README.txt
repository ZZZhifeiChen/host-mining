########workflow########

1, Users should import database downloaded from NCBI or other sources 
into "database" folder. Each species in database requires a ".gff" 
file and ".faa" file. These two files must have the same file name, 
such as:
"GCA_000007085.1_ASM708v1.faa"
"GCA_000007085.1_ASM708v1.gff"

2, Enter the working path set by user. The working path must contain 
four files: 
"screening.py"(the python script), 
"database"(the folder storing ".gff" and ".faa" files of each species), 
"diamond"(Diamond v2.1.4 software of Linux version), 
"input.xlsx"(Excel table storing input sequences information). 
The sample files of “database”, “diamond”, “input.xlsx” are displayed.

3, Run the "screening.py" script using Python. The script workflow is as 
follows: 
   Firstly, use Diamond software to establish the ".dmnd" database files 
of each species, and save them into the newly created "DB" folder.

   Secondly, read the input sequences information from "input.xlsx" file 
as BLASTp query sequences, and save them into the newly created "protein_seqs.fasta" 
file. Then, use Diamond software to perform BLASTp on the ".dmnd" database 
files in the "DB" folder with the sequences in the "protein_seqs.fasta" file. 
Save the BLASTp results in the newly created "blast_output" folder.

   Thirdly, analyze the results in the "blast_output" folder, and score each
 species based on two aspects: "gene number" and "gene distance". Then obtain 
the total score according to these two aspects. Conserve the alignment results 
and scores of each species in the “Dataframe” format data.

   Finaly, Sort the total scores of the "output.xlsx", and output the Top50 
species into the "output.xlsx" result file. Delete the intermediate process 
files such as the "DB" folder, "blast_output", "protein_seqs.fasta", etc.

#################################

########"output.xlsx" file illustration########

The "output.xlsx" file is sorted by the total score of each species and 
contains Top50 species information. 

1, The first column contains "Species ID", which refers to the file name of 
corresponding species in the "database" folder.

2, The second column contains "aligned genes", which refers to the query 
sequences name and the BLASTp aligned protein id of species in the database. 
The query sequences name is displayed between the "|" character, and the BLASTp 
aligned protein ids are displayed at intervals of " ".

3, The third column contains "gene number", which refers to the numbers of 
BLASTp aligned protein.

4, The forth column contains "number score", which refers to the scores of 
species in the "gene number" aspects. If there are n query sequences and the 
number of BLASTp aligned protein of the species is m, the "number score" of 
the species is m/n * 100.

5, The fifth column contains "gene distance", which refers to the minimum 
distance of the genes corresponding to the BLASTp aligned proteins.

6, The sixth column contains "distance score", which refers to the scores of 
species in the "gene distance" aspects. We firstly sort the species according 
to the calculated "gene distance", and score them from minimum to maximum.

7, The seventh column contains "total score", which refers to the final scores 
of each species. The weight of "number score" is 1, and the weight of "distance 
score" is 1/n. Then add them up to get "total score".

#################################

########Note########

1, This script can only be used on computers with Linux operating system. 
The “input.xlsx” can be written on the windows operating system.

2, When selecting the database, it is recommended to select species files with 
complete genome assembly on the NCBI website to ensure the integrity of the 
species genome information.

3, Before using the "screening.py" python script, please make sure that third-party 
libraries such as "os", "subprocess", "shutil", "pandas", "itertools", and 
"biopython" are installed in the Python environment.

4, The cells in the "Protein seq" column of the "input.xlsx" file can contain line 
break character "\n".

###################################