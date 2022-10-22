import csv
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

PLM_DATABASE = "../Acetylation.elm"
OUTPUT = "data/acetylation.faa"

recs = []
with open(PLM_DATABASE, "r+") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        seq = Seq(row[4])
        id = f"{row[0]}|{row[1]}|{row[2]}"
        rec = SeqRecord(seq, id=id, description=f"{row[3]} [{row[5]}]", annotations={"Species": row[5], "Uniprot Accession": row[1], "Position": row[2], "Type": row[3], "PMIDs": row[6].split(";")})
        recs.append(rec)

SeqIO.write(recs, OUTPUT, "fasta")
