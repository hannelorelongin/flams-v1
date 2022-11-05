import csv
from zipfile import ZipFile
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

# Dummy protein modifications database
# Can be used to generate BLASTDB from local FASTA zip files.


def get_fasta(descriptor, location):
    # Read ZIP-compressed FASTA from path designated by descriptor and write to location
    with ZipFile(descriptor, mode="r") as myzip:
        # Extract the single txt file and return as UTF-8 string
        plm = myzip.read(myzip.namelist()[0]).decode("UTF-8")

    with open(location, "a") as out:
        SeqIO.write(_convert_plm_to_fasta(plm), out, "fasta")


def _convert_plm_to_fasta(plm):
    recs = []
    reader = csv.reader(StringIO(plm), delimiter="\t")
    for row in reader:
        seq = Seq(row[6])
        id = f"{row[0]}|{row[1]}|{row[2]}"
        rec = SeqRecord(
            seq,
            id=id,
            description=f"{row[3]} [{row[5]}]",
            annotations={
                "Species": row[5],
                "Uniprot Accession": row[1],
                "Position": row[2],
                "Type": row[3],
                "PMIDs": row[6].split(";"),
            },
        )
        recs.append(rec)
    return recs
