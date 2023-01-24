import csv
import requests
from zipfile import ZipFile
from io import BytesIO, StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

# CPLM Database
# Compendium of Protein Lysine Modifications, Version 4.0
# Database information available at http://cplm.biocuckoo.cn/

# Fallback to data files hosted by ourselves because the original CPLM host has issues
URL = "http://cplm.biocuckoo.cn/Download/{0}.zip"
# URL_working = "https://kasgel.fi/flams/dbs/CPLMv4/{0}.zip"

## Testing if we can switch back to using CPLM download links


def get_fasta(descriptor, location):
    # HTTP request with stream. This way, we get the size of the file first and can begin downloading it in chunks.
    req = requests.get(URL.format(descriptor), stream=True)

    # Raise an exception if HTTP request failed.
    req.raise_for_status()

    size_in_mb = int(req.headers.get("content-length")) / 1048576
    print(
        f"Downloading CPMLv4 {descriptor} Database, please wait. Size: {size_in_mb:.1f} MB"
    )

    with ZipFile(BytesIO(req.content)) as myzip:
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
