import csv
from io import BytesIO, StringIO
from zipfile import ZipFile
from urllib.request import urlopen
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

URL = "http://cplm.biocuckoo.cn/Download/{0}.zip"


def get_fasta(label, location):
    zip = urlopen(URL.format(label))
    size_in_mb = int(zip.getheader('content-length')) / 1048576
    print("Downloading CPMLv4 {0} Database, please wait. Size: {1:.1f} MB".format(label, size_in_mb))

    myzip = ZipFile(BytesIO(zip.read()))

    # Extract the single txt file and return as UTF-8 string
    plm = myzip.read(myzip.namelist()[0]).decode('UTF-8')

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