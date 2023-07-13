#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: hannelorelongin, kasgel
"""

import csv
import logging
import requests
from io import BytesIO, StringIO
from zipfile import ZipFile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


""" cplmv4
This script downloads the different contents of the CPLM database, and saves them under the FASTA format.
Script developed to work with CPLM database version 4.
"""

URL = "http://cplm.biocuckoo.cn/Download/{0}.zip"

def get_fasta(descriptor, location):
    """
    This function downloads the entries of the CPLM database for a modification type (specified by $descriptor),
    and saves this information as a multi-FASTA file, with one sequence record per modification, in $location.

    Parameters
    ----------
    descriptor: str
        Label that identifies the modification type,
        reflecting the spelling of the modification on the CPLM website
        (must be one from .setup.ModificationDatabase)
    location: str
        Name of output file

    """
    # HTTP request with stream. This way, we get the size of the file first and can begin downloading it in chunks.
    req = requests.get(URL.format(descriptor), stream=True)

    # Raise an exception if HTTP request failed.
    req.raise_for_status()

    size_in_mb = int(req.headers.get("content-length")) / 1048576

    logging.info(f"Downloading CPMLv4 {descriptor} Database, please wait. Size: {size_in_mb:.1f} MB")

    with ZipFile(BytesIO(req.content)) as myzip:
        # Extract the single txt file and return as UTF-8 string
        plm = myzip.read(myzip.namelist()[0]).decode("UTF-8")
        # SeqIO can not write greek letter beta, immediately change this
        if descriptor == "β-Hydroxybutyrylation":
            plm = plm.replace('β','beta')

    with open(location, "a") as out:
        SeqIO.write(_convert_plm_to_fasta(plm), out, "fasta")

    logging.info(f"Converted and stored CPMLv4 {descriptor} Database entries as FASTA entries for the local {descriptor} BLAST database format.")


def _convert_plm_to_fasta(plm):
    """
    This function converts the string containing all entries of the CPLM database for a specific modification
    to a list of sequence records of the FASTA format (one sequence record per modification).

    Parameters
    ----------
    plm: str
        Content of the text file detailing all entries of the CPLM database for a specific modification

    """
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
