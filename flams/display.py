#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Retro212, hannelorelongin, annkamsk
"""

import csv
import logging

""" setup
This script deals with returning the results of FLAMS to the user in a tsv file.
"""

def display_result(output_filename, blast_records):
    """
    This function creates a tsv file containing all conserved modification sites, based on a specific FLAMS run.

    Parameters
    ----------
    output_filename: str
        Output file name
    blast_records: array
        Array containing BLAST records that met search criteria of FLAMS run.

    """
    logging.info(f"Writing .tsv output file with all conserved lysine modifications.")
    with open(output_filename, "w") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        tsv_writer.writerow(
            [
                "Uniprot ID",
                "Modification",
                "Lysine location",
                "Lysine Window",
                "Species",
            ]
        )
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    ## Parsing header of format PLMD_ID | UniProt_ID | lysinePosition modificationType [species]
                    headerSplitPipe = (alignment.title).split("|")  # split up header into list, seperated by pipe
                    plmd_id = headerSplitPipe[0]
                    uniprot_id = headerSplitPipe[1]
                    pos_type_speciesSplitBracket = headerSplitPipe[2].split("[")
                    lysine_location = int(pos_type_speciesSplitBracket[0].split()[0])
                    modification_type = pos_type_speciesSplitBracket[0].split()[1]
                    species = pos_type_speciesSplitBracket[1][:-1]
                    tsv_writer.writerow(
                        [
                            uniprot_id,
                            modification_type,
                            lysine_location,
                            _getSequenceWindow(hsp, lysine_location),
                            species,
                        ]
                    )

def _getSequenceWindow(hsp,lysine_location):
    """
    This function creates generates the sequence window around the modified lysine.
    If the modified lysine is not near the end (neither in the query as well as aligned sequence),
        it simply returns the window containing the 5 amino acids before and after the modified lysine.
    However, if the modified lysine is near either the start or the end of the aligned sequence, the sequence window can only contain part of this window,
        and this function makes sure this limit is respected.

    Parameters
    ----------
    hsp: hsp
        High Scoring partner, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    lysine_location: int
        Position of lysine in the aligned protein that is known to be modified

    """
    sequence = hsp.sbjct.replace("-","")
    protSize = len(sequence)
    modPos = lysine_location - hsp.sbjct_start
    lysineWindowMax = modPos+6
    lysineWindowMin = modPos-5
    if modPos + 6 > protSize:
        lysineWindowMax = protSize
    if modPos - 6 < 0:
        lysineWindowMin = 0
    windowString = (str(lysineWindowMin+hsp.sbjct_start) + "-" + sequence[lysineWindowMin:lysineWindowMax] + "-" + str(lysineWindowMax+hsp.sbjct_start-1))
    return windowString
