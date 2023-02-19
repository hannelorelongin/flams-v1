#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Retro212, annkamsk
"""

import csv

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
    with open(output_filename, "w") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        tsv_writer.writerow(
            [
                "Uniprot ID",
                "Species",
                "Modification",
                "Lysine location",
                "Lysine Window",
            ]
        )
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    pubmedID = (alignment.title).split("|")  # from title pubmed_Id
                    Species = pubmedID[2].split("[")  # isolates species name
                    lysine = Species[0]
                    lysine_no = lysine.split()  # gets lysine modification
                    Species_name = Species[1].split("]")
                    lysine_location = int(
                        lysine_no[0]
                    )  # gets lysine no from the title alignment
                    no_dashes = hsp.sbjct.replace(
                        "-", ""
                    )  # dashes in alignment are replaced for getting lysine window
                    mod_pos = lysine_location - (
                        hsp.sbjct_start
                    )  # lysine location calculated relative to blast file
                    lysine_window = no_dashes[
                        mod_pos - 5 : mod_pos + 6
                    ]  # lysine window isolated (Range set to 5)
                    start = str(mod_pos - 5 + hsp.sbjct_start)
                    end = str(mod_pos + 6 + hsp.sbjct_start)
                    window = (
                        start + "-" + lysine_window + "-" + end
                    )  # displaying the exact window location within the sequence (gaps are not counted)
                    tsv_writer.writerow(
                        [
                            pubmedID[1],
                            Species_name[0],
                            lysine_no[1],
                            lysine_location,
                            window,
                        ]
                    )
