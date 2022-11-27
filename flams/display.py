import csv


def display_result(output_filename, blast_records):
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec129
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
                    ]  # lysine window islaoted (Range set to 5)
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
