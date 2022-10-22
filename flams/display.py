import csv


def display_result(blast_records):
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec129
    with open(r"D:\git\flams\output.tsv", "w") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        tsv_writer.writerow(["ID", "length", "e-value"])
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    tsv_writer.writerow(
                        [alignment.title, alignment.length, hsp.expect]
                    )
