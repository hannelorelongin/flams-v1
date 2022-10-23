import csv


def display_result(blast_records):
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec129
    with open(r"output.tsv", "w") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        tsv_writer.writerow(["Uniprot ID", "Species", "lysine location"])
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    pubmedID = (alignment.title).split("|")
                    Species_name = pubmedID[4].split("[")
                    lysine = Species_name[0]
                    lysine_no = lysine.split()
                    Species_name = Species_name[1].split("]")
                    lysine_location = int(lysine_no[0])
                    # print(lysine_location)
                    # print(alignment.title)
                    # hit_seq= [x for x in hsp.sbjct]
                    # mod_pos= lysine_location - (hsp.sbjct_start-1) +1
                    # if lysine_location <= hsp.sbjct_end and lysine_location >= hsp.sbjct_start:
                    #     print(hsp.sbjct[mod_pos-5:mod_pos+5])
                    tsv_writer.writerow([pubmedID[3], Species_name, lysine_location])
