from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Record import Blast, Alignment
from Bio.Blast import NCBIXML
import flams.databases.setup
import re


def run_blast(input, modifications, lysine_pos, lysine_range, evalue=0.01, num_threads=1, **kwargs):
    # For each modification, run blast and flatten results to an array
    results = []
    for m in modifications:
        result = _run_blast(input, m, lysine_pos, lysine_range, evalue, num_threads)
        for r in result:
            results.append(r)
    return results


def _run_blast(input, modification, lysine_pos, lysine_range, evalue, num_threads=1):
    # Get BLASTDB path for selected modification + get a temporary path for output
    BLASTDB = flams.databases.setup.get_blastdb_for_modification(modification)
    BLAST_OUT = f'{flams.databases.setup.get_data_dir()}/temp.xml'

    # Run BLAST
    blast_exec = NcbiblastpCommandline(
        query=input,
        db=BLASTDB,
        evalue=evalue,
        outfmt=5,
        out=BLAST_OUT,
        num_threads=num_threads,
    )
    blast_exec()

    with open(BLAST_OUT) as handle:
        blast_records = list(NCBIXML.parse(handle))

    return [_filter_blast(i, lysine_pos, lysine_range, evalue) for i in blast_records]


def _filter_blast(blast_record, lysine_pos, lysine_range, evalue):
    # Create new Blast Record where we append filtered matches.
    filtered = Blast()
    # We will here assume that the input FASTA contains only one sequence (array is of length 1).
    for a in blast_record.alignments:
        # We will append matching High Scoring partners here, which will then be added to the 'filtered' BLAST frame
        filtered_hsps = []
        for hsp in a.hsps:
            if hsp.expect < evalue:

                # Parse FASTA title so we can pick out the PTM position
                title_spl = _parse_fasta_title(a.title)

                mod_pos = int(title_spl[2])

                if not _check_ptm_is_within_match(hsp, mod_pos):
                    continue

                query_pos = lysine_pos

                if not _check_user_query_is_within_match(hsp, query_pos):
                    continue

                mod_pos, query_pos, limit_low, limit_high = _standardise_positions(
                    hsp, mod_pos, query_pos, lysine_range
                )

                # Check 3. Check if mod_pos is within range of low and high
                if limit_low <= mod_pos <= limit_high:
                    # WEE! we have a match.
                    filtered_hsps.append(hsp)

        # If some HSPS matched, let's append that to the filtered BLAST frame for future processing.
        if filtered_hsps:
            new_alignment = Alignment()
            new_alignment.title = a.title
            new_alignment.hsps = filtered_hsps
            filtered.alignments.append(new_alignment)

    # Display results expects an array of BLAST records.
    return filtered


def _parse_fasta_title(title):
    # Get lysine modification position from alignment title.
    # User input is: lysine_pos, lysine_range

    # Parse PTM modification from the alignment title
    # Example: PLMD-7244|P25665|304 Acetylation [Escherichia coli (strain K12)]
    title_spl = re.split(r"\||\s", title, maxsplit=4)

    # title_spl[0] contains PLMD id
    # title_spl[1] contains Uniprot ID
    # title_spl[2] contains modification position
    # title_spl[3] contains type of modification
    # title_spl[4] contains name or organism.
    return title_spl


# Check 1. mod_pos must be within the match of the subject
def _check_ptm_is_within_match(hsp, mod_pos):
    return hsp.sbjct_start <= mod_pos <= hsp.sbjct_end


# Check 2. Our user queried position must be within the match of the query
def _check_user_query_is_within_match(hsp, query_pos):
    return hsp.query_start <= query_pos <= hsp.query_end


def _standardise_positions(hsp, mod_pos, lysine_pos, lysine_range):
    # Standardise the position of the found PTM to local alignment
    mod_pos = mod_pos - (hsp.sbjct_start - 1)

    # Standardise the position of the query to local alignment and set range
    query_pos = lysine_pos - (hsp.query_start - 1)
    limit_low = query_pos - lysine_range
    limit_high = query_pos + lysine_range

    return mod_pos, query_pos, limit_low, limit_high
