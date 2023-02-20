#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kasgel, hannelorelongin, annkamsk
"""

from dataclasses import dataclass
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Record import Blast, Alignment
from Bio.Blast import NCBIXML
from flams.utils import get_data_dir
from pathlib import Path
import flams.databases.setup
import re
import os


def run_blast(
    input,
    modifications,
    lysine_pos,
    lysine_range=0,
    evalue=0.01,
    num_threads=1,
    **kwargs,
):
    # For each modification, run blast and flatten results to an array
    results = []
    input = input.absolute()
    for m in modifications:
        result = _run_blast(input, m, lysine_pos, lysine_range, evalue, num_threads)
        for r in result:
            results.append(r)
    return results


@dataclass
class ModificationHeader:
    plmd_id: str
    uniprot_id: str
    position: int
    modification: str
    species: str

    @staticmethod
    def parse(title: str) -> "ModificationHeader":
        # Parse modification from the alignment title of structure:
        # {PLMD id}|{Uniprot ID}|{modification position} {type of modification} [{species}]
        # Example: PLMD-7244|P25665|304 Acetylation [Escherichia coli (strain K12)]
        regex = (
            r"(?P<plmd_id>\S+)\|"
            r"(?P<uniprot_id>\S+)\|"
            r"(?P<position>\d+) (?P<modification>[A-Za-z1-9-]+) \[(?P<species>.+)\]"
        )
        vars = re.match(regex, title).groupdict()
        vars["position"] = int(vars["position"])
        return ModificationHeader(**vars)


def _run_blast(input, modification, lysine_pos, lysine_range, evalue, num_threads=1):
    # Get BLASTDB name for selected modification + get a temporary path for output
    BLASTDB = flams.databases.setup.get_blastdb_name_for_modification(modification)
    BLAST_OUT = "temp.xml"

    # Adjust working directory conditions and convert input file into absolute path
    os.chdir(get_data_dir())

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


def _filter_blast(blast_record, lysine_position, lysine_range, evalue) -> Blast:
    # Create new Blast Record where we append filtered matches.

    filtered = Blast()

    for a in blast_record.alignments:
        # Parse FASTA title where Posttranslational modification info is stored
        mod = ModificationHeader.parse(a.title)

        # Append matching High Scoring partners here, which will then be added to the 'filtered' BLAST frame
        filter1_hsps = [] ## Filter1: filters out all hsps which do not contain the modification (both in query and hit)
        filter2_hsps = [] ## Filter2: filters out hsps that do not contain CONSERVED modification

        for hsp in a.hsps:
            if hsp.expect < evalue and _is_modHit_in_alignment(hsp, mod.position) and _is_modQuery_in_alignment(hsp, lysine_position):
                # WEE! we have a match.
                filter1_hsps.append(hsp)

        for hsp in filter1_hsps:
            # To assess whether a hsp contains a conserved modification, we need to
            # (1) find the location of the query modification in the aligned query
            if hsp.query.find('-') == -1:
            # (2) find out if the aligned position (+- range) in the hit is a lysine
                if len(_findKs_in_alignedHit(hsp, lysine_position, lysine_range)) != 0:
            # (3) if this aligned position is a lysine, was this the lysine carrying the modification
                    _add_conservedModK_to_listConsHsp(hsp, lysine_position, lysine_range, mod, filter2_hsps)
            # (1) find the location of the query modification in the aligned query
            elif (hsp.query_start + hsp.query.find('-') + 1) > lysine_position:
            # (2) find out if the aligned position (+- range) in the hit is a lysine
                if len(_findKs_in_alignedHit(hsp, lysine_position, lysine_range)) != 0:
            # (3) if this aligned position is a lysine, was this the lysine carrying the modification
                    _add_conservedModK_to_listConsHsp(hsp, lysine_position, lysine_range, mod, filter2_hsps)
            # (1) find the location of the query modification in the aligned query
            else:
            #    should adapt lysine position here to match number of gaps before
                countGapBefore = hsp.query[0:lysine_position+1].count("-")
                newSeq = hsp.query[0:lysine_position+1].replace("-","") + hsp.query[lysine_position+1:len(hsp.query)]
                while newSeq[0:lysine_position+1].find('-') != -1:
                    newSeq = newSeq[0:lysine_position+1].replace("-","") + newSeq[lysine_position+1:len(newSeq)]
                    countGapBefore += 1
            # (2) find out if the aligned position (+- range) in the hit is a lysine
                if len(_findKs_in_alignedHit(hsp, lysine_position + countGapBefore, lysine_range))  != 0:
            # (3) if this aligned position is a lysine, was this the lysine carrying the modification
                    _add_conservedModK_to_listConsHsp(hsp, lysine_position + countGapBefore, lysine_range, mod, filter2_hsps)

        # If some HSPS matched, let's append that to the filtered BLAST frame for future processing.
        if filter2_hsps:
            new_alignment = Alignment()
            new_alignment.title = a.title
            new_alignment.hsps = filter2_hsps
            filtered.alignments.append(new_alignment)

    return filtered

def _is_modHit_in_alignment(hsp, mod_pos) -> bool:
    return hsp.sbjct_start <= mod_pos <= hsp.sbjct_end

def _is_modQuery_in_alignment(hsp, query_pos) -> bool:
    return hsp.query_start <= query_pos <= hsp.query_end

def _findKs_in_alignedHit(hsp, lysine_position, lysine_range):
    rangeK = []
    for i in range(-lysine_range, lysine_range + 1):
        ##need to check that we do not try to access an index out of range for this subject
        if (lysine_position - hsp.query_start + i <= len(hsp.sbjct) - 1) and (lysine_position - hsp.query_start + i >= 0):
            if hsp.sbjct[lysine_position - hsp.query_start + i] == "K":
                rangeK.append(i)
    return rangeK

def _add_conservedModK_to_listConsHsp(hsp, lysine_pos, lysine_range, modification, listHsp):
    for i in _findKs_in_alignedHit(hsp, lysine_pos, lysine_range):
        indexKhit = lysine_pos - hsp.query_start + i
        numGapUntilK = hsp.sbjct[0:lysine_pos - hsp.query_start + i].count('-')
        coordKOriginalSubject = indexKhit - numGapUntilK + hsp.sbjct_start
        if modification.position == coordKOriginalSubject:
            listHsp.append(hsp)
