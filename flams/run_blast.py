#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: annkamsk, hannelorelongin, kasgel, MaartenLangen
"""

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Record import Blast, Alignment

from .databases import setup as db_setup
from .utils import get_data_dir

""" run_blast
This script contains all functions necessary to search through the proteins stored in the CPLM database, with BLAST,
and to retrieve those that contain conserved lysine modifications.
"""

def run_blast(
    input,
    modifications,
    lysine_pos,
    lysine_range=0,
    evalue=0.01,
    num_threads=1,
    **kwargs,
):
    """
    This function runs the BLAST search and the following filter steps for each modification.
    Ultimately, it only returns conserved (within the specified range) protein modifications for similar proteins.
    It flattens the results to an array.

    Parameters
    ----------
    input: fasta
        Sequence file of query protein
    modifications: str
        Space-seperated list of modifications (which are keys to any of the ModificationType's stored in the MODIFICATIONS dictionary)
    lysine_pos: int
        Position of lysine in query that is under investigation for conservation
    lysine_range: int (default: 0)
        Error margin for conservation of lysine_pos
    evalue: float (default: 0.01)
        BLAST parameter, e-value for BLAST run
    num_threads: int (default: 1)
        BLAST parameter, number of threads that can be used by BLAST

    """
    results = []
    input = input.absolute()
    for m in modifications:
        result = _run_blast(input, m, lysine_pos, lysine_range, evalue, num_threads)
        for r in result:
            results.append(r)
    return results


@dataclass
class ModificationHeader:
    """
    This dataclass consists of the different components contained in the header of each modification entry, and a function to parse it.

    Parameters
    ----------
    plmd_id: str
        PLMD ID for each modification
    uniprot_id: str
        UniProt ID for protein containing the modification
    position: int
        Position at which the modification was detected
    modification: str
        Post-translational modification found at $position in protein with $uniprot_id
    species: str
        Species that encodes the protein containing the modification

    """
    plmd_id: str
    uniprot_id: str
    position: int
    modification: str
    species: str

    @staticmethod
    def parse(title: str) -> "ModificationHeader":

        regex = (
            r"(?P<plmd_id>\S+)\|"
            r"(?P<uniprot_id>\S+)\|"
            r"(?P<position>\d+) (?P<modification>[A-Za-z1-9-]+) \[(?P<species>.+)\]"
        )
        vars = re.match(regex, title).groupdict()
        vars["position"] = int(vars["position"])
        return ModificationHeader(**vars)


def _run_blast(input, modification, lysine_pos, lysine_range, evalue, num_threads=1):
    """
    This function runs the BLAST search and the following filter steps for 1 modification.
    Ultimately, it only returns conserved (within the specified range) protein modifications for similar proteins.

    Parameters
    ----------
    input: fasta
        Sequence file of query protein
    modification: str
        Modification for which you search (which is the key to any of the ModificationType's stored in the MODIFICATIONS dictionary)
    lysine_pos: int
        Position of lysine in query that is under investigation for conservation
    lysine_range: int
        Error margin for conservation of lysine_pos
    evalue: float
        BLAST parameter, e-value for BLAST run
    num_threads: int (default: 1)
        BLAST parameter, number of threads that can be used by BLAST

    """
    # Get BLASTDB name for selected modification + get a temporary path for output
    BLASTDB = db_setup.get_blastdb_name_for_modification(modification)
    BLAST_OUT = "temp.xml"

    # Adjust working directory conditions
    os.chdir(get_data_dir())

    logging.info(f"Running BLAST search for {input} against local {modification} BLAST database.")
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

    logging.info(f"Filtering results of BLAST search for {input} against local {modification} BLAST database.")
    return [_filter_blast(i, lysine_pos, lysine_range, evalue) for i in blast_records]


def _filter_blast(blast_record, lysine_position, lysine_range, evalue) -> Blast:
    """
    This function filters the BLAST results.
    First, it filters out any BLAST results where the alignment does not contain:
    - the queried lysine_pos (in the protein query)
    - the modification position (in the aligned protein)
    Then, it filters out results where the queried modified lysine does not align with the modified lysine in the aligned protein.
    Ultimately, it only returns a BLAST record containing conserved (within range) protein modifications for similar proteins.

    Parameters
    ----------
    blast_record: Blast
        Blast record containing all similar proteins to the queried one, that are in the specific modification database
    lysine_pos: int
        Position of lysine in query that is under investigation for conservation
    lysine_range: int
        Error margin for conservation of lysine_pos
    evalue: float
        BLAST parameter, e-value for BLAST run

    """
    # Create new Blast Record where we append filtered matches.
    filtered = Blast()

    for a in blast_record.alignments:
        # Parse FASTA title where post-translational modification info is stored
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
    """
    This function asserts that the aligned hit does contain its modification in the aligned portion of the protein.

    Parameters
    ----------
    hsp: hsp
        High Scoring partners, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    mod_pos: int
        Position of lysine in the aligned protein that is known to be modified

    """
    return hsp.sbjct_start <= mod_pos <= hsp.sbjct_end

def _is_modQuery_in_alignment(hsp, query_pos) -> bool:
    """
    This function asserts that the aligned portion of the query protein contains the modification being queried.

    Parameters
    ----------
    hsp: hsp
        High Scoring partners, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    query_pos: int
        Position of lysine in query that is under investigation for conservation

    """
    return hsp.query_start <= query_pos <= hsp.query_end

def _findKs_in_alignedHit(hsp, lysine_position, lysine_range):
    """
    This function finds the relative positions of Ks in the neighbourhood of the position of the residue aligned to the lysine being queried.
    It returns a list of relative positions, all within the lysine_range.

    Parameters
    ----------
    hsp: hsp
        High Scoring partners, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    lysine_pos: int
        Position of lysine in query that is under investigation for conservation
    lysine_range: int
        Error margin for conservation of lysine_pos

    """
    rangeK = []
    for i in range(-lysine_range, lysine_range + 1):
        ##need to check that we do not try to access an index out of range for this subject
        if (lysine_position - hsp.query_start + i <= len(hsp.sbjct) - 1) and (lysine_position - hsp.query_start + i >= 0):
            if hsp.sbjct[lysine_position - hsp.query_start + i] == "K":
                rangeK.append(i)
    return rangeK

def _add_conservedModK_to_listConsHsp(hsp, lysine_pos, lysine_range, modification, listHsp):
    """
    This function adds the hsps of modification database entries with conserved modified lysines to a list, namely the listHsp.

    Parameters
    ----------
    hsp: hsp
        High Scoring partners, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    lysine_pos: int
        Position of lysine in query that is under investigation for conservation
    lysine_range: int
        Error margin for conservation of lysine_pos
    modification: ModificationHeader
        Modification for which you search
    listHsp: list
        List that will be used to append hsps of modification database entries with conserved modified lysines to

    """
    for i in _findKs_in_alignedHit(hsp, lysine_pos, lysine_range):
        indexKhit = lysine_pos - hsp.query_start + i
        numGapUntilK = hsp.sbjct[0:lysine_pos - hsp.query_start + i].count('-')
        coordKOriginalSubject = indexKhit - numGapUntilK + hsp.sbjct_start
        if modification.position == coordKOriginalSubject:
            listHsp.append(hsp)
