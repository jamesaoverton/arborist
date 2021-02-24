#!/usr/bin/env python3

import csv
import sqlite3
import sys

from argparse import ArgumentParser
from collections import defaultdict
from helpers import copy_database, get_curie


EUKARYOTES = ["NCBITaxon:4751", "NCBITaxon:33208", "NCBITaxon:58024"]
METAZOANS = [
    "NCBITaxon:6854",
    "NCBITaxon:6657",
    "NCBITaxon:50557",
    "NCBITaxon:6447",
    "NCBITaxon:6231",
    "NCBITaxon:6157",
    "NCBITaxon:7742",
]


def get_all_ancestors(child_parent, node, limit, ancestors=None):
    if not ancestors:
        ancestors = []
    if node in child_parent:
        parent = child_parent[node]
        if parent == limit:
            return ancestors
        ancestors.append(parent)
        return get_all_ancestors(child_parent, parent, limit, ancestors=ancestors)
    return ancestors


def get_top_ancestor(child_parent, node, limit):
    parent = child_parent[node]
    if parent == limit:
        return node
    return get_top_ancestor(child_parent, parent, limit)


def get_descendants_and_ranks(cur, child_parent, ranks, node):
    # Get the rank of this node
    cur.execute(
        "SELECT object FROM statements WHERE predicate = 'ncbitaxon:has_rank' AND subject = ?",
        (node,),
    )
    res = cur.fetchone()
    if res:
        ranks[node] = res[0]
    # Get the children and maybe iterate
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE object = ? AND predicate = 'rdfs:subClassOf'""",
        (node,),
    )
    for row in cur.fetchall():
        child_parent[row[0]] = node
        get_descendants_and_ranks(cur, child_parent, ranks, row[0])


def get_descendants(cur, node, limit, descendants):
    # Get the children and maybe iterate
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE object = ? AND predicate = 'rdfs:subClassOf'""",
        (node,),
    )
    for row in cur.fetchall():
        if row[0] == limit:
            return
        descendants.append(row[0])
        get_descendants(cur, row[0], limit, descendants)


def move_to_other(cur, parent_tax_id, parent_tax_label, others, precious_others=None):
    cur.execute(
        f"""INSERT INTO statements (stanza, subject, predicate, object, value) VALUES
        ('iedb-taxon:{parent_tax_id}-other',
         'iedb-taxon:{parent_tax_id}-other',
         'rdfs:subClassOf',
         'NCBITaxon:{parent_tax_id}',
         null),
        ('iedb-taxon:{parent_tax_id}-other',
         'iedb-taxon:{parent_tax_id}-other',
         'rdfs:label',
         null,
         'Other {parent_tax_label}');"""
    )
    for o in others:
        # Create map of parent -> child and ranks
        child_parent = {}
        ranks = {}
        get_descendants_and_ranks(cur, child_parent, ranks, o)

        # Find all species-level
        species = [x for x, y in ranks.items() if y == "NCBITaxon:species"]

        # These get bumped up to 'other' then all extra nodes get deleted
        species_str = ", ".join([f"'{x}'" for x in species])
        cur.execute(
            f"""UPDATE statements SET object = 'iedb-taxon:{parent_tax_id}-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({species_str})"""
        )

        # Find nodes to move to 'other organism'
        other_organisms = set()
        for s in species:
            if s not in child_parent:
                continue
            ancestors = get_all_ancestors(
                child_parent, s, "NCBITaxon:" + parent_tax_id
            )
            if not ancestors:
                continue
            other_organisms.add(ancestors[-1])
        o_str = ", ".join([f"'{x}'" for x in other_organisms])
        cur.execute(
            f"""UPDATE statements SET object = 'iedb-taxon:0100026-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({o_str})"""
        )

    # Finally, move the 'precious' others to the correct level
    if precious_others:
        po_str = ", ".join([f"'{x}'" for x in precious_others])
        cur.execute(
            f"""UPDATE statements SET object = 'iedb-taxon:{parent_tax_id}-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({po_str})"""
        )


def move_up(cur, top_level_id, top_level_label, rank, precious=None, extras=None):
    child_parent = {}
    ranks = {}
    top_level = get_curie(top_level_id)
    get_descendants_and_ranks(cur, child_parent, ranks, top_level)

    # Find all nodes of the given rank
    at_rank = [x for x, y in ranks.items() if y == "NCBITaxon:" + rank]
    if not extras:
        extras = []
    at_rank.extend(extras)

    # Sometimes rank-level nodes may be under an extra, make sure these aren't moved
    keep_in_place = set()
    for f in at_rank:
        ancestors = get_all_ancestors(child_parent, f, top_level)
        if not set(ancestors).isdisjoint(set(extras)):
            # Remove if they are a descendant of an "extra"
            keep_in_place.add(f)
    at_rank = list(set(at_rank) - keep_in_place)

    # Bump all nodes of given rank to top-level
    at_rank_str = ", ".join([f"'{x}'" for x in at_rank])
    cur.execute(
        f"""UPDATE statements SET object = '{top_level}'
        WHERE predicate = 'rdfs:subClassOf' AND subject IN ({at_rank_str})"""
    )

    # Find nodes to remove (ancestors to limit) - excluding at_rank under extras
    other_organisms = set()
    precious_others = set()
    for f in at_rank:
        if f not in child_parent:
            continue
        ancestors = get_all_ancestors(child_parent, f, top_level)
        if not ancestors:
            continue
        if not set(ancestors).isdisjoint(set(extras)):
            # Extras may not be of given rank
            # Make sure we don't accidentally remove them
            continue
        move = ancestors[-1]
        # Check for a descendant that is in precious and make sure to move it to 'other'
        descendants = []
        get_descendants(cur, move, f, descendants)
        if not set(precious).isdisjoint(set(descendants)):
            for x in list(set(precious) & set(descendants)):
                precious_others.add(x)
        # Otherwise, move the last of the ancestors to 'other organism'
        other_organisms.add(ancestors[-1])
    o_str = ", ".join([f"'{x}'" for x in other_organisms])
    cur.execute(
        f"""UPDATE statements SET object = 'iedb-taxon:0100026-other'
        WHERE predicate = 'rdfs:subClassOf' AND subject IN ({o_str})"""
    )

    # Find non-rank level nodes under top-level
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE predicate = 'rdfs:subClassOf' AND object = ?""",
        (top_level,),
    )
    non_at_rank = []
    for row in cur.fetchall():
        if row[0] in extras:
            continue
        r = ranks.get(row[0], "")
        if r != "NCBITaxon:" + rank:
            non_at_rank.append(row[0])

    # Move all species-level terms to "Other" then delete ancestors
    if non_at_rank or precious_others:
        move_to_other(cur, top_level_id, top_level_label, non_at_rank, precious_others=precious_others)


def organize_archeobacterium(cur, precious):
    # move species to top-level
    move_up(cur, "2157", "Archeobacterium", "species", precious=precious)


def organize_bacterium(cur, precious):
    # move phylum to top-level
    move_up(cur, "2", "Bacterium", "phylum", precious=precious)


def organize_eukaryote(cur, precious):
    # First move EUKARYOTES to top leve
    e_str = ", ".join([f"'{x}'" for x in EUKARYOTES])
    cur.execute(
        f"""UPDATE statements SET object = 'NCBITaxon:2759'
        WHERE predicate = 'rdfs:subClassOf' AND subject IN ({e_str})"""
    )

    # Then get all the children of eukaryote
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE predicate = 'rdfs:subClassOf' AND object = 'NCBITaxon:2759';"""
    )
    all_eukaryotes = []
    for row in cur.fetchall():
        all_eukaryotes.append(row[0])

    # Move any non-EUKARYOTES to other
    other_eukaryotes = [x for x in all_eukaryotes if x not in EUKARYOTES]
    if other_eukaryotes:
        cur.execute(
            """INSERT INTO statements (stanza, subject, predicate, object, value) VALUES
            ('iedb-taxon:2759-other',
             'iedb-taxon:2759-other',
             'rdfs:subClassOf',
             'NCBITaxon:2759',
             null),
            ('iedb-taxon:2759-other',
             'iedb-taxon:2759-other',
             'rdfs:label',
             null,
             'Other Eukaryote');"""
        )
        oe_string = ", ".join([f"'{x}'" for x in other_eukaryotes])
        cur.execute(
            f"""UPDATE statements SET object = 'iedb-taxon:2759-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({oe_string})"""
        )

    organize_fungus(cur, precious)
    organize_metazoan(cur)
    organize_seed_plant(cur, precious)


def organize_fungus(cur, precious):
    move_up(cur, "4751", "Fungus", "phylum", precious=precious)


def organize_metazoan(cur):
    # Get direct children of metazoan
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE predicate = 'rdfs:subClassOf' AND object = 'NCBITaxon:33208';"""
    )
    all_metazoans = []
    for row in cur.fetchall():
        all_metazoans.append(row[0])

    # Set the METAZOANS to be direct children of Metazoan
    met_string = ", ".join([f"'{x}'" for x in METAZOANS])
    cur.execute(
        f"""UPDATE statements SET object = 'NCBITaxon:33208'
        WHERE predicate = 'rdfs:subClassOf' AND subject IN ({met_string});"""
    )

    # Handle other metazoans
    other_metazoans = [x for x in all_metazoans if x not in METAZOANS]
    if other_metazoans:
        move_to_other(cur, "33208", "Metazoan (multicellular animal)", other_metazoans)
    organize_vertebrate(cur)


def organize_vertebrate(cur):
    pass


def organize_seed_plant(cur, precious):
    extras = [
        "iedb-taxon:10002037",
        "iedb-taxon:10002038",
        "iedb-taxon:10002036",
        "iedb-taxon:10002035",
    ]
    move_up(cur, "58024", "Spermatophyte (seed plant)", "family", extras=extras, precious=precious)


def organize_virus(cur, precious):
    # Move family to top-level, keep DNA & RNA viruses
    move_up(cur, "10239", "Virus", "family", extras=["iedb-taxon:10002039", "iedb-taxon:10002040"], precious=precious)


def main():
    parser = ArgumentParser()
    parser.add_argument("db")
    parser.add_argument("precious")
    parser.add_argument("output")
    args = parser.parse_args()

    precious = []
    with open(args.precious, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            precious.append(get_curie(row[0]))

    copy_database(args.db, args.output)
    with sqlite3.connect(args.output) as conn:
        cur = conn.cursor()
        cur.execute(
            """INSERT INTO statements (stanza, subject, predicate, object, value) VALUES
            ('iedb-taxon:0100026-other',
             'iedb-taxon:0100026-other',
             'rdfs:subClassOf',
             'OBI:0100026',
             null),
            ('iedb-taxon:0100026-other',
             'iedb-taxon:0100026-other',
             'rdfs:label',
             null,
             'Other Organism');"""
        )
        organize_archeobacterium(cur, precious)
        organize_bacterium(cur, precious)
        organize_eukaryote(cur, precious)
        organize_virus(cur, precious)


if __name__ == "__main__":
    main()
