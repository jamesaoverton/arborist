#!/usr/bin/env python3

import csv
import sqlite3
import sys

from argparse import ArgumentParser
from collections import defaultdict
from helpers import copy_database


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


def get_all_ancestors(child_parent, ancestors, node):
    ancestors.append(node)
    if node in child_parent:
        get_all_ancestors(child_parent, ancestors, child_parent[node])


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


def move_to_other(cur, parent_tax_id, parent_tax_label, others):
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

        # Find nodes to delete
        for s in species:
            if s not in child_parent:
                continue
            ancestors = []
            get_all_ancestors(child_parent, ancestors, child_parent[s])
            ancestors_str = ", ".join([f"'{x}'" for x in ancestors])
            cur.execute(f"DELETE FROM statements WHERE subject IN ({ancestors_str})")


def move_up(cur, top_level_id, top_level_label, rank, extras=None):
    child_parent = {}
    ranks = {}
    top_level = "NCBITaxon:" + top_level_id
    get_descendants_and_ranks(cur, child_parent, ranks, top_level)

    # Find all nodes of the given rank
    families = [x for x, y in ranks.items() if y == "NCBITaxon:" + rank]
    if not extras:
        extras = []
    families.extend(extras)

    # Sometimes rank-level nodes may be under an extra, make sure these aren't moved
    keep_in_place = set()
    for f in families:
        ancestors = []
        get_all_ancestors(child_parent, ancestors, child_parent[f])
        if not set(ancestors).isdisjoint(set(extras)):
            # Remove if they are a descendant of an "extra"
            keep_in_place.add(f)
    families = list(set(families) - keep_in_place)

    # Bump all nodes of given rank to top-level
    families_str = ", ".join([f"'{x}'" for x in families])
    cur.execute(
        f"""UPDATE statements SET object = 'NCBITaxon:{top_level_id}'
        WHERE predicate = 'rdfs:subClassOf' AND subject IN ({families_str})"""
    )

    # Find nodes to delete - excluding families under extras
    for f in families:
        if f not in child_parent:
            continue
        ancestors = []
        get_all_ancestors(child_parent, ancestors, child_parent[f])
        if not set(ancestors).isdisjoint(set(extras)):
            # Extras may not be of given rank
            # Make sure we don't accidentally remove them
            continue
        if top_level in ancestors:
            ancestors.remove(top_level)
        ancestors_str = ", ".join([f"'{x}'" for x in ancestors])
        cur.execute(f"DELETE FROM statements WHERE subject IN ({ancestors_str})")

    # Find non-rank level nodes under top-level
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE predicate = 'rdfs:subClassOf' AND object = ?""",
        (top_level,),
    )
    non_families = []
    for row in cur.fetchall():
        if row[0] in extras:
            continue
        r = ranks.get(row[0], "")
        if r != "NCBITaxon:" + rank:
            non_families.append(row[0])

    # Move all species-level terms to "Other" then delete ancestors
    if non_families:
        move_to_other(cur, top_level_id, top_level_label, non_families)


def organize_archeobacterium(cur):
    # Create map of parent -> child and ranks
    child_parent = {}
    ranks = {}
    get_descendants_and_ranks(cur, child_parent, ranks, "NCBITaxon:2157")

    # Find all species-level
    species = [x for x, y in ranks.items() if y == "NCBITaxon:species"]

    # These get bumped up to 'other' then all extra nodes get deleted
    species_str = ", ".join([f"'{x}'" for x in species])
    cur.execute(
        f"""UPDATE statements SET object = 'NCBITaxon:2157'
        WHERE predicate = 'rdfs:subClassOf' AND subject IN ({species_str})"""
    )

    # Find nodes to delete
    for s in species:
        ancestors = []
        get_all_ancestors(child_parent, ancestors, child_parent[s])
        if "NCBITaxon:2157" in ancestors:
            ancestors.remove("NCBITaxon:2157")
        ancestors_str = ", ".join([f"'{x}'" for x in ancestors])
        cur.execute(f"DELETE FROM statements WHERE subject IN ({ancestors_str})")


def organize_bacterium(cur):
    # move phylum to top-level
    move_up(cur, "2", "Bacterium", "phylum")


def organize_eukaryote(cur):
    # First get all the children of eukaryote
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE predicate = 'rdfs:subClassOf' AND object = 'NCBITaxon:2759';"""
    )
    all_eukaryotes = []
    for row in cur.fetchall():
        all_eukaryotes.append(row[0])

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

    organize_fungus(cur)
    organize_metazoan(cur)
    organize_seed_plant(cur)


def organize_fungus(cur):
    move_up(cur, "4751", "Fungus", "phylum")


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


def organize_seed_plant(cur):
    extras = [
        "iedb-taxon:10002037",
        "iedb-taxon:10002038",
        "iedb-taxon:10002036",
        "iedb-taxon:10002035",
    ]
    move_up(cur, "58024", "Spermatophyte (seed plant)", "family", extras=extras)


def organize_virus(cur):
    # Move family to top-level, keep DNA & RNA viruses
    move_up(cur, "10239", "Virus", "family", extras=["iedb-taxon:10002039", "iedb-taxon:10002040"])


def main():
    parser = ArgumentParser()
    parser.add_argument("db")
    parser.add_argument("output")
    args = parser.parse_args()

    copy_database(args.db, args.output)
    with sqlite3.connect(args.output) as conn:
        cur = conn.cursor()
        organize_archeobacterium(cur)
        organize_bacterium(cur)
        organize_eukaryote(cur)
        organize_virus(cur)


if __name__ == "__main__":
    main()
