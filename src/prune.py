#!/usr/bin/env python3

import csv
import sqlite3
import sys

from argparse import ArgumentParser
from collections import defaultdict
from helpers import (
    copy_database,
    get_child_ancestors,
    get_cumulative_counts,
    get_curie,
    get_descendants,
    get_descendants_and_ranks,
)


def get_parent(all_removed, child_parents, node):
    p = child_parents[node]
    if p in all_removed:
        return get_parent(all_removed, child_parents, p)
    return p


def find_collapse_nodes(collapse_nodes, precious, cuml_counts, child_parents, prev_nodes, node):
    if node in collapse_nodes:
        # The rest of this line has already been checked
        if len(prev_nodes) > 1:
            # Collapse previous nodes to bottom
            collapse_to = prev_nodes.pop(0)
            for pn in prev_nodes:
                collapse_nodes[pn] = collapse_to
        return collapse_nodes

    # Compare current and previous counts
    prev_count = cuml_counts[prev_nodes[-1]]
    cur_count = cuml_counts[node]
    if prev_count == cur_count and node not in precious:
        # Continue checking next level parent
        parent = child_parents.get(node)
        if not parent or parent == "OBI:0100026":
            return
        prev_nodes.append(node)
        find_collapse_nodes(
            collapse_nodes, precious, cuml_counts, child_parents, prev_nodes, parent
        )

    # Start again
    if len(prev_nodes) > 1:
        # Collapse previous nodes to bottom
        collapse_to = prev_nodes.pop(0)
        for pn in prev_nodes:
            collapse_nodes[pn] = collapse_to

    # If length is only 1 do not collapse
    parent = child_parents.get(node)
    if not parent or parent == "OBI:0100026":
        return

    find_collapse_nodes(collapse_nodes, precious, cuml_counts, child_parents, [node], parent)


def find_start(precious, child_parents, node):
    if node in precious and node != "NCBITaxon:1":
        return find_start(precious, child_parents, child_parents[node])
    return node


def clean_collapse_nodes(collapse_nodes, node):
    replace = collapse_nodes[node]
    if replace in collapse_nodes:
        return clean_collapse_nodes(collapse_nodes, replace)
    return replace


def move_species_up(cur, precious, nodes, limit=20):
    for tax_id in nodes:
        child_parent = {}
        ranks = {}
        get_descendants_and_ranks(cur, child_parent, ranks, tax_id)

        # Find all at species-level
        species = [x for x, y in ranks.items() if y == "NCBITaxon:species"]

        # Add any precious terms
        for child in child_parent.keys():
            if child in precious:
                species.append(child)

        if len(species) < limit:
            # Find list of all descendants, and if there is a precious node in descendants, use that
            descendants = []
            get_descendants(cur, tax_id, species, descendants)
            precious_descendants = set(descendants).intersection(set(precious))

            replace = None
            if precious_descendants:
                # if len(precious_descendants) > 1:
                    # print("More than one precious under " + tax_id)
                    # print(precious_descendants)

                # Move this term to replace the current tax_id
                replace = list(precious_descendants)[0]
                cur.execute(
                    """SELECT object FROM statements
                    WHERE subject = ? AND predicate = 'rdfs:subClassOf'""",
                    (tax_id,),
                )

                parent = cur.fetchone()[0]
                cur.execute(
                    """UPDATE statements SET object = ?
                    WHERE subject = ? AND predicate = 'rdfs:subClassOf'""",
                    (parent, replace),
                )

            # Get direct children of this node
            cur.execute(
                """SELECT DISTINCT subject FROM statements
                WHERE predicate = 'rdfs:subClassOf' AND object = ?""",
                (tax_id,),
            )

            # Move the non-species to other organism
            move_to_other = [x[0] for x in cur.fetchall() if x[0] != replace]
            move_str = ", ".join([f"'{x}'" for x in move_to_other])
            cur.execute(
                f"""UPDATE statements SET object = 'iedb-taxon:0100026-other'
                WHERE subject IN ({move_str}) AND predicate = 'rdfs:subClassOf'"""
            )
            if replace:
                # Get rid of the old tax_id
                cur.execute(
                    """UPDATE statements SET object = 'iedb-taxon:0100026-other'
                    WHERE subject = ? AND predicate = 'rdfs:subClassOf'""",
                    (tax_id,),
                )

            # Move the species back to the node
            species_str = ", ".join([f"'{x}'" for x in species])
            if replace:
                tax_id = replace
            cur.execute(
                f"""UPDATE statements SET object = '{tax_id}'
                WHERE subject IN ({species_str}) AND predicate = 'rdfs:subClassOf'"""
            )
            continue
        """else:
            # Work down until we find the node with < limit species
            cur.execute(
                ""SELECT DISTINCT subject FROM statements
                WHERE predicate = 'rdfs:subClass' AND object = ?"",
                (tax_id,)
            )
            children = [x[0] for x in cur.fetchall()]
            move_species_up(cur, precious, children)"""


def prune(cur, precious, cuml_counts, child_parents):
    # Start from bottom level nodes
    cur.execute(
        """SELECT DISTINCT stanza FROM statements
           WHERE object = 'owl:Class'
           AND stanza NOT IN (SELECT object FROM statements
            WHERE predicate = 'rdfs:subClassOf')"""
    )
    start = []
    for row in cur.fetchall():
        start.append(row[0])

    collapse_nodes = {}
    i = 0
    for s in start:
        # Make sure we start with a non-precious node
        s2 = find_start(precious, child_parents, s)
        if s2 == "NCBITaxon:1":
            continue
        parent = child_parents[s2]
        # From there, move up the tree to find nodes to collapse
        this_collapse = {}
        find_collapse_nodes(this_collapse, precious, cuml_counts, child_parents, [s2], parent)
        collapse_nodes.update(this_collapse)

    clean = {}
    for remove in collapse_nodes.keys():
        replace = clean_collapse_nodes(collapse_nodes, remove)
        clean[remove] = replace

    # print(f"Collapsing {len(clean)} nodes...")
    for remove, replace in clean.items():
        new_parent = get_parent(collapse_nodes.keys(), child_parents, replace)
        cur.execute("UPDATE statements SET object = ? WHERE object = ?", (replace, remove))
        cur.execute("DELETE FROM statements WHERE stanza = ?", (remove,))
        cur.execute(
            """INSERT INTO statements (stanza, subject, predicate, object)
            VALUES (?, ?, 'rdfs:subClassOf', ?)""",
            (replace, replace, new_parent),
        )
        if remove in cuml_counts:
            del cuml_counts[remove]

    # print("Removing circular logic...")
    cur.execute("DELETE FROM statements WHERE object = stanza")

    print("two")
    cur.execute("SELECT object FROM statements WHERE predicate = 'rdfs:subClassOf' AND stanza = 'NCBITaxon:5654'")
    print(cur.fetchone())

    ranks = ["superorder", "order", "suborder", "superfamily", "family", "subfamily"]
    for r in ranks:
        # print(f"Moving species up for {r}...")
        cur.execute(f"SELECT DISTINCT subject FROM statements WHERE object = 'NCBITaxon:{r}'")
        at_rank = [x[0] for x in cur.fetchall()]
        move_species_up(cur, precious, at_rank)

    print("three")
    cur.execute("SELECT object FROM statements WHERE predicate = 'rdfs:subClassOf' AND stanza = 'NCBITaxon:5654'")
    print(cur.fetchone())


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="Database to add counts to")
    parser.add_argument("precious", help="List of taxa to keep")
    parser.add_argument("counts", help="TSV containing ID -> epitope count")
    parser.add_argument("child_parents", help="TSV containing ID -> parent")
    parser.add_argument("output", help="Output database")
    args = parser.parse_args()

    precious = []
    with open(args.precious, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            precious.append(get_curie(row[0]))

    child_parents = {}
    with open(args.child_parents, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row[0] == row[1]:
                continue
            parent = row[1]
            child_parents[row[0]] = parent

    child_ancestors = defaultdict(set)
    for child in child_parents.keys():
        if child not in child_ancestors:
            child_ancestors[child] = set()
        get_child_ancestors(child_ancestors, child_parents, child, child)

    count_map = {}
    with open(args.counts, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row[0] == "NULL":
                continue
            count_map[get_curie(row[0])] = int(row[1])

    copy_database(args.db, args.output)
    with sqlite3.connect(args.output) as conn:
        cur = conn.cursor()

        cuml_counts = get_cumulative_counts(cur, count_map, child_ancestors)
        prune(cur, precious, cuml_counts, child_parents)


if __name__ == "__main__":
    main()
