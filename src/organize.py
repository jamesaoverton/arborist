#!/usr/bin/env python3

import csv
import sqlite3

from argparse import ArgumentParser
from collections import defaultdict
from helpers import (
    copy_database,
    get_all_ancestors,
    get_curie,
    get_descendants,
    get_descendants_and_ranks,
    move_precious_to_other,
    move_rank_to_other,
)


def get_top_ancestor(child_parent, node, limit):
    parent = child_parent[node]
    if parent == limit:
        return node
    return get_top_ancestor(child_parent, parent, limit)


def move_up(cur, top_level_id, top_level_label, rank, precious=None, extras=None):
    # Init empty lists if these are not included
    if not precious:
        precious = []
    if not extras:
        extras = []

    # Get all descendants and their ranks under this top level
    child_parent = {}
    ranks = {}
    top_level = get_curie(top_level_id)
    get_descendants_and_ranks(cur, child_parent, ranks, top_level)

    # Find all nodes of the given rank
    at_rank = [x for x, y in ranks.items() if y == "NCBITaxon:" + rank]
    if not extras:
        extras = []
    at_rank.extend(extras)

    # Sometimes rank-level nodes may be under an extra/precious, make sure these aren't moved
    keep_in_place = set()
    # If a rank-level term has a precious ancestor,
    # we want to keep that ancestor as the direct child instead
    replace_as_child = {}
    for taxa in at_rank:
        ancestors = reversed(get_all_ancestors(child_parent, taxa, top_level))
        if not set(ancestors).isdisjoint(set(extras)) or set(ancestors).isdisjoint(set(precious)):
            # Remove if they are a descendant of an extra or precious node
            keep_in_place.add(taxa)
            # Find which to move directly under top level
            precious_ancestors = list(set(ancestors).intersection(set(precious)))
            if precious_ancestors:
                replace_as_child[taxa] = precious_ancestors[0]
    at_rank = set(at_rank) - keep_in_place
    for taxa, replacement in replace_as_child.items():
        if taxa in at_rank:
            at_rank.remove(taxa)
        at_rank.add(replacement)

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
        if not set(ancestors).isdisjoint(set(extras)) or not set(ancestors).isdisjoint(
                set(precious)):
            # Extras & precious may not be of given rank
            # Make sure we don't accidentally remove them
            continue
        move = ancestors[-1]
        # Check for a descendant that is in precious and make sure to move it to 'other'
        descendants = []
        get_descendants(cur, move, [f], descendants)
        if not set(precious).isdisjoint(set(descendants)):
            for x in list(set(precious) & set(descendants)):
                precious_others.add(x)
        # Otherwise, move the last of the ancestors to 'other organism'
        other_organisms.add(move)

    # Find non-rank level nodes under top-level
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE predicate = 'rdfs:subClassOf' AND object = ?""",
        (top_level,),
    )
    non_at_rank = []
    for row in cur.fetchall():
        tax_id = row[0]
        if tax_id in extras or tax_id in precious:
            continue
        r = ranks.get(tax_id, "")
        if r != "NCBITaxon:" + rank:
            non_at_rank.append(tax_id)

    # Move all species-level or precious terms to "Other" then delete ancestors
    if non_at_rank or precious_others:
        move_precious_to_other(cur, precious, top_level_id, top_level_label, non_at_rank)
        # move_rank_to_other(cur, top_level_id, top_level_label, non_at_rank, precious=precious)

    o_str = ", ".join([f"'{x}'" for x in other_organisms])
    cur.execute(
        f"""UPDATE statements SET object = 'iedb-taxon:0100026-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({o_str})"""
    )


def organize(cur, top_level, precious):
    for curie, details in top_level.items():
        parent = get_curie(details["Parent ID"])

        # First, rehome this node
        cur.execute(
            "UPDATE statements SET object = ? WHERE subject = ? AND predicate = 'rdfs:subClassOf'",
            (parent, curie),
        )

        rank = details.get("Child Rank", "").strip()
        if rank == "":
            # Nothing to do, everything stays as-is
            continue

        if rank == "manual":
            # Everything NOT in this set gets moved to other
            cur.execute(
                """SELECT DISTINCT subject FROM statements
                WHERE object = ? AND predicate = 'rdfs:subClassOf'""",
                (curie,),
            )
            others = []
            for row in cur.fetchall():
                if row[0] not in top_level:
                    others.append(row[0])
            other_rank = details.get("Other Rank", "")
            if other_rank.strip() == "":
                other_rank = "species"
            move_rank_to_other(
                cur, details["ID"], details["Label"], others, rank=other_rank, precious=precious
            )
            continue

        # Otherwise, move all of given rank to top-level
        extras = [
            get_curie(x) for x in details.get("Extra Nodes", "").split(", ") if x.strip() != ""
        ]
        move_up(cur, details["ID"], details["Label"], rank, extras=extras, precious=precious)


def get_line(top_structure, line, node):
    children = top_structure.get(node)
    if not children:
        return line
    line.extend(children)
    for c in children:
        get_line(top_structure, line, c)


def main():
    parser = ArgumentParser()
    parser.add_argument("db")
    parser.add_argument("top_level")
    parser.add_argument("precious")
    parser.add_argument("output")
    args = parser.parse_args()

    precious = []
    with open(args.precious, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            precious.append(get_curie(row[0]))

    # ID -> Details
    top_level_unordered = {}
    # Parent -> Children
    top_structure = defaultdict(set)
    with open(args.top_level, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tax_id = get_curie(row["ID"])
            parent_id = get_curie(row["Parent ID"])
            if parent_id not in top_structure:
                top_structure[parent_id] = set()
            top_structure[parent_id].add(tax_id)
            top_level_unordered[tax_id] = row

    # Sort top level by structure (starting with children of cellular organism)
    orgs = top_structure["OBI:0100026"]
    full_line = []
    for o in orgs:
        line = [o]
        get_line(top_structure, line, o)
        full_line.extend(line)
    # Go from lowest -> highest level
    full_line.reverse()
    top_level = {node: top_level_unordered[node] for node in full_line}

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
        organize(cur, top_level, precious)


if __name__ == "__main__":
    main()
