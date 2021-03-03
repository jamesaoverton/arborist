#!/usr/bin/env python3

import csv
import sqlite3
import sys

from argparse import ArgumentParser
from collections import defaultdict
from helpers import copy_database, get_curie


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


def move_to_other(
    cur, parent_tax_id, parent_tax_label, others, rank="species", precious_others=None
):
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

    if rank == "none":
        # Just set the others to be children of "Other" and return
        others_str = ", ".join([f"'{x}'" for x in others])
        cur.execute(
            f"""UPDATE statements SET object = 'iedb-taxon:{parent_tax_id}-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({others_str})"""
        )
        return

    for o in others:
        # Create map of parent -> child and ranks
        child_parent = {}
        ranks = {}
        get_descendants_and_ranks(cur, child_parent, ranks, o)

        # Find all at rank level
        at_rank = [x for x, y in ranks.items() if y == "NCBITaxon:" + rank]

        # These get bumped up to 'other' then all extra nodes get deleted
        at_rank_str = ", ".join([f"'{x}'" for x in at_rank])
        cur.execute(
            f"""UPDATE statements SET object = 'iedb-taxon:{parent_tax_id}-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({at_rank_str})"""
        )

        # Find nodes to move to 'other organism'
        other_organisms = set()
        for s in at_rank:
            if s not in child_parent:
                continue
            ancestors = get_all_ancestors(child_parent, s, "NCBITaxon:" + parent_tax_id)
            if not ancestors:
                continue
            other_organisms.add(ancestors[-1])
        if not other_organisms:
            other_organisms.add(o)
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
        move_to_other(
            cur, top_level_id, top_level_label, non_at_rank, precious_others=precious_others
        )


def organize(cur, top_level, precious):
    for curie, details in top_level.items():
        parent = get_curie(details["Parent ID"])

        # First, rehome this node
        cur.execute(
            "UPDATE statements SET object = ? WHERE subject = ? AND predicate = 'rdfs:subClassOf'",
            (parent, curie),
        )

        rank = details.get("Rank", "").strip()
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
            other_rank = details.get("Others", "")
            if other_rank.strip() == "":
                other_rank = "species"
            move_to_other(cur, details["ID"], details["Label"], others, rank=other_rank)
            continue

        # Otherwise, move all of given rank to top-level
        extras = [
            get_curie(x) for x in details.get("Extra Nodes", "").split(", ") if x.strip() != ""
        ]
        move_up(cur, details["ID"], details["Label"], rank, extras=extras, precious=precious)


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

    top_level = {}
    with open(args.top_level, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            top_level[get_curie(row["ID"])] = row

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
