import csv
import sqlite3

from argparse import ArgumentParser
from collections import defaultdict
from helpers import (
    copy_database,
    create_other,
    get_child_ancestors,
    get_count_map,
    get_cumulative_counts,
    get_curie,
)


def get_collapse(data, collapse, prev_nodes, threshold=0.99):
    child_parents = data["child_parents"]
    counts = data["counts"]
    precious = data["precious"]

    last_node = prev_nodes[-1]
    current_node = child_parents.get(last_node)
    if not current_node or current_node == "NCBITaxon:1" or current_node.endswith("other"):
        return

    prev_count = counts.get(last_node, 0)
    cur_count = counts.get(current_node, 0)
    if current_node not in precious and prev_count / cur_count > threshold:
        # Continue to check the next one until we hit one with less than threshold
        prev_nodes.append(current_node)
        get_collapse(data, collapse, prev_nodes, threshold=threshold)

    # End with prev node
    if prev_nodes[0] != last_node:
        prev_nodes.append(last_node)
        collapse.append(prev_nodes)

    # Continue with this node as the start node
    prev_nodes = [current_node]
    get_collapse(data, collapse, prev_nodes, threshold=threshold)


def prune(cur, data, threshold=0.99):
    child_parents = data["child_parents"]
    counts = data["counts"]

    # Start from bottom nodes
    cur.execute(
        """SELECT DISTINCT stanza FROM statements
           WHERE object = 'owl:Class'
           AND stanza NOT IN (SELECT object FROM statements
            WHERE predicate = 'rdfs:subClassOf')"""
    )
    for res in cur.fetchall():
        s = res[0]
        if counts.get(s, 0) == 0:
            continue
        # Go up until we find a parent that does not have > 99% of epitopes
        collapse = []
        get_collapse(data, collapse, [s], threshold=threshold)
        for c in collapse:
            first_node = c[0]
            last_node = c[-1]
            parent_node = child_parents.get(last_node)
            if not parent_node:
                continue
            # Create the "other" class for the parent
            cur.execute(
                "SELECT value FROM statements WHERE stanza = ? AND predicate = 'rdfs:label'",
                (parent_node,),
            )
            res = cur.fetchone()
            if res:
                parent_label = res[0]
            else:
                parent_label = parent_node
            parent_tax = parent_node.split(":")[1]
            create_other(cur, parent_tax, parent_label)

            cur.execute(
                "UPDATE statements SET object = ? WHERE stanza = ? AND predicate = 'rdfs:subClassOf'",
                (parent_node, first_node),
            )
            cur.execute(
                "UPDATE statements SET object = ? WHERE stanza = ? AND predicate = 'rdfs:subClassOf'",
                (f"iedb-taxon:{parent_tax}-other", last_node),
            )


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

    count_map = get_count_map(args.counts)

    data = {
        "child_ancestors": child_ancestors,
        "child_parents": child_parents,
        "precious": precious,
    }

    copy_database(args.db, args.output)
    with sqlite3.connect(args.output) as conn:
        cur = conn.cursor()
        cuml_counts = get_cumulative_counts(count_map, child_ancestors)
        data["counts"] = cuml_counts
        prune(cur, data)


if __name__ == "__main__":
    main()
