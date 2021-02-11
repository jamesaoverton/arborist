#!/usr/bin/env python3

import csv
import sqlite3
import sys

from argparse import ArgumentParser
from collections import defaultdict

def build_counts(tax_id, child_parent, add_counts, count):
    if tax_id in add_counts:
        cur_count = add_counts[tax_id]
        count = cur_count + count
        add_counts[tax_id] = count
    else:
        add_counts[tax_id] = count
    parent = child_parent.get(tax_id, "OBI:0100026")
    if parent == "OBI:0100026":
        return add_counts
    return build_counts(parent, child_parent, add_counts, count)


def get_curie(tax_id):
    if tax_id.startswith("OBI:"):
        return tax_id
    if len(tax_id) == 8 and tax_id.startswith("100"):
        return "iedb-taxon:" + tax_id
    return "NCBITaxon:" + tax_id


def get_cumulative_counts(cur, count_map, child_ancestors):
    cuml_counts = {}
    for child, ancestors in child_ancestors.items():
        count = count_map.get(child, 0)
        if child in cuml_counts:
            cur_count = cuml_counts[child]
            cuml_counts[child] = cur_count + count
        else:
            cuml_counts[child] = count

        for a in ancestors:
            if a in cuml_counts:
                cur_count = cuml_counts[a]
                cuml_counts[a] = cur_count + count
            else:
                cuml_counts[a] = count
    return cuml_counts


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
        find_collapse_nodes(collapse_nodes, precious, cuml_counts, child_parents, prev_nodes, parent)

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
    if node in precious and node != "OBI:0100026":
        return find_start(precious, child_parents, child_parents[node])
    return node


def clean_collapse_nodes(collapse_nodes, node):
    replace = collapse_nodes[node]
    if replace in collapse_nodes:
        return clean_collapse_nodes(collapse_nodes, replace)
    return replace


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
        parent = child_parents[s2]
        # From there, move up the tree to find nodes to collapse
        find_collapse_nodes(collapse_nodes, precious, cuml_counts, child_parents, [s2], parent)

    clean = {}
    for remove in collapse_nodes.keys():
        replace = clean_collapse_nodes(collapse_nodes, remove)
        clean[remove] = replace

    print(f"Collapsing {len(clean)} nodes...")
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

    print("Removing circular logic...")
    cur.execute("DELETE FROM statements WHERE object = stanza")

    return cuml_counts


def get_child_ancestors(child_ancestors, child_parents, child, node):
    p = child_parents.get(node)
    if not p or p == node:
        return
    child_ancestors[child].add(p)
    get_child_ancestors(child_ancestors, child_parents, child, p)


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="Database to add counts to")
    parser.add_argument("precious", help="List of taxa to keep")
    parser.add_argument("counts", help="TSV containing ID -> epitope count")
    parser.add_argument("child_parents", help="TSV containing ID -> [aremt")
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

    with sqlite3.connect(args.db) as conn:
        # Get stanzas from source database
        cur = conn.cursor()
        cur.execute("SELECT * FROM statements")
        rows = cur.fetchall()
        insert = []
        for r in rows:
            vals = []
            for itm in r:
                if not itm:
                    vals.append("null")
                else:
                    itm = itm.replace("'", "''")
                    vals.append(f"'{itm}'")
            insert.append(", ".join(vals))

        # print(f"Inserting {len(insert)} statements into new database...")
        insert = ", ".join([f"({x})" for x in insert])

        # Insert all into target database then add counts
        with sqlite3.connect(args.output) as conn_new:
            cur_new = conn_new.cursor()
            cur_new.execute(
                """CREATE TABLE statements (stanza TEXT,
                                                        subject TEXT,
                                                        predicate TEXT,
                                                        object TEXT,
                                                        value TEXT,
                                                        datatype TEXT,
                                                        language TEXT)"""
            )
            cur_new.execute("INSERT INTO statements VALUES " + insert)
            cur_new.execute("CREATE INDEX stanza_idx ON statements (stanza)")
            cur_new.execute("CREATE INDEX subject_idx ON statements (subject)")
            cur_new.execute("CREATE INDEX object_idx ON statements (object)")
            cur_new.execute("ANALYZE")
            count_map = {}
            with open(args.counts, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if row[0] == "NULL":
                        continue
                    count_map[get_curie(row[0])] = int(row[1])

            cuml_counts = get_cumulative_counts(cur_new, count_map, child_ancestors)
            cuml_counts = prune(cur_new, precious, cuml_counts, child_parents)

            print("Adding epitope counts...")
            for tax_id, count in cuml_counts.items():
                cur_new.execute(
                    f"""UPDATE statements SET value = value || ' ({count})'
                    WHERE stanza = '{tax_id}'
                      AND subject = '{tax_id}'
                      AND predicate = 'rdfs:label';"""
                )


if __name__ == "__main__":
    main()
