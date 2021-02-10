#!/usr/bin/env python3

import csv
import sqlite3

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


def find_collapse_nodes(collapse_nodes, active_taxa, cuml_counts, parent_children, prev_nodes, node):
    if node not in active_taxa:
        prev_count = cuml_counts[prev_nodes[0]]
        cur_count = cuml_counts[node]
        if prev_count == cur_count:
            # If count is the same, we keep going until it isn't
            prev_nodes.append(node)
            for c in parent_children[node]:
                find_collapse_nodes(collapse_nodes, active_taxa, cuml_counts, parent_children, prev_nodes, c)

    # Otherwise, we want to remove everything except the last node
    if len(prev_nodes) > 1:
        collapse_to = prev_nodes.pop()
        for n in prev_nodes:
            collapse_nodes[n] = collapse_to

    # Then continue with the next level
    for c in parent_children[node]:
        find_collapse_nodes(collapse_nodes, active_taxa, cuml_counts, parent_children, [node], c)


def get_parent(all_removed, child_parents, node):
    p = child_parents[node]
    if p in all_removed:
        return get_parent(all_removed, child_parents, p)
    return p


def prune(cur, active_taxa, cuml_counts, parent_children, child_parents):
    top = parent_children["OBI:0100026"]
    collapse_nodes = {}
    for t in top:
        for t_child in parent_children[t]:
            find_collapse_nodes(collapse_nodes, active_taxa, cuml_counts, parent_children, [t], t_child)

    print(f"Collapsing {len(collapse_nodes)} nodes...")
    cur.execute(
        """INSERT INTO statements (stanza, subject, predicate, object)
        VALUES ('iedb-taxon:subsumes', 'iedb-taxon:subsumes', 'rdf:type', 'owl:AnnotationProperty')""")
    for remove, replace in collapse_nodes.items():
        new_parent = get_parent(collapse_nodes.keys(), child_parents, replace)
        cur.execute("UPDATE statements SET object = ? WHERE object = ?", (replace, remove))
        cur.execute("DELETE FROM statements WHERE stanza = ?", (remove,))
        cur.execute(
            """INSERT INTO statements (stanza, subject, predicate, value)
            VALUES (?, ?, 'iedb-taxon:subsumes', ?)""",
            (replace, replace, remove))
        cur.execute(
            """INSERT INTO statements (stanza, subject, predicate, object)
            VALUES (?, ?, 'rdfs:subClassOf', ?)""",
            (replace, replace, new_parent))
        if remove in cuml_counts:
            del cuml_counts[remove]

    print("Removing extras...")
    cur.execute("DELETE FROM statements WHERE object = stanza")

    return cuml_counts


def get_child_ancestors(child_ancestors, child_parents, child, node):
    p = child_parents.get(node)
    if not p:
        return
    child_ancestors[child].add(p)
    get_child_ancestors(child_ancestors, child_parents, child, p)


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="Database to add counts to")
    parser.add_argument("active_taxa", help="IEDB active taxa")
    parser.add_argument("counts", help="TSV containing ID -> epitope count")
    parser.add_argument("child_parents", help="TSV containing ID -> [aremt")
    parser.add_argument("output", help="Output database")
    args = parser.parse_args()

    active_taxa = []
    with open(args.active_taxa, "r") as f:
        for line in f:
            if not line.strip():
                continue
            active_taxa.append(get_curie(line.strip()))

    child_parents = {}
    with open(args.child_parents, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            parent = row[1]
            if parent == "NCBITaxon:1":
                parent = "OBI:0100026"
            child_parents[row[0]] = parent

    child_ancestors = defaultdict(set)
    for child in child_parents.keys():
        if child not in child_ancestors:
            child_ancestors[child] = set()
        get_child_ancestors(child_ancestors, child_parents, child, child)

    parent_children = defaultdict(set)
    for child, parent in child_parents.items():
        if child not in parent_children:
            parent_children[child] = set()
        if parent not in parent_children:
            parent_children[parent] = set()
        parent_children[parent].add(child)

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
            cuml_counts = prune(cur_new, active_taxa, cuml_counts, parent_children, child_parents)

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
