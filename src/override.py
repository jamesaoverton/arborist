#!/usr/bin/env python3

import csv
import sqlite3

from argparse import ArgumentParser
from collections import defaultdict
from helpers import clean_no_epitopes, get_child_ancestors, get_child_parents, get_count_map, get_cumulative_counts, get_curie


def update_names(cur, names):
    # Replace cellular organisms with OBI organism
    cur.execute(
        "UPDATE statements SET stanza = 'OBI:0100026' WHERE stanza = 'NCBITaxon:131567';"
    )
    cur.execute(
        "UPDATE statements SET subject = 'OBI:0100026' WHERE subject = 'NCBITaxon:131567';"
    )
    cur.execute(
        """UPDATE statements SET value = 'Organism'
        WHERE stanza = 'OBI:0100026'
          AND subject = 'OBI:0100026'
          AND predicate = 'rdfs:label';"""
    )
    cur.execute(
        "UPDATE statements SET object = 'OBI:0100026' WHERE object = 'NCBITaxon:131567';"
    )
    bnode_id = 1
    with open(names, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tax_id = row["Taxon ID"]
            new_label = row["Label"]
            label_source = row["Label Source"]
            cur.execute(
                """UPDATE statements SET value = ?
                WHERE stanza = ? AND subject = ? AND predicate = 'rdfs:label'""",
                (new_label, tax_id, tax_id),
            )
            cur.execute(
                """INSERT INTO statements (stanza, subject, predicate, value) VALUES
                (?, ?, 'oio:hasLabelSource', ?)""", (tax_id, f"_:bnode{bnode_id:08}", label_source))
            bnode_id += 1
            for syn in row["Synonyms"].split(","):
                if syn.strip() == "":
                    continue
                cur.execute(
                    """INSERT INTO statements (stanza, subject, predicate, value) VALUES
                    (?, ?, "oio:hasExactSynonym", ?)""",
                    (tax_id, tax_id, syn),
                )


def update_parents(cur, parents):
    with open(parents, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tax_id = get_curie(row["Taxon ID"])
            parent_id = get_curie(row["Parent ID"])
            cur.execute(
                """
                UPDATE statements
                SET object = ?
                WHERE stanza = ? AND subject = ? AND predicate = 'rdfs:subClassOf';
            """,
                (parent_id, tax_id, tax_id),
            )


def update(source, target, precious, counts, names, parents):
    with sqlite3.connect(source) as conn:
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
        insert = ", ".join([f"({x})" for x in insert])

        # Insert all into target database then run updates
        with sqlite3.connect(target) as conn_new:
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
            cur_new.execute(f"INSERT INTO statements VALUES " + insert)

            # Override labels with IEDB labels
            print("updating labels...")
            update_names(cur_new, names)

            # Override parents
            print("updating parents...")
            update_parents(cur_new, parents)

            # Get updated child->ancestors
            child_parents = get_child_parents(cur_new)
            child_ancestors = defaultdict(set)
            for child in child_parents.keys():
                if child not in child_ancestors:
                    child_ancestors[child] = set()
                get_child_ancestors(child_ancestors, child_parents, child, child)

            # Use child->ancestors to get updated cumulative epitope counts
            count_map = get_count_map(counts)
            cuml_counts = get_cumulative_counts(count_map, child_ancestors)

            precious_terms = []
            with open(precious, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    precious_terms.append(get_curie(row[0]))

            # Clean up zero-epitope terms
            clean_no_epitopes(cur_new, cuml_counts, precious_terms)


def main():
    parser = ArgumentParser()
    parser.add_argument("db", type=str, help="Existing NCBITaxon database to update")
    parser.add_argument("name_overrides", type=str, help="Label overrides")
    parser.add_argument("parent_overrides", type=str, help="Parent taxa overrides")
    parser.add_argument("precious")
    parser.add_argument("child_parents", type=str, help="Child parent map from database to update")
    parser.add_argument("counts", type=str, help="Epitope counts")
    parser.add_argument("output", type=str, help="Output database")
    args = parser.parse_args()
    update(args.db, args.output, args.precious, args.counts, args.name_overrides, args.parent_overrides)


if __name__ == "__main__":
    main()
