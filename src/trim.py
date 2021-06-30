#!/usr/bin/env python3

import csv
import logging
import sqlite3

from argparse import ArgumentParser
from helpers import get_curie


def add_iedb_taxa(cur, iedb_taxa):
    with open(iedb_taxa, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tax_id = get_curie(row["Taxon ID"])
            label = row["Label"]
            cur.execute(
                """
                INSERT INTO statements (stanza, subject, predicate, object, value) VALUES
                (?, ?, "rdf:type", "owl:Class", null),
                (?, ?, "rdfs:label", null, ?);
            """,
                (tax_id, tax_id, tax_id, tax_id, label),
            )
            for pid in row["Parent IDs"].split("|"):
                cur.execute(
                    """
                    INSERT INTO statements (stanza, subject, predicate, object) VALUES
                    (?, ?, "rdfs:subClassOf", ?);
                """,
                    (tax_id, tax_id, get_curie(pid)),
                )


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="NCBITaxon database")
    parser.add_argument("active_taxa", help="Active tax IDs used in IEDB")
    parser.add_argument("iedb_taxa", help="IEDB custom taxa")
    parser.add_argument("output", help="Output database")
    args = parser.parse_args()

    weights = {}
    with sqlite3.connect(args.db) as conn:
        cur = conn.cursor()
        cur.execute(
            """SELECT DISTINCT stanza FROM statements
            WHERE stanza LIKE 'iedb-taxon:%' OR stanza LIKE 'NCBITaxon:%'"""
        )
        for row in cur.fetchall():
            tax_id = row[0]
            try:
                int(tax_id.split(":")[1])
            except ValueError:
                # Not a tax ID
                continue
            weights[tax_id] = 0

        active_tax_ids = set()
        with open(args.active_taxa, "r") as f:
            for line in f:
                tax_id = get_curie(line.strip())
                active_tax_ids.add(tax_id)

        # Add the parents of IEDB taxa
        with open(args.iedb_taxa, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                for pid in row["Parent IDs"].split("|"):
                    active_tax_ids.add(get_curie(pid))

        # print(f"Retrieving ancestors for {len(active_tax_ids)} active taxa...")

        for act_tax in active_tax_ids:
            weights[act_tax] = 1
            cur.execute(
                """WITH RECURSIVE active(node) AS (
                VALUES (?)
                UNION
                 SELECT object AS node
                FROM statements
                WHERE predicate = 'rdfs:subClassOf'
                  AND object = ?
                UNION
                SELECT object AS node
                FROM statements, active
                WHERE active.node = statements.stanza
                  AND statements.predicate = 'rdfs:subClassOf'
                  AND statements.object NOT LIKE '_:%'
              )
              SELECT * FROM active""",
                (act_tax, act_tax),
            )
            for row in cur.fetchall():
                parent_id = row[0]
                if parent_id == act_tax:
                    continue
                weights[parent_id] = 1

        active_nodes = [x for x, y in weights.items() if y > 0]
        # print(f"Filtering for {len(active_nodes)} active nodes...")

        # use f-string because we don't know how many values we have
        active_nodes = ", ".join([f"'{x}'" for x in active_nodes])
        cur.execute(f"SELECT * FROM statements WHERE stanza IN ({active_nodes})")
        insert = []
        for r in cur.fetchall():
            vals = []
            for itm in r:
                if not itm:
                    vals.append("null")
                else:
                    itm = itm.replace("'", "''")
                    vals.append(f"'{itm}'")
            insert.append(", ".join(vals))
        insert = ", ".join([f"({x})" for x in insert])

        # print(f"Adding {len(rows)} stanzas to new database...")
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

            # Add the IEDB taxa
            add_iedb_taxa(cur_new, args.iedb_taxa)

            # Check for active taxa not in database
            cur_new.execute("SELECT DISTINCT stanza FROM statements WHERE object = 'owl:Class'")
            existing_ids = set([x[0] for x in cur_new.fetchall()])
            missing = set(active_tax_ids) - existing_ids
            if missing:
                logging.error(
                    f"{len(missing)} active taxa missing from NCBITaxonomy:\n- "
                    + "\n- ".join(missing)
                )


if __name__ == "__main__":
    main()
