#!/usr/bin/env python3

import sqlite3

from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="NCBITaxon database")
    parser.add_argument("active_taxa", help="Active tax IDs used in IEDB")
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
                tax_int = int(tax_id.split(":")[1])
            except ValueError:
                # Not a tax ID
                continue
            weights[tax_id] = 0

        active_tax_ids = []
        with open(args.active_taxa, "r") as f:
            for line in f:
                tax_id = line.strip()
                if len(tax_id) == 8 and tax_id.startswith("100"):
                    tax_id = "iedb-taxon:" + tax_id
                elif not tax_id.startswith("OBI:"):
                    tax_id = "NCBITaxon:" + tax_id
                active_tax_ids.append(tax_id)

        # print(f"Retrieving ancestors for {len(active_tax_ids)} active taxa...")

        for act_tax in active_tax_ids:
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
                weights[parent_id] = 1

        active_nodes = [x for x, y in weights.items() if y > 0]
        # print(f"Filtering for {len(active_nodes)} active nodes...")

        # use f-string because we don't know how many values we have
        active_nodes = ", ".join([f"'{x}'" for x in active_nodes])
        cur.execute(f"SELECT * FROM statements WHERE stanza IN ({active_nodes})")
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


if __name__ == "__main__":
    main()
