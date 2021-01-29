#!/usr/bin/env python3

import sqlite3

from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="NCBITaxon database")
    parser.add_argument("active_taxa", help="Active tax IDs used in IEDB")
    args = parser.parse_args()

    weights = {}
    with sqlite3.connect(args.db) as conn:
        cur = conn.cursor()
        cur.execute("SELECT DISTINCT stanza FROM statements WHERE stanza LIKE 'iedb-taxon:%' OR stanza LIKE 'NCBITaxon:%'")
        for row in cur.fetchall():
            tax_id = row[0]
            try:
                tax_int = int(tax_id.split(":")[1])
            except ValueError:
                # Not a tax ID
                continue
            weights[tax_id] = 0

        with open(args.active_taxa, "r") as f:
            for line in f:
                tax_id = line.strip()
                if len(tax_id) == 8 and tax_id.startswith("100"):
                    tax_id = "iedb-taxon:" + tax_id
                elif not tax_id.startswith("OBI:"):
                    tax_id = "NCBITaxon:" + tax_id
                cur.execute(
                f"""WITH RECURSIVE ancestors(parent) AS (
                    VALUES ('{tax_id}')
                    UNION
                     SELECT object AS parent
                    FROM statements
                    WHERE predicate = 'rdfs:subClassOf'
                      AND object = '{tax_id}'
                    UNION
                    SELECT object AS parent
                    FROM statements, ancestors
                    WHERE ancestors.parent = statements.stanza
                      AND statements.predicate = 'rdfs:subClassOf'
                      AND statements.object NOT LIKE '_:%'
                  )
                  SELECT * FROM ancestors"""
                )
                for row in cur.fetchall():
                    parent_id = row[0]
                    weights[parent_id] = 1

        remove = [x for x, y in weights.items() if y == 0]
        remove = ", ".join([f"'{x}'" for x in remove])
        print(f"Removing {len(remove)} inactive nodes...")
        cur.execute(f"DELETE FROM statements WHERE stanza IN ({remove})")
        cur.execute(f"DELETE FROM statements WHERE object IN ({remove})")


if __name__ == '__main__':
    main()


