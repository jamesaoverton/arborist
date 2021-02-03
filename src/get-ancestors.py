import csv
import sqlite3

from argparse import ArgumentParser


def get_curie(tax_id):
    if tax_id.startswith("OBI:"):
        return tax_id
    if len(tax_id) == 8 and tax_id.startswith("100"):
        return "iedb-taxon:" + tax_id
    return "NCBITaxon:" + tax_id


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="Database to get nodes from")
    parser.add_argument("counts", help="TSV containing epitope counts")
    parser.add_argument("output", help="Output TSV")
    args = parser.parse_args()

    with sqlite3.connect(args.db) as conn:
        cur = conn.cursor()
        ids = []
        with open(args.counts, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                ids.append(get_curie(row[0]))

        child_ancestors = {}
        print(f"Getting ancestors for {len(ids)} taxa...")
        for tax_id in ids:
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
                (tax_id, tax_id),
            )
            ancestors = []
            for row in cur.fetchall():
                ancestors.append(row[0])
            child_ancestors[tax_id] = ancestors

        with open(args.output, "w") as f:
            for child, ancestors in child_ancestors.items():
                f.write(f"{child}\t{','.join(ancestors)}\n")


if __name__ == "__main__":
    main()
