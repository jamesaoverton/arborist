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

        child_parent = {}
        print(f"Getting ancestors for {len(ids)} taxa...")
        for tax_id in ids:
            cur.execute(
                """WITH RECURSIVE ancestors(parent, child) AS (
                VALUES (?, NULL)
                UNION
                -- The non-blank parents of all of the parent terms extracted so far:
                SELECT object AS parent, subject AS child
                FROM statements, ancestors
                WHERE ancestors.parent = statements.stanza
                  AND statements.predicate = 'rdfs:subClassOf'
                  AND statements.object NOT LIKE '_:%'
                )
                SELECT * FROM ancestors""",
                (tax_id,),
            )
            for row in cur.fetchall():
                if not row[1]:
                    continue
                child_parent[row[1]] = row[0]

        with open(args.output, "w") as f:
            for child, parent in child_parent.items():
                f.write(f"{child}\t{parent}\n")


if __name__ == "__main__":
    main()
