import sqlite3

from argparse import ArgumentParser
from helpers import get_child_parents


def get_curie(tax_id):
    if tax_id.startswith("OBI:"):
        return tax_id
    if len(tax_id) == 8 and tax_id.startswith("100"):
        return "iedb-taxon:" + tax_id
    return "NCBITaxon:" + tax_id


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="Database to get nodes from")
    parser.add_argument("output", help="Output TSV")
    args = parser.parse_args()

    with sqlite3.connect(args.db) as conn:
        cur = conn.cursor()
        child_parents = get_child_parents(cur)

        with open(args.output, "w") as f:
            for child, parent in child_parents.items():
                f.write(f"{child}\t{parent}\n")


if __name__ == "__main__":
    main()
