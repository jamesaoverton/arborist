import csv
import sqlite3

from argparse import ArgumentParser
from collections import defaultdict
from helpers import copy_database, get_child_ancestors, get_cumulative_counts, get_curie


def main():
    parser = ArgumentParser()
    parser.add_argument("db")
    parser.add_argument("counts")
    parser.add_argument("child_parents")
    parser.add_argument("output")
    args = parser.parse_args()

    copy_database(args.db, args.output)
    with sqlite3.connect(args.output) as conn:
        cur = conn.cursor()
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

        count_map = {}
        with open(args.counts, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            for row in reader:
                if row[0] == "NULL":
                    continue
                count_map[get_curie(row[0])] = int(row[1])
        cuml_counts = get_cumulative_counts(count_map, child_ancestors)

        for tax_id, count in cuml_counts.items():
            cur.execute(
                f"""UPDATE statements SET value = value || ' ({count})'
                WHERE stanza = '{tax_id}'
                  AND subject = '{tax_id}'
                  AND predicate = 'rdfs:label';"""
            )


if __name__ == '__main__':
    main()
