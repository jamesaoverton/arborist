#!/usr/bin/env python3

import csv
import sqlite3

from argparse import ArgumentParser


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


def add_counts(cur, count_map, child_ancestors):
    # print("Building counts for ancestors...")
    add_counts = {}
    for child, ancestors in child_ancestors.items():
        count = count_map.get(child, 0)
        if child in add_counts:
            cur_count = add_counts[child]
            add_counts[child] = cur_count + count
        else:
            add_counts[child] = count
        for a in ancestors:
            if a in add_counts:
                cur_count = add_counts[a]
                add_counts[a] = cur_count + count
            else:
                add_counts[a] = count

    # print(f"Adding counts to {len(add_counts)} nodes...")
    for tax_id, count in add_counts.items():
        cur.execute(
            f"""UPDATE statements SET value = value || ' ({count})'
            WHERE stanza = '{tax_id}'
              AND subject = '{tax_id}'
              AND predicate = 'rdfs:label';"""
        )


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="Database to add counts to")
    parser.add_argument("counts", help="TSV containing ID -> epitope count")
    parser.add_argument("child_ancestors", help="TSV containing ID -> ancestors")
    parser.add_argument("output", help="Output database")
    args = parser.parse_args()

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

            child_ancestors = {}
            with open(args.child_ancestors, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    child_ancestors[row[0]] = row[1].split(",")

            add_counts(cur_new, count_map, child_ancestors)


if __name__ == "__main__":
    main()
