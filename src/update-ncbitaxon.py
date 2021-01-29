#!/usr/bin/env python3

import csv
import sqlite3

from argparse import ArgumentParser


def get_id(tax):
    if tax == "OBI:0100026":
        return tax
    if tax.startswith("1") and len(tax) == 8:
        return "iedb-taxon:" + tax
    return "NCBITaxon:" + tax


def add_nodes(cur, nodes):
    with open(nodes, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tax_id = get_id(row["Taxon ID"])
            parent_id = get_id(row["Parent ID"])
            label = row["Label"]
            cur.execute("""
                INSERT INTO statements (stanza, subject, predicate, object, value) VALUES
                (?, ?, "rdf:type", "owl:Class", null),
                (?, ?, "rdfs:subClassOf", ?, null),
                (?, ?, "rdfs:label", null, ?);
            """, (tax_id, tax_id, tax_id, tax_id, parent_id, tax_id, tax_id, label))


def update_names(cur, names):
    # Replace cellular organisms with OBI organism
    cur.execute("""
        UPDATE statements
        SET stanza = "OBI:0100026"
        WHERE stanza = "NCBITaxon:131567";
    """)
    cur.execute("""
        UPDATE statements
        SET subject = "OBI:0100026"
        WHERE subject = "NCBITaxon:131567";
    """)
    cur.execute("""
        UPDATE statements
        SET value = "Organism"
        WHERE stanza = "OBI:0100026"
          AND subject = "OBI:0100026"
          AND predicate = "rdfs:label";
    """)
    cur.execute("""
        UPDATE statements
        SET object = "OBI:0100026"
        WHERE object = "NCBITaxon:131567";
    """)
    with open(names, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tax_id = get_id(row["Taxon ID"])
            new_label = row["Label"]
            cur.execute("""
                UPDATE statements
                SET value = ?
                WHERE stanza = ? AND subject = ? AND predicate = 'rdfs:label';
            """, (new_label, tax_id, tax_id))
            for syn in row["IEDB Synonyms"].split(","):
                if syn.strip() == "":
                    continue
                cur.execute("""
                    INSERT INTO statements (stanza, subject, predicate, value) VALUES
                    (?, ?, "oboInOwl:hasExactSynonym", ?);
                """, (tax_id, tax_id, syn))


def update_parents(cur, parents):
    with open(parents, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tax_id = get_id(row["Taxon ID"])
            parent_id = get_id(row["Parent ID"])
            cur.execute("""
                UPDATE statements
                SET object = ?
                WHERE stanza = ? AND subject = ? AND predicate = 'rdfs:subClassOf';
            """, (parent_id, tax_id, tax_id))


def update(source, target, names, nodes, parents):
    with sqlite3.connect(source) as conn:
        # Get stanzas from source database
        cur = conn.cursor()
        cur.execute(
            """SELECT stanza, subject, predicate, object, value, datatype, language FROM statements"""
        )
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
            cur_new.execute("""CREATE TABLE statements (stanza TEXT,
                                                        subject TEXT,
                                                        predicate TEXT,
                                                        object TEXT,
                                                        value TEXT,
                                                        datatype TEXT,
                                                        language TEXT)""")
            cur.execute(f"INSERT INTO statements VALUES " + insert)
            print("updating labels...")
            update_names(cur_new, names)
            print("adding IEDB taxa...")
            add_nodes(cur_new, nodes)
            print("updating parents...")
            update_parents(cur_new, parents)


def main():
    parser = ArgumentParser()
    parser.add_argument("db", type=str, help="Existing NCBITaxon database to update")
    parser.add_argument("name_overrides", type=str, help="Label overrides")
    parser.add_argument("iedb_taxa", type=str, help="IEDB taxa nodes to add")
    parser.add_argument("parent_overrides", type=str, help="Parent taxa overrides")
    parser.add_argument("output", type=str, help="Output database")
    args = parser.parse_args()
    update(args.db, args.output, args.name_overrides, args.iedb_taxa, args.parent_overrides)


if __name__ == '__main__':
    main()
