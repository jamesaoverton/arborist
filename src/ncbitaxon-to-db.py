#!/usr/bin/env python3

import io
import sqlite3

from argparse import ArgumentParser
from collections import defaultdict
from datetime import date
from zipfile import ZipFile

oboInOwl = {
    "SynonymTypeProperty": "synonym_type_property",
    "hasAlternativeId": "has_alternative_id",
    "hasBroadSynonym": "has_broad_synonym",
    "hasDbXref": "database_cross_reference",
    "hasExactSynonym": "has_exact_synonym",
    "hasOBOFormatVersion": "has_obo_format_version",
    "hasOBONamespace": "has_obo_namespace",
    "hasRelatedSynonym": "has_related_synonym",
    "hasScope": "has_scope",
    "hasSynonymType": "has_synonym_type"
}

exact_synonym = "oboInOwl:hasExactSynonym"
related_synonym = "oboInOwl:hasRelatedSynonym"
broad_synonym = "oboInOwl:hasBroadSynonym"

predicates = {
    "acronym": broad_synonym,
    "anamorph": related_synonym,
    "blast name": related_synonym,
    "common name": exact_synonym,
    "equivalent name": exact_synonym,
    "genbank acronym": broad_synonym,
    "genbank anamorph": related_synonym,
    "genbank common name": exact_synonym,
    "genbank synonym": related_synonym,
    "in-part": related_synonym,
    "misnomer": related_synonym,
    "misspelling": related_synonym,
    "synonym": related_synonym,
    "scientific name": exact_synonym,
    "teleomorph": related_synonym,
}

ranks = [
    "class",
    "cohort",
    "family",
    "forma",
    "genus",
    "infraclass",
    "infraorder",
    "kingdom",
    "order",
    "parvorder",
    "phylum",
    "section",
    "series",
    "species group",
    "species subgroup",
    "species",
    "subclass",
    "subcohort",
    "subfamily",
    "subgenus",
    "subkingdom",
    "suborder",
    "subphylum",
    "subsection",
    "subspecies",
    "subtribe",
    "superclass",
    "superfamily",
    "superkingdom",
    "superorder",
    "superphylum",
    "tribe",
    "varietas",
]

nodes_fields = [
    "tax_id",  # node id in GenBank taxonomy database
    "parent_tax_id",  # parent node id in GenBank taxonomy database
    "rank",  # rank of this node (superkingdom, kingdom, ...)
    "embl_code",  # locus-name prefix; not unique
    "division_id",  # see division.dmp file
    "inherited_div_flag",  # (1 or 0) 1 if node inherits division from parent
    "genetic_code_id",  # see gencode.dmp file
    "inherited_GC_flag",  # (1 or 0) 1 if node inherits genetic code from parent
    "mitochondrial_genetic_code_id",  # see gencode.dmp file
    "inherited_MGC_flag",  # (1 or 0) 1 if node inherits mitochondrial gencode from parent
    "GenBank_hidden_flag",  # (1 or 0) 1 if name is suppressed in GenBank entry lineage
    "hidden_subtree_root_flag",  # (1 or 0) 1 if this subtree has no sequence data yet
    "comments",  # free-text comments and citations
]


def escape_literal(text):
    return text.replace("'", "''")


def label_to_id(text):
    return text.replace(" ", "_").replace("-", "_")


def split_line(line):
    """Split a line from a .dmp file"""
    return [x.replace("\t", " ").strip() for x in line.split("\t|")]


def create_tables(cur):
    cur.execute("""CREATE TABLE IF NOT EXISTS statements (
                     stanza TEXT,
                     subject TEXT,
                     predicate TEXT,
                     object TEXT,
                     value TEXT,
                     datatype TEXT,
                     language TEXT
                   );""")
    cur.execute("""
        INSERT OR IGNORE INTO statements (stanza, subject, predicate, object, value, datatype) VALUES
        ("rdfs:label", "rdfs:label", "rdf:type", "owl:AnnotationProperty", null, null),
        ("rdfs:comment", "rdfs:comment", "rdf:type", "owl:AnnotationProperty", null, null),
        ("obo:IAO_0000115", "obo:IAO_0000115", "rdf:type", "owl:AnnotationProperty", null, null),
        ("obo:IAO_0000115", "obo:IAO_0000115", "rdfs:label", null, "definition", "xsd:string"),
        ("ncbitaxon:has_rank", "ncbitaxon:has_rank", "rdf:type", "owl:AnnotationProperty", null, null),
        ("ncbitaxon:has_rank", "ncbitaxon:has_rank", "rdfs:label", null, "has_rank", "xsd:string"),
        ("oboInOwl:hasExactSynonym", "oboInOwl:hasExactSynonym", "rdf:type", "owl:AnnotationProperty", null, null),
        ("oboInOwl:hasRelatedSynonym", "oboInOwl:hasRelatedSynonym", "rdf:type", "owl:AnnotationProperty", null, null),
        ("oboInOwl:hasBroadSynonym", "oboInOwl:hasBroadSynonym", "rdf:type", "owl:AnnotationProperty", null, null),
        ("<http://purl.obolibrary.org/obo/NCBITaxon#_taxonomic_rank>", "<http://purl.obolibrary.org/obo/NCBITaxon#_taxonomic_rank>", "rdf:type", "owl:Class", null, null),
        ("<http://purl.obolibrary.org/obo/NCBITaxon#_taxonomic_rank>", "<http://purl.obolibrary.org/obo/NCBITaxon#_taxonomic_rank>", "rdfs:label", null, "taxonomic rank", "xsd:string");
    """)
    
    for predicate, label in oboInOwl.items():
        predicate = "oboInOwl:" + predicate
        cur.execute(f"""
            INSERT OR IGNORE INTO statements (stanza, subject, predicate, object, value, datatype) VALUES
            ("{predicate}", "{predicate}", "rdf:type", "owl:AnnotationProperty", null, null),
            ("{predicate}", "{predicate}", "rdfs:label", null, "{label}", "xsd:string");
        """)

    for label, parent in predicates.items():
        predicate = "ncbitaxon:" + label_to_id(label)
        parent = parent.replace("oboInOwl", "oio")
        query = f"""
            INSERT OR IGNORE INTO statements (stanza, subject, predicate, object, value, datatype) VALUES
            ("{predicate}", "{predicate}", "rdf:type", "owl:AnnotationProperty", null, null),
            ("{predicate}", "{predicate}", "rdfs:label", null, "{label}", "xsd:string"),
            ("{predicate}", "{predicate}", "oboInOwl:hasScope", null, "{parent}", "xsd:string"),
            ("{predicate}", "{predicate}", "rdfs:subPropertyOf", "oboInOwl:SynonymTypeProperty", null, null);
        """
        cur.execute(query)


def add_node(cur, node, label, merged, synonyms, citations):
    tax_id = node["tax_id"]
    curie = "NCBITaxon:" + tax_id
    label = escape_literal(label)
    query = f"""INSERT INTO statements (stanza, subject, predicate, object, value, datatype) VALUES
        ("{curie}", "{curie}", "rdf:type", "owl:Class", null, null),
        ("{curie}", "{curie}", "rdfs:label", null, '{label}', "xsd:string"),
        ("{curie}", "{curie}", "oboInOwl:hasOBONamespace", null, "ncbi_taxonomy", "xsd:string")"""

    parent_tax_id = node["parent_tax_id"]
    if parent_tax_id and parent_tax_id != "" and parent_tax_id != tax_id:
        query += f""",
        ("{curie}", "{curie}", "rdfs:subClassOf", "NCBITaxon:{parent_tax_id}", null, null)"""

    rank = node["rank"]
    if rank and rank != "" and rank != "no rank":
        # if rank not in ranks:
        # print(f"WARN: Unrecognized rank '{rank}'")
        rank = label_to_id(rank)
        if rank in ["species_group", "species_subgroup"]:
            value = f"<http://purl.obolibrary.org/obo/NCBITaxon#_{rank}>"
        else:
            value = "NCBITaxon:" + rank
        query += f""",
        ("{curie}", "{curie}", "ncbitaxon:has_rank", "{value}", null, null)"""

    gc_id = node["genetic_code_id"]
    if gc_id:
        query += f""",
        ("{curie}", "{curie}", "oboInOwl:hasDbXref", null, "GC_ID:{gc_id}", "xsd:string")"""

    for merge in merged:
        query += f""",
        ("{curie}", "{curie}", "oboInOwl:hasAlternativeId", null, "NCBITaxon:{merge}", "xsd:string")"""

    for pubmed_id in citations:
        query += f""",
        ("{curie}", "{curie}", "oboInOwl:hasDbXref", null, "PMID:{pubmed_id}", "xsd:string")"""

    query += ";"
    cur.execute(query)


def add_synonyms(cur, curie, synonyms):
    query_values = []
    for synonym, unique, name_class in synonyms:
        if name_class in predicates:
            synonym = escape_literal(synonyms)
            predicate = predicates[name_class]
            synonym_type = label_to_id(name_class)
            # TODO - add blank nodes
            query_values.append(f"""("{curie}", "{curie}", "{predicate}", '{synonym}', "xsd:string")""")
    if query_values:
        query = "INSERT INTO statements (stanza, subject, predicate, value, datatype) VALUES "
        query += ", ".join(query_values)
        query += ";"
        cur.execute(query)


def convert(taxdmp_path, database):
    scientific_names = defaultdict(list)
    labels = {}
    synonyms = defaultdict(list)
    merged = defaultdict(list)
    citations = defaultdict(list)
    with sqlite3.connect(database) as conn:
        cur = conn.cursor()
        create_tables(cur)

        with ZipFile(taxdmp_path) as taxdmp:
            # get names
            print("retrieving names...")
            with taxdmp.open("names.dmp") as dmp:
                for line in io.TextIOWrapper(dmp):
                    tax_id, name, unique, name_class, _ = split_line(line)
                    if name_class == "scientific name":
                        labels[tax_id] = name
                        scientific_names[name].append([tax_id, unique])
                    else:
                        synonyms[tax_id].append([name, unique, name_class])

            # use unique name only if there's a conflict
            for name, values in scientific_names.items():
                tax_ids = [x[0] for x in values]
                if len(tax_ids) > 1:
                    uniques = [x[1] for x in values]
                    if len(tax_ids) != len(set(uniques)):
                        print("WARN: Duplicate unique names", tax_ids, uniques)
                    for tax_id, unique in values:
                        labels[tax_id] = unique
                        synonyms[tax_id].append([name, unique, "scientific name"])

            print("retrieving merged...")
            with taxdmp.open("merged.dmp") as dmp:
                for line in io.TextIOWrapper(dmp):
                    old_tax_id, new_tax_id, _ = split_line(line)
                    merged[new_tax_id].append(old_tax_id)

            print("retrieving citations...")
            with taxdmp.open("citations.dmp") as dmp:
                for line in io.TextIOWrapper(dmp):
                    cit_id, cit_key, pubmed_id, medline_id, url, text, tax_id_list, _ = split_line(line)
                    if medline_id == "0":
                        continue
                    for tax_id in tax_id_list.split():
                        citations[tax_id].append(medline_id)

            print("adding nodes...")
            with taxdmp.open("nodes.dmp") as dmp:
                for line in io.TextIOWrapper(dmp):
                    node = {}
                    fields = split_line(line)
                    for i in range(0, min(len(fields), len(nodes_fields))):
                        node[nodes_fields[i]] = fields[i]
                    tax_id = node["tax_id"]
                    add_node(cur, node, labels[tax_id], merged[tax_id], synonyms[tax_id], citations[tax_id])

        for label in ranks:
            rank = label_to_id(label)
            if rank in ["species_group", "species_subgroup"]:
                iri = f"<http://purl.obolibrary.org/obo/NCBITaxon#_{rank}>"
            else:
                iri = "NCBITaxon:" + rank
            cur.execute(f"""
                INSERT INTO statements (stanza, subject, predicate, object, value, datatype) VALUES
                ("{iri}", "{iri}", "rdf:type", "owl:Class", null, null),
                ("{iri}", "{iri}", "rdfs:label", null, "{label}", "xsd:string"),
                ("{iri}", "{iri}", "rdfs:subClassOf", "<http://purl.obolibrary.org/obo/NCBITaxon#_taxonomic_rank>", null, null),
                ("{iri}", "{iri}", "oboInOwl:hasOBONamespace", null, "ncbi_taxonomy", "xsd:string");
            """)



def main():
    parser = ArgumentParser()
    parser.add_argument("taxdmp", type=str, help="The taxdmp.zip file to read")
    parser.add_argument("db", type=str, help="The output database to write to")
    args = parser.parse_args()

    convert(args.taxdmp, args.db)


if __name__ == '__main__':
    main()
