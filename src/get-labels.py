import csv
import logging
import re
import sqlite3
import string

from argparse import ArgumentParser
from helpers import get_curie


def bare(s):
    """Given a string, make it lowercase and remove useless bits
  so we can judge that labels are too similar."""
    parts = set(s.strip().lower().replace(string.punctuation, "").split())
    parts.discard("strain")
    parts.discard("str")
    return " ".join(sorted(parts))


def clean(s):
    return re.sub(r"\s+<.*>", "", s.strip())


def get_new_label(cur, labels, tax_id):
    # First check if tax_id is in labels (tax_id -> label)
    if tax_id in labels:
        return labels[tax_id]["Label"], "IEDB", labels[tax_id]["Synonyms"]

    cur.execute("SELECT value FROM statements WHERE stanza = ? AND predicate = 'rdfs:label'", (tax_id,))
    res = cur.fetchone()
    if not res:
        return None, None
    base_label = res[0]
    label = clean(base_label)
    bare_label = bare(label)
    source = "NCBI Taxonomy scientific name"

    # query for exact synonyms & their synonym types
    cur.execute(
        """SELECT s1.value, s2.object FROM statements s1
           JOIN statements s2 ON s1.stanza = s2.stanza
           WHERE s1.stanza = ?
           AND s1.predicate = 'oio:hasExactSynonym'
           AND s2.subject LIKE '_:%'
           AND s2.predicate = 'oio:hasSynonymType'""",
        (tax_id,)
    )
    exact_syns = {}
    for res in cur.fetchall():
        syn = clean(res[0])
        if bare(syn) == bare_label:
            continue
        source = res[1][10:].replace("_", " ")
        exact_syns[source] = syn

    # if exact syn has_synonym_type 'scientific name', override label with this
    if "scientific name" in exact_syns:
        label = clean(exact_syns["scientific name"])

    # if label starts with [ look for related synonym, override label with this
    if label.startswith("["):
        bare_label = bare(label)
        cur.execute("SELECT value FROM statements WHERE stanza = ? AND predicate = 'oio:hasRelatedSynonym'", (tax_id,))
        res = cur.fetchone()
        if res:
            syn = clean(res[0])
            if bare(syn) != bare_label:
                label = syn
                source = "NCBI Taxonomy equivalent name"

    common_name = exact_syns.get("common name")
    genbank_name = exact_syns.get("genbank common name")
    equivalent_name = exact_syns.get("equivalent name")

    if common_name and bare(common_name) != bare_label:
        # if exact syn has_synonym_type 'common name' append this to the label in parentheses
        label += f" ({common_name})"
        source = "NCBI Taxonomy scientific name (NCBI Taxonomy common name)"

    elif genbank_name and bare(genbank_name) != bare_label:
        # ... or for 'genbank common name'
        label += f" ({genbank_name})"
        source = "NCBI Taxonomy scientific name (GenBank common name)"

    elif equivalent_name and bare(equivalent_name) != bare_label:
        # 6. ... or for 'equivalent name'
        label += f" ({equivalent_name})"
        source = "NCBI Taxonomy scientific name (NCBI Taxonomy equivalent name)"

    if label != base_label:
        return label, source, ""
    return None, None, None


def main():
    parser = ArgumentParser()
    parser.add_argument("db", help="NCBITaxon database")
    parser.add_argument("labels", help="LJI SoT ncbi_taxa sheet with preferred labels")
    args = parser.parse_args()

    preferred_labels = {}
    with open(args.labels, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            preferred_labels[get_curie(row["Taxon ID"])] = {"Label": row["Label"], "Synonyms": row["IEDB Synonyms"]}

    new_labels = []
    with sqlite3.connect(args.db) as conn:
        cur = conn.cursor()
        cur.execute("SELECT DISTINCT stanza FROM statements WHERE stanza LIKE 'NCBITaxon:%' AND object = 'owl:Class'")
        for res in cur.fetchall():
            tax_id = res[0]
            label, source, synonyms = get_new_label(cur, preferred_labels, tax_id)
            if not label:
                continue
            new_labels.append([tax_id, label, source, synonyms])

    print("Taxon ID\tLabel\tLabel Source\tSynonyms")
    for detail in new_labels:
        print("\t".join(detail))


if __name__ == '__main__':
    main()
