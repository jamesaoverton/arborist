import csv
import sqlite3

from argparse import ArgumentParser
from collections import defaultdict
from helpers import (
    clean_no_epitopes,
    clean_others,
    copy_database,
    get_child_ancestors,
    get_child_parents,
    get_count_map,
    get_cumulative_counts,
    get_curie,
    move_precious_to_other
)


def rehome(cur, precious, count_map, parent_id, threshold=0.01):
    # Check if all children are already species/subspecies
    cur.execute(
        """SELECT DISTINCT s2.object FROM statements s1
        JOIN statements s2 ON s1.stanza = s2.stanza
        WHERE s1.predicate = 'rdfs:subClassOf'
        AND s1.object = ?
        AND s2.predicate = 'ncbitaxon:has_rank'""",
        (parent_id,),
    )
    child_ranks = [x[0] for x in cur.fetchall()]
    if "NCBITaxon:species" in child_ranks:
        child_ranks.remove("NCBITaxon:species")
    if "NCBITaxon:subspecies" in child_ranks:
        child_ranks.remove("NCBITaxon:subspecies")
    if not child_ranks:
        # Continue, do not go to next level because it won't change
        return

    parent_count = count_map[parent_id]
    # Get parent label
    cur.execute(
        "SELECT value FROM statements WHERE stanza = ? AND predicate = 'rdfs:label'", (parent_id,)
    )
    res = cur.fetchone()
    if res:
        parent_label = res[0]
    else:
        parent_label = parent_id
    # Get direct children of parent ID
    cur.execute(
        "SELECT stanza FROM statements WHERE predicate = 'rdfs:subClassOf' AND object = ?",
        (parent_id,),
    )
    others = []

    for res in cur.fetchall():
        term_id = res[0]
        if term_id.endswith("other"):
            continue
        count = count_map.get(term_id, 0)
        # TODO - what about if it only has one epitope?
        try:
            count / parent_count
        except ZeroDivisionError as e:
            print(term_id)
            print(parent_id)
            raise e
        if count / parent_count < threshold and not term_id.startswith("iedb-taxon"):
            # Move all species to other
            others.append(term_id)
        else:
            # Go to next level
            rehome(cur, precious, count_map, term_id, threshold=threshold)
    if others:
        move_precious_to_other(cur, precious, parent_id.split(":")[1], parent_label, others)


def main():
    parser = ArgumentParser()
    parser.add_argument("db")
    parser.add_argument("precious")
    parser.add_argument("counts")
    parser.add_argument("child_parents")
    parser.add_argument("output")
    args = parser.parse_args()

    precious = []
    with open(args.precious, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            precious.append(get_curie(row[0]))

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

    count_map = get_count_map(args.counts)

    copy_database(args.db, args.output)
    with sqlite3.connect(args.output) as conn:
        cur = conn.cursor()
        cuml_counts = get_cumulative_counts(count_map, child_ancestors)
        # Make sure we aren't rehoming for anything we gave manual structure to
        for taxa in [
            "NCBITaxon:2",
            "NCBITaxon:10239",
            "NCBITaxon:4751",
            "NCBITaxon:58024",
            "NCBITaxon:6854",
            "NCBITaxon:6657",
            "NCBITaxon:50557",
            "NCBITaxon:6447",
            "NCBITaxon:6231",
            "NCBITaxon:6157",
        ]:
            rehome(cur, precious, cuml_counts, taxa)

        # Get the child-ancestors again
        child_parents = get_child_parents(cur)
        child_ancestors = defaultdict(set)
        for child in child_parents.keys():
            if child not in child_ancestors:
                child_ancestors[child] = set()
            get_child_ancestors(child_ancestors, child_parents, child, child)

        cuml_counts = get_cumulative_counts(count_map, child_ancestors)
        clean_no_epitopes(cur, cuml_counts, precious)
        clean_others(cur, precious)



if __name__ == "__main__":
    main()
