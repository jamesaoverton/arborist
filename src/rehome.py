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

    # First pass over results to check for special case:
    # where "other" node is the only node above threshold
    # In this case, if we rehome, "other" will be the only child of this node
    under_threshold = []
    all_children = []
    for res in cur.fetchall():
        term_id = res[0]
        all_children.append(term_id)
        count = count_map.get(term_id, 0)
        try:
            count / parent_count
        except ZeroDivisionError:
            print(term_id)
            continue
        if count / parent_count < threshold and not term_id.startswith("iedb-taxon"):
            under_threshold.append(term_id)

    # Remaining is everything that will not be moved to other
    remaining = list(set(all_children) - set(under_threshold))
    if len(remaining) == 1 and remaining[0].endswith("-other"):
        print("Special - " + parent_id)
        # Jump to rehoming the next level down
        for term_id in under_threshold:
            rehome(cur, precious, count_map, term_id, threshold=threshold)
    else:
        for term_id in under_threshold:
            # Move terms under threshold to other
            if term_id.endswith("other"):
                continue
            others.append(term_id)
        for term_id in remaining:
            # Go to next level for each not under threshold
            if term_id.endswith("other"):
                continue
            rehome(cur, precious, count_map, term_id, threshold=threshold)
    if others:
        # Find precious nodes in descendants of new other terms & move these to the other node
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
        # TODO: this should get the manual structure nodes from top-level sheet
        # - start at top level and then go down until we find non-manual node
        # - that is where we want to start rehoming
        # - skip rehoming for terms that have other children that are precious and go to next level
        for taxa in [
            "NCBITaxon:2",      # bacterium
            "NCBITaxon:10239",  # virus
            "NCBITaxon:4751",   # fungus
            "NCBITaxon:58024",  # spermatophyte
            "NCBITaxon:6854",   # arachnid
            "NCBITaxon:6657",   # crustacean
            "NCBITaxon:50557",  # insect
            "NCBITaxon:6447",   # mollusc
            "NCBITaxon:6231",   # nematode
            "NCBITaxon:6157",   # platyhelminth
        ]:
            # Check for direct 'other' children that are precious
            # This is a weird scenario because we may end up moving
            children = []
            rehome_children = False
            cur.execute(
                "SELECT stanza FROM statements WHERE predicate = 'rdfs:subClassOf' AND object = ?",
                (taxa,)
            )
            for res in cur.fetchall():
                child_id = res[0]
                if child_id.endswith("-other"):
                    if child_id in precious:
                        rehome_children = True
                else:
                    # Only rehome under non-other nodes
                    children.append(child_id)

            if rehome_children:
                for c in children:
                    rehome(cur, precious, cuml_counts, c)
            else:
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
