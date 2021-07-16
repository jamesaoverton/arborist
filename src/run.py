import csv
import logging
import re
import sqlite3
import string
import sys

from argparse import ArgumentParser, FileType
from collections import defaultdict
from helpers import (
    clean_no_epitopes,
    create_other,
    get_all_ancestors,
    get_child_ancestors,
    get_child_parents,
    get_cumulative_counts,
    get_curie,
    get_descendants,
    get_descendants_and_ranks,
    move_precious_to_other,
    move_rank_to_other,
)


def add_iedb_taxa(cur, iedb_taxa):
    """Add IEDB taxa to the target database.

    :param cur: database connection cursor (target database)
    :param iedb_taxa: map of IEDB tax ID to details (label and parents)
    """
    for tax_id, details in iedb_taxa.items():
        label = details["Label"]
        cur.execute(
            """INSERT INTO statements (stanza, subject, predicate, object, value) VALUES
            (?, ?, "rdf:type", "owl:Class", null),
            (?, ?, "rdfs:label", null, ?);""",
            (tax_id, tax_id, tax_id, tax_id, label),
        )
        for pid in details["Parent IDs"]:
            cur.execute(
                """INSERT INTO statements (stanza, subject, predicate, object) VALUES
                (?, ?, "rdfs:subClassOf", ?);""",
                (tax_id, tax_id, pid),
            )
        if details["Synonyms"] != "NULL":
            for synonym in details["Synonyms"].split(";"):
                cur.execute(
                    """INSERT INTO statements (stanza, subject, predicate, value) VALUES
                    (?, ?, 'oio:hasExactSynonym', ?)""", (tax_id, tax_id, synonym))
        if details["Rank"] != "NULL":
            cur.execute(
                """INSERT INTO statements (stanza, subject, predicate, value) VALUES
                (?, ?, 'ONTIE:0003617', ?)""", (tax_id, tax_id, details["Rank"])
            )


def add_indexes(conn):
    cur = conn.cursor()
    cur.execute("CREATE INDEX idx_stanza ON statements (stanza)")
    cur.execute("CREATE INDEX idx_subject ON statements (subject)")
    cur.execute("CREATE INDEX idx_predicate ON statements (predicate)")
    cur.execute("CREATE INDEX idx_object ON statements (object)")
    cur.execute("CREATE INDEX idx_value ON statements (value)")
    cur.execute("ANALYZE")


def bare_label(s):
    """Given a string, make it lowercase and remove useless bits so we can judge that labels are too similar.

    :param s:
    :return:
    """
    parts = set(s.strip().lower().replace(string.punctuation, "").split())
    parts.discard("strain")
    parts.discard("str")
    return " ".join(sorted(parts))


def clean_label(s):
    """

    :param s:
    :return:
    """
    return re.sub(r"\s+<.*>", "", s.strip())


def fix_ranks(cur):
    """Replace ncbitaxon:has_rank with ONITE:0003617 (has taxonomic rank).

    :param cur: database connection cursor
    """
    cur.execute("SELECT DISTINCT stanza, object FROM statements WHERE predicate = 'ncbitaxon:has_rank'")
    insert = []
    for res in cur.fetchall():
        term_id = res[0]
        rank = res[1]
        if not rank.startswith("NCBITaxon:"):
            continue
        rank = rank.replace("NCBITaxon:", "")
        insert.append(", ".join([f"'{term_id}'", f"'{term_id}'", "'ONTIE:0003617'", f"'{rank}'"]))
    insert = ", ".join([f"({x})" for x in insert])
    cur.execute(f"INSERT INTO statements (stanza, subject, predicate, value) VALUES {insert}")
    cur.execute("DELETE FROM statements WHERE predicate = 'ncbitaxon:has_rank'")


def get_active_nodes(cur, active_taxa, iedb_taxa):
    active_tax_ids = set(active_taxa)
    # Add the parents of IEDB taxa
    for details in iedb_taxa.values():
        for pid in details["Parent IDs"]:
            active_tax_ids.add(pid)

    # print(f"Retrieving ancestors for {len(active_tax_ids)} active taxa...")

    all_taxa = set()

    for act_tax in active_tax_ids:
        all_taxa.add(act_tax)
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
            (act_tax, act_tax),
        )
        for row in cur.fetchall():
            parent_id = row[0]
            if parent_id == act_tax:
                continue
            all_taxa.add(parent_id)
    return all_taxa


def get_all_labels(conn, label_overrides):
    """Add automatically-chosen 'best' labels for all taxa in the database
    if they do not already have an IEDB label override.

    :param conn: database connection
    :param label_overrides: IEDB label overrides
    :return: IEDB label overrides + automatically picked best labels
    """
    new_labels = {}
    cur = conn.cursor()
    cur.execute(
        "SELECT DISTINCT stanza FROM statements WHERE stanza LIKE 'NCBITaxon:%' AND object = 'owl:Class'"
    )
    for res in cur.fetchall():
        tax_id = res[0]
        label, source, synonyms = get_best_label(cur, label_overrides, tax_id)
        if not label:
            continue
        new_labels[tax_id] = {
                "Label": label,
                "Label Source": source,
                "IEDB Synonyms": synonyms,
            }
    return new_labels


def get_best_label(cur, label_overrides, tax_id):
    """Select the best label for a taxon term based on its synonyms. If an IEDB label override is provided, use that.

    :param cur: database connection cursor
    :param label_overrides: IEDB label overrides
    :param tax_id: ID of taxon to select label for
    :return: best label, label source, synonyms
    """
    # First check if tax_id is in labels (tax_id -> label)
    if tax_id in label_overrides:
        return label_overrides[tax_id]["Label"], "IEDB", label_overrides[tax_id]["IEDB Synonyms"]

    cur.execute(
        "SELECT value FROM statements WHERE stanza = ? AND predicate = 'rdfs:label'", (tax_id,)
    )
    res = cur.fetchone()
    if not res:
        return None, None
    base_label = res[0]
    label = clean_label(base_label)
    bare = bare_label(label)
    source = "NCBI Taxonomy scientific name"

    # query for exact synonyms & their synonym types
    cur.execute(
        """SELECT s1.value, s2.object FROM statements s1
           JOIN statements s2 ON s1.stanza = s2.stanza
           WHERE s1.stanza = ?
           AND s1.predicate = 'oio:hasExactSynonym'
           AND s2.subject LIKE '_:%'
           AND s2.predicate = 'oio:hasSynonymType'""",
        (tax_id,),
    )
    exact_syns = {}
    for res in cur.fetchall():
        syn = clean_label(res[0])
        if bare_label(syn) == bare:
            continue
        source = res[1][10:].replace("_", " ")
        exact_syns[source] = syn

    # if exact syn has_synonym_type 'scientific name', override label with this
    if "scientific name" in exact_syns:
        label = clean_label(exact_syns["scientific name"])

    # if label starts with [ look for related synonym, override label with this
    if label.startswith("["):
        bare = bare_label(label)
        cur.execute(
            "SELECT value FROM statements WHERE stanza = ? AND predicate = 'oio:hasRelatedSynonym'",
            (tax_id,),
        )
        res = cur.fetchone()
        if res:
            syn = clean_label(res[0])
            if bare_label(syn) != bare:
                label = syn
                source = "NCBI Taxonomy equivalent name"

    common_name = exact_syns.get("common name")
    genbank_name = exact_syns.get("genbank common name")
    equivalent_name = exact_syns.get("equivalent name")

    if common_name and bare_label(common_name) != bare:
        # if exact syn has_synonym_type 'common name' append this to the label in parentheses
        label += f" ({common_name})"
        source = "NCBI Taxonomy scientific name (NCBI Taxonomy common name)"

    elif genbank_name and bare_label(genbank_name) != bare:
        # ... or for 'genbank common name'
        label += f" ({genbank_name})"
        source = "NCBI Taxonomy scientific name (GenBank common name)"

    elif equivalent_name and bare_label(equivalent_name) != bare:
        # 6. ... or for 'equivalent name'
        label += f" ({equivalent_name})"
        source = "NCBI Taxonomy scientific name (NCBI Taxonomy equivalent name)"

    if label != base_label:
        return label, source, ""
    return None, None, None


def get_collapse(cuml_counts, precious, child_parents, collapse, prev_nodes, threshold=0.99):
    """Get the full list of nodes to collapse between a lower-level and upper-level based on the
    threshold. We collapse terms that have less than the threshold percentage of epitopes and move
    the species to "other" for the upper-level node.

    :param cuml_counts: map of ID -> cumulative epitope count
    :param precious:
    :param child_parents: map of child -> parent
    :param collapse: list to add collapsed nodes to
    :param prev_nodes: previous nodes in this iteration
    :param threshold: threshold for percentage of epitopes
    """
    last_node = prev_nodes[-1]
    current_node = child_parents.get(last_node)
    if not current_node or current_node == "NCBITaxon:1" or current_node.endswith("other"):
        return

    prev_count = cuml_counts.get(last_node, 0)
    cur_count = cuml_counts.get(current_node, 0)
    if current_node not in precious and prev_count / cur_count > threshold:
        # Continue to check the next one until we hit one with less than threshold
        prev_nodes.append(current_node)
        get_collapse(cuml_counts, precious, child_parents, collapse, prev_nodes, threshold=threshold)

    # End with prev node
    if prev_nodes[0] != last_node:
        prev_nodes.append(last_node)
        collapse.append(prev_nodes)

    # Continue with this node as the start node
    prev_nodes = [current_node]
    get_collapse(cuml_counts, precious, child_parents, collapse, prev_nodes, threshold=threshold)


def get_start_nodes(cur, top_level, node, start_nodes):
    """

    :param cur:
    :param top_level:
    :param node:
    :param start_nodes:
    :return:
    """
    cur.execute(
        "SELECT stanza FROM statements WHERE predicate = 'rdfs:subClassOf' AND object = ?", (node,)
    )
    for res in cur.fetchall():
        term_id = res[0]
        if term_id in top_level:
            if top_level[term_id]["Child Rank"] != "manual":
                start_nodes.append(term_id)
            else:
                get_start_nodes(cur, top_level, term_id, start_nodes)


def get_top_level_line(top_structure, line, node):
    """Get a line of descendants from the top level structure.

    :param top_structure: map of node -> children
    :param line: line of descendants to add to
    :param node: current node
    """
    children = top_structure.get(node)
    if not children:
        return line
    line.extend(children)
    for c in children:
        get_top_level_line(top_structure, line, c)


def insert_taxa(source_conn, target_conn, active_taxa, iedb_taxa):
    """Insert triples about NCBI & IEDB taxa into a new database.
    This includes active taxa (with epitopes) and their ancestors.

    :param source_conn: database connection for NCBITaxon input
    :param target_conn: database connection for IEDB tree output
    :param active_taxa: list of active taxa IDs (with epitopes)
    :param iedb_taxa: map of IEDB tax ID to details (label and parents)
    """
    source_cur = source_conn.cursor()
    active_nodes = get_active_nodes(source_cur, active_taxa, iedb_taxa)
    # use f-string because we don't know how many values we have
    active_nodes = ", ".join([f"'{x}'" for x in active_nodes])
    source_cur.execute(f"SELECT * FROM statements WHERE stanza IN ({active_nodes})")

    # Create a list of INSERT statements & combine into one string
    insert = []
    for r in source_cur.fetchall():
        vals = []
        for itm in r:
            if not itm:
                vals.append("null")
            else:
                itm = itm.replace("'", "''")
                vals.append(f"'{itm}'")
        insert.append(", ".join(vals))
    insert = ", ".join([f"({x})" for x in insert])

    # Create tables in target database
    target_cur = target_conn.cursor()
    target_cur.execute(
        """CREATE TABLE statements (stanza TEXT,
                                    subject TEXT,
                                    predicate TEXT,
                                    object TEXT,
                                    value TEXT,
                                    datatype TEXT,
                                    language TEXT)"""
    )
    target_cur.execute("INSERT INTO statements VALUES " + insert)

    # Add the IEDB taxa
    add_iedb_taxa(target_cur, iedb_taxa)

    # Check for active taxa not in database
    target_cur.execute("SELECT DISTINCT stanza FROM statements WHERE object = 'owl:Class'")
    existing_ids = set([x[0] for x in target_cur.fetchall()])
    missing = set(active_taxa) - existing_ids
    if missing:
        logging.error(
            f"{len(missing)} active taxa missing from NCBITaxonomy:\n- " + "\n- ".join(missing)
        )


def move_up(cur, top_level_id, top_level_label, rank, precious=None, extras=None):
    # Init empty lists if these are not included
    if not precious:
        precious = []
    if not extras:
        extras = []

    # Get all descendants and their ranks under this top level
    child_parent = {}
    ranks = {}
    top_level = get_curie(top_level_id)
    get_descendants_and_ranks(cur, child_parent, ranks, top_level)

    # Find all nodes of the given rank
    at_rank = [x for x, y in ranks.items() if y == "NCBITaxon:" + rank]
    if not extras:
        extras = []
    at_rank.extend(extras)

    # Sometimes rank-level nodes may be under an extra/precious, make sure these aren't moved
    keep_in_place = set()
    # If a rank-level term has a precious ancestor,
    # we want to keep that ancestor as the direct child instead
    replace_as_child = {}
    for taxa in at_rank:
        ancestors = reversed(get_all_ancestors(child_parent, taxa, top_level))
        if not set(ancestors).isdisjoint(set(extras)) or not set(ancestors).isdisjoint(
            set(precious)
        ):
            # Remove if they are a descendant of an extra or precious node
            keep_in_place.add(taxa)
            # Find which to move directly under top level
            precious_ancestors = list(set(ancestors).intersection(set(precious)))
            if precious_ancestors:
                replace_as_child[taxa] = precious_ancestors[0]

    at_rank = set(at_rank) - keep_in_place
    for taxa, replacement in replace_as_child.items():
        if taxa in at_rank:
            at_rank.remove(taxa)
        at_rank.add(replacement)

    # Bump all nodes of given rank to top-level
    at_rank_str = ", ".join([f"'{x}'" for x in at_rank])
    cur.execute(
        f"""UPDATE statements SET object = '{top_level}'
        WHERE predicate = 'rdfs:subClassOf' AND subject IN ({at_rank_str})"""
    )

    # Find nodes to remove (ancestors to limit) - excluding at_rank under extras
    other_organisms = set()
    precious_others = set()
    for f in at_rank:
        if f not in child_parent:
            continue
        ancestors = get_all_ancestors(child_parent, f, top_level)
        if not ancestors:
            continue
        if not set(ancestors).isdisjoint(set(extras)) or not set(ancestors).isdisjoint(
            set(precious)
        ):
            # Extras & precious may not be of given rank
            # Make sure we don't accidentally remove them
            continue
        move = ancestors[-1]
        # Check for a descendant that is in precious and make sure to move it to 'other'
        descendants = []
        get_descendants(cur, move, [f], descendants)
        if not set(precious).isdisjoint(set(descendants)):
            for x in list(set(precious) & set(descendants)):
                precious_others.add(x)
        # Otherwise, move the last of the ancestors to 'other organism'
        other_organisms.add(move)

    # Find non-rank level nodes under top-level
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE predicate = 'rdfs:subClassOf' AND object = ?""",
        (top_level,),
    )
    non_at_rank = []
    for row in cur.fetchall():
        tax_id = row[0]
        if tax_id in extras or tax_id in precious:
            continue
        r = ranks.get(tax_id, "")
        if r != "NCBITaxon:" + rank:
            non_at_rank.append(tax_id)

    # Move all species-level or precious terms to "Other" then delete ancestors
    if non_at_rank or precious_others:
        move_precious_to_other(cur, precious, top_level_id, top_level_label, non_at_rank)
        # move_rank_to_other(cur, top_level_id, top_level_label, non_at_rank, precious=precious)

    o_str = ", ".join([f"'{x}'" for x in other_organisms])
    cur.execute(
        f"""UPDATE statements SET object = 'iedb-taxon:0100026-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({o_str})"""
    )


def organize(conn, top_level, precious):
    cur = conn.cursor()
    for curie, details in top_level.items():
        parent = get_curie(details["Parent ID"])

        # First, rehome this node
        cur.execute(
            "UPDATE statements SET object = ? WHERE subject = ? AND predicate = 'rdfs:subClassOf'",
            (parent, curie),
        )

        rank = details.get("Child Rank", "").strip()
        if rank == "":
            # Nothing to do, everything stays as-is
            continue

        if rank == "manual":
            # Everything NOT in this set gets moved to other
            cur.execute(
                """SELECT DISTINCT subject FROM statements
                WHERE object = ? AND predicate = 'rdfs:subClassOf'""",
                (curie,),
            )
            others = []
            for row in cur.fetchall():
                if row[0] not in top_level:
                    others.append(row[0])
            other_rank = details.get("Other Rank", "")
            if other_rank.strip() == "":
                other_rank = "species"
            move_rank_to_other(
                cur, details["ID"], details["Label"], others, rank=other_rank, precious=precious
            )
            continue

        # Otherwise, move all of given rank to top-level
        extras = [
            get_curie(x) for x in details.get("Extra Nodes", "").split(", ") if x.strip() != ""
        ]
        move_up(cur, details["ID"], details["Label"], rank, extras=extras, precious=precious)


def override(conn, label_overrides, parent_overrides):
    """

    :param conn:
    :param label_overrides:
    :param parent_overrides:
    :return:
    """
    cur = conn.cursor()

    # Replace cellular organisms with OBI organism
    cur.execute("UPDATE statements SET stanza = 'OBI:0100026' WHERE stanza = 'NCBITaxon:131567';")
    cur.execute("UPDATE statements SET subject = 'OBI:0100026' WHERE subject = 'NCBITaxon:131567';")
    cur.execute(
        """UPDATE statements SET value = 'Organism'
        WHERE stanza = 'OBI:0100026'
          AND subject = 'OBI:0100026'
          AND predicate = 'rdfs:label';"""
    )
    cur.execute("UPDATE statements SET object = 'OBI:0100026' WHERE object = 'NCBITaxon:131567';")
    bnode_id = 1

    # Override labels with IEDB labels
    for tax_id, row in label_overrides.items():
        new_label = row["Label"]
        label_source = row["Label Source"]
        cur.execute(
            """UPDATE statements SET value = ?
            WHERE stanza = ? AND subject = ? AND predicate = 'rdfs:label'""",
            (new_label, tax_id, tax_id),
        )
        # Add label source as annotation on label
        cur.execute(
            """INSERT INTO statements (stanza, subject, predicate, value) VALUES
            (?, ?, 'oio:hasLabelSource', ?)""",
            (tax_id, f"_:bnode{bnode_id:08}", label_source),
        )
        bnode_id += 1
        # Add synonyms if they exist
        for syn in row["IEDB Synonyms"].split(","):
            if syn.strip() == "":
                continue
            cur.execute(
                """INSERT INTO statements (stanza, subject, predicate, value) VALUES
                (?, ?, "oio:hasExactSynonym", ?)""",
                (tax_id, tax_id, syn),
            )

    # Override parents
    for tax_id, row in parent_overrides.items():
        parent_id = get_curie(row["Parent ID"])
        cur.execute(
            """
            UPDATE statements
            SET object = ?
            WHERE stanza = ? AND subject = ? AND predicate = 'rdfs:subClassOf';
        """,
            (parent_id, tax_id, tax_id),
        )


def parse_top_level(top_level_file):
    # ID -> Details
    top_level_unordered = {}
    # Parent -> Children
    top_structure = defaultdict(set)
    reader = csv.DictReader(top_level_file, delimiter="\t")
    for row in reader:
        tax_id = get_curie(row["ID"])
        parent_id = get_curie(row["Parent ID"])
        if parent_id not in top_structure:
            top_structure[parent_id] = set()
        top_structure[parent_id].add(tax_id)
        top_level_unordered[tax_id] = row

    # Sort top level by structure (starting with children of cellular organism)
    orgs = top_structure["OBI:0100026"]
    full_line = []
    for o in orgs:
        line = [o]
        get_top_level_line(top_structure, line, o)
        full_line.extend(line)
    non_orgs = top_structure["NCBITaxon:1"]
    for no in non_orgs:
        line = [no]
        get_top_level_line(top_structure, line, no)
        full_line.extend(line)

    # Go from lowest -> highest level
    full_line.reverse()
    return {node: top_level_unordered[node] for node in full_line}


def get_collapse_terms(counts, precious, parent_child, current_node, collapse_terms: list, threshold=0.99):
    collapse_terms.append(current_node)

    if current_node in precious:
        # if current node has epitopes, we collapse to the current node
        return

    parent = parent_child.get(current_node)
    if not parent:
        return

    # Check this number of epitopes vs. its parent's number of epitopes
    current_count = counts.get(current_node, 0)
    parent_count = counts.get(parent, 0)
    if current_count == parent_count:
        get_collapse_terms(counts, precious, parent_child, parent, collapse_terms, threshold=threshold)


def collapse(cur, counts, precious, top_level, child_parents, current_node, threshold=0.99):
    # Go up until we find a parent that does not have > 99% of epitopes
    collapse_terms = []
    get_collapse_terms(counts, precious, child_parents, current_node, collapse_terms, threshold=threshold)

    if len(collapse_terms) > 1:
        first_node = collapse_terms[0]
        last_node = collapse_terms[-1]
        parent_node = child_parents.get(last_node)
        if parent_node:
            cur.execute(
                "UPDATE statements SET object = ? WHERE stanza = ? AND predicate = 'rdfs:subClassOf'",
                (parent_node, first_node),
            )
            if parent_node not in top_level:
                collapse(cur, counts, precious, top_level, child_parents, parent_node, threshold=threshold)
        else:
            cur.execute(
                "UPDATE statements SET object = ? WHERE stanza = ? AND predicate = 'rdfs:subClassOf'",
                (last_node, first_node),
            )
    else:
        parent_node = child_parents.get(current_node)
        if parent_node and parent_node not in top_level:
            collapse(cur, counts, precious, top_level, child_parents, parent_node, threshold=threshold)


def prune(conn, counts, top_level, threshold=0.99):
    cur = conn.cursor()
    child_parents = get_child_parents(cur)
    child_ancestors = defaultdict(set)
    for child in child_parents.keys():
        if child not in child_ancestors:
            child_ancestors[child] = set()
        get_child_ancestors(child_ancestors, child_parents, child, child)
    cuml_counts = get_cumulative_counts(counts, child_ancestors)

    precious = set(counts.keys())
    precious.update(set(top_level.keys()))

    # Start from bottom nodes
    cur.execute(
        """SELECT DISTINCT stanza FROM statements
           WHERE object = 'owl:Class'
           AND stanza NOT IN (SELECT object FROM statements
            WHERE predicate = 'rdfs:subClassOf')"""
    )
    for res in cur.fetchall():
        s = res[0]
        if counts.get(s, 0) == 0:
            continue
        collapse(cur, cuml_counts, precious, top_level, child_parents, s, threshold=threshold)


def rehome(conn, counts, precious, top_level, threshold=0.01):
    """

    :param conn:
    :param counts:
    :param precious:
    :param top_level:
    :param threshold:
    :return:
    """
    cur = conn.cursor()
    child_parents = get_child_parents(cur)
    child_ancestors = defaultdict(set)
    for child in child_parents.keys():
        if child not in child_ancestors:
            child_ancestors[child] = set()
        get_child_ancestors(child_ancestors, child_parents, child, child)
    cuml_counts = get_cumulative_counts(counts, child_ancestors)
    # Make sure we aren't rehoming for anything we gave manual structure to
    # TODO: this should get the manual structure nodes from top-level sheet
    # - start at top level and then go down until we find non-manual node
    # - that is where we want to start rehoming
    # - skip rehoming for terms that have other children that are precious and go to next level
    for taxa in [
        "NCBITaxon:2",  # bacterium
        "NCBITaxon:10239",  # virus
        "NCBITaxon:4751",  # fungus
        "NCBITaxon:58024",  # spermatophyte
        "NCBITaxon:6854",  # arachnid
        "NCBITaxon:6657",  # crustacean
        "NCBITaxon:50557",  # insect
        "NCBITaxon:6447",  # mollusc
        "NCBITaxon:6231",  # nematode
        "NCBITaxon:6157",  # platyhelminth
    ]:
        # Check for direct 'other' children that are precious
        # This is a weird scenario because we may end up moving
        children = []
        flag = False
        cur.execute(
            "SELECT stanza FROM statements WHERE predicate = 'rdfs:subClassOf' AND object = ?",
            (taxa,),
        )
        for res in cur.fetchall():
            child_id = res[0]
            if child_id.endswith("-other"):
                if child_id in precious:
                    flag = True
            else:
                # Only rehome under non-other nodes
                children.append(child_id)

        if flag:
            for c in children:
                rehome_children(cur, precious, cuml_counts, c, threshold=threshold)
        else:
            rehome_children(cur, precious, cuml_counts, taxa, threshold=threshold)


def rehome_children(cur, precious, cuml_counts, parent_id, threshold=0.01):
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

    parent_count = cuml_counts[parent_id]
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
        count = cuml_counts.get(term_id, 0)
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
        # print("Special - " + parent_id)
        # Jump to rehoming the next level down
        for term_id in under_threshold:
            rehome_children(cur, precious, cuml_counts, term_id, threshold=threshold)
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
            rehome_children(cur, precious, cuml_counts, term_id, threshold=threshold)
    if others:
        # Find precious nodes in descendants of new other terms & move these to the other node
        move_precious_to_other(cur, precious, parent_id.split(":")[1], parent_label, others)


def main():
    parser = ArgumentParser()
    parser.add_argument("ncbitaxonomy", help="Path to NCBITaxonomy SQLite database")
    parser.add_argument(
        "counts",
        help="Path to counts TSV from IEDB (active taxa & eptiope counts)",
        type=FileType("r"),
    )
    parser.add_argument("iedb_taxa", help="Path to iedb_taxa TSV", type=FileType("r", encoding="latin1"))
    parser.add_argument(
        "ncbi_taxa", help="Path to ncbi_taxa TSV with updated labels", type=FileType("r")
    )
    parser.add_argument(
        "taxon_parents", help="Path to taxon_parents TSV with manual parents", type=FileType("r")
    )
    parser.add_argument(
        "top_level",
        help="Path to top_level TSV with stable top level structure",
        type=FileType("r"),
    )
    parser.add_argument("output", help="Output database")
    args = parser.parse_args()

    # Read in counts
    counts = {}
    reader = csv.reader(args.counts, delimiter="\t")
    next(reader)
    for row in reader:
        curie = get_curie(row[0])
        if curie == "NCBITaxon:694009":
            # Manually remove old COVID term from counts, we don't want this in the tree
            continue
        counts[curie] = int(row[1])

    # Read in custom IEDB taxa
    iedb_taxa = {}
    reader = csv.reader(args.iedb_taxa, delimiter="\t")
    next(reader)
    for row in reader:
        iedb_taxa[get_curie(row[0])] = {"Label": row[1], "Parent IDs": [get_curie(x) for x in row[2].split(",")], "Rank": row[3], "Synonyms": row[4]}

    # Read in label overrides from ncbi_taxa sheet
    reader = csv.DictReader(args.ncbi_taxa, delimiter="\t")
    label_overrides = {}
    for row in reader:
        label_overrides[get_curie(row["Taxon ID"])] = row

    # Read in manual parents from taxon_parents sheet
    reader = csv.DictReader(args.taxon_parents, delimiter="\t")
    parent_overrides = {}
    for row in reader:
        parent_overrides[get_curie(row["Taxon ID"])] = row

    # Read in stable top level
    top_level = parse_top_level(args.top_level)

    precious = []
    precious.extend(counts.keys())
    precious.extend(label_overrides.keys())

    with sqlite3.connect(args.output) as target_conn:
        # Copy the taxa from source to target (and add IEDB taxa)
        with sqlite3.connect(args.ncbitaxonomy) as source_conn:
            print("Inserting taxa into new database...")
            insert_taxa(source_conn, target_conn, precious, iedb_taxa)

        target_conn.execute(
            """INSERT INTO statements (stanza, subject, predicate, object, value) VALUES
            ('iedb-taxon:0100026-other', 'iedb-taxon:0100026-other', 'rdf:type', 'owl:Class', null),
            ('iedb-taxon:0100026-other', 'iedb-taxon:0100026-other', 'rdfs:label', null, 'Other')"""
        )

        # Add indexes
        print("Adding indexes...")
        add_indexes(target_conn)

        # Update label overrides to include best labels from synonyms
        print("Retrieving new labels...")
        label_overrides = get_all_labels(target_conn, label_overrides)

        # Override hierarchy with manual labels and parents
        print("Adding IEDB overrides...")
        override(target_conn, label_overrides, parent_overrides)

        # Organize hierarchy with stable top level
        print("Organizing stable top level...")
        organize(target_conn, top_level, precious)

        # Prune unnecessary intermediate nodes based on epitope percentage threshold (>99%)
        print("Pruning intermediate nodes...")
        prune(target_conn, counts, top_level)

        # Rehome nodes to "other" based on epitope percentage threshold (<1%)
        print("Moving nodes to 'other'...")
        rehome(target_conn, counts, precious, top_level)

        # Get updated child->ancestors
        cur = target_conn.cursor()
        child_parents = get_child_parents(cur)
        child_ancestors = defaultdict(set)
        for child in child_parents.keys():
            if child not in child_ancestors:
                child_ancestors[child] = set()
            get_child_ancestors(child_ancestors, child_parents, child, child)

        # Use child->ancestors to get updated cumulative epitope counts
        cuml_counts = get_cumulative_counts(counts, child_ancestors)

        # Clean up zero-epitope terms
        print("Cleaning zero-epitope terms...")
        clean_no_epitopes(cur, cuml_counts)

        # Replace ncbitaxon:has_rank with ONTIE property
        fix_ranks(cur)

        target_cur = target_conn.cursor()
        # Check for active taxa not in database
        target_cur.execute("SELECT DISTINCT stanza FROM statements WHERE object = 'owl:Class'")
        existing_ids = set([x[0] for x in target_cur.fetchall()])
        missing = set(counts.keys()) - existing_ids
        if missing:
            logging.error(
                f"{len(missing)} active taxa missing from NCBITaxonomy:\n- " + "\n- ".join(missing)
            )


if __name__ == "__main__":
    main()
