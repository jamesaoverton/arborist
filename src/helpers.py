import csv
import sqlite3


def clean_no_epitopes(cur, counts, precious):
    # Get bottom-level terms (are not object of subclass statement)
    cur.execute(
        """SELECT DISTINCT stanza FROM statements WHERE stanza NOT IN
        (SELECT object FROM statements WHERE predicate = 'rdfs:subClassOf')"""
    )
    remove = set()
    for res in cur.fetchall():
        term_id = res[0]
        if term_id in precious or counts.get(term_id, 0) > 0:
            continue
        # TODO - Get the last ancestor that has no epitopes
        remove.add(get_term_to_remove(cur, counts, precious, term_id))
    for term_id in remove:
        cur.execute(
            """UPDATE statements SET object = 'iedb-taxon:0100026-other'
            WHERE predicate = 'rdfs:subClassOf' AND stanza = ?""",
            (term_id,),
        )


def clean_others(cur, precious):
    cur.execute(
        """SELECT DISTINCT stanza FROM statements
        WHERE stanza LIKE '%-other' AND stanza IS NOT 'iedb-taxon:0100026-other'"""
    )
    for res in cur.fetchall():
        other_id = res[0]
        if other_id in precious:
            # This 'other' node has epitopes, do nothing
            continue

        # Get the precious descendants before moving anything around
        precious_descendants = get_precious_descendants(cur, precious, other_id)

        # Check if this is the ONLY child of the parent class
        cur.execute(
            f"""SELECT s2.stanza, s2.object FROM statements s1 
            JOIN statements s2 ON s1.object = s2.object
            WHERE s1.stanza = ? AND s1.predicate = 'rdfs:subClassOf'
            AND s2.predicate = 'rdfs:subClassOf'""",
            (other_id,),
        )
        res = cur.fetchall()
        siblings = [x[0] for x in res]
        parent = res[0][1]
        if len(siblings) == 1:
            # Get rid of the other node and bump up all terms
            cur.execute(
                """UPDATE statements SET object = 'iedb-taxon:0100026-other'
                WHERE predicate = 'rdfs:subClassOf' AND stanza = ?""",
                (other_id,),
            )
            move_to = parent
        else:
            # Otherwise move the direct children to Other Organism
            move_to = other_id
            cur.execute(
                "SELECT stanza FROM statements WHERE object = ? AND predicate = 'rdfs:subClassOf'",
                (other_id,),
            )
            children = [x[0] for x in cur.fetchall()]
            child_str = ", ".join([f"'{x}'" for x in children])
            cur.execute(
                f"""UPDATE statements SET object = 'iedb-taxon:0100026-other'
                        WHERE predicate = 'rdfs:subClassOf' AND stanza in ({child_str})"""
            )

        # Move the precious descendants back to this other term OR its parent
        if precious_descendants:
            desc_str = ", ".join([f"'{x}'" for x in precious_descendants])
            cur.execute(
                f"""UPDATE statements SET object = '{move_to}'
                WHERE predicate = 'rdfs:subClassOf' AND stanza IN ({desc_str})"""
            )


def copy_database(input_db, output_db):
    with sqlite3.connect(input_db) as conn:
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

        # Insert all into target database then organize top levels
        with sqlite3.connect(output_db) as conn_new:
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
            cur_new.execute(f"INSERT INTO statements VALUES {insert}")
            cur_new.execute("CREATE INDEX stanza_idx ON statements (stanza)")
            cur_new.execute("CREATE INDEX subject_idx ON statements (subject)")
            cur_new.execute("CREATE INDEX object_idx ON statements (object)")
            cur_new.execute("ANALYZE")


def create_other(cur, parent_tax, parent_label):
    # Create other node if it does not exist
    other_id = f"iedb-taxon:{parent_tax}-other"
    cur.execute("SELECT * FROM statements WHERE stanza = ?", (other_id,))
    res = cur.fetchone()
    if not res:
        # Create a new class for "other" node
        parent_id = get_curie(parent_tax)
        other_label = "Other " + parent_label
        cur.execute(
            f"""INSERT INTO statements (stanza, subject, predicate, object, value) VALUES
                    (?, ?, 'rdfs:subClassOf', ?, null),
                    (?, ?, 'rdfs:label', null, ?);""",
            (other_id, other_id, parent_id, other_id, other_id, other_label),
        )


def get_all_ancestors(child_parent, node, limit, ancestors=None):
    if not ancestors:
        ancestors = []
    if node in child_parent:
        parent = child_parent[node]
        if parent == limit:
            return ancestors
        ancestors.append(parent)
        return get_all_ancestors(child_parent, parent, limit, ancestors=ancestors)
    return ancestors


def get_child_ancestors(child_ancestors, child_parents, child, node):
    p = child_parents.get(node)
    if not p or p == node:
        return
    child_ancestors[child].add(p)
    get_child_ancestors(child_ancestors, child_parents, child, p)


def get_child_parents(cur):
    ids = []
    cur.execute(
        """SELECT DISTINCT stanza FROM statements
        WHERE object = 'owl:Class'
        AND stanza NOT IN (SELECT object FROM statements
         WHERE predicate = 'rdfs:subClassOf')"""
    )
    for row in cur.fetchall():
        ids.append(row[0])

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
            parent_id = row[0]
            if not row[1]:
                continue
            child_parent[row[1]] = parent_id
    return child_parent


def get_count_map(f):
    count_map = {}
    with open(f, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row[0] == "NULL":
                continue
            count_map[get_curie(row[0])] = int(row[1])
    return count_map


def get_cumulative_counts(count_map, child_ancestors):
    cuml_counts = {}
    for child, ancestors in child_ancestors.items():
        count = count_map.get(child, 0)
        if child in cuml_counts:
            cur_count = cuml_counts.get(child, 0)
            cuml_counts[child] = cur_count + count
        else:
            cuml_counts[child] = count

        for a in ancestors:
            if a in cuml_counts:
                cur_count = cuml_counts[a]
                cuml_counts[a] = cur_count + count
            else:
                cuml_counts[a] = count
    return cuml_counts


def get_curie(tax_id):
    if tax_id.startswith("OBI:"):
        return tax_id
    if tax_id.endswith("other"):
        return "iedb-taxon:" + tax_id
    if len(tax_id) == 8 and tax_id.startswith("100"):
        return "iedb-taxon:" + tax_id
    return "NCBITaxon:" + tax_id


def get_descendants(cur, node, limits, descendants, only_limit=False):
    # Get the children and maybe iterate
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE object = ? AND predicate = 'rdfs:subClassOf'""",
        (node,),
    )
    for row in cur.fetchall():
        tax_id = row[0]
        if tax_id in limits:
            if only_limit:
                descendants.append(tax_id)
            return
        if not only_limit:
            descendants.append(tax_id)
        get_descendants(cur, tax_id, limits, descendants)


def get_descendants_and_ranks(cur, child_parent, ranks, node):
    # Get the rank of this node
    cur.execute(
        "SELECT object FROM statements WHERE predicate = 'ncbitaxon:has_rank' AND subject = ?",
        (node,),
    )
    res = cur.fetchone()
    if res:
        ranks[node] = res[0]
    # Get the children and maybe iterate
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE object = ? AND predicate = 'rdfs:subClassOf'""",
        (node,),
    )
    for row in cur.fetchall():
        child_parent[row[0]] = node
        get_descendants_and_ranks(cur, child_parent, ranks, row[0])


def get_precious_descendants(cur, precious, node):
    descendants = []
    get_descendants(cur, node, [], descendants)
    precious_descendants = set(descendants).intersection(precious)

    child_parent = {}
    ranks = {}
    get_descendants_and_ranks(cur, child_parent, ranks, node)

    # Determine if a term in species has a parent in species (maybe strain subclass of species?)
    # Do not double up, just move the parent
    remove = set()
    for p in precious_descendants:
        ancestors = get_all_ancestors(child_parent, p, node)
        if set(ancestors).intersection(precious_descendants):
            remove.add(p)
    return precious_descendants - remove


def get_term_to_remove(cur, counts, precious, term_id):
    cur.execute(
        "SELECT object FROM statements WHERE predicate = 'rdfs:subClassOf' AND stanza = ?",
        (term_id,),
    )
    res = cur.fetchone()
    if res:
        parent_id = res[0]
        if parent_id in precious or counts.get(parent_id, 0) > 0:
            return term_id
        return get_term_to_remove(cur, counts, precious, parent_id)
    return term_id


def move_rank_to_other(cur, parent_tax_id, parent_tax_label, others, rank="species", precious=None):
    create_other(cur, parent_tax_id, parent_tax_label)
    if rank == "none":
        # Just set the others to be children of "Other" and return
        others_str = ", ".join([f"'{x}'" for x in others])
        cur.execute(
            f"""UPDATE statements SET object = 'iedb-taxon:{parent_tax_id}-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({others_str})"""
        )
        return

    precious_others = None
    for o in others:
        other_id = f"iedb-taxon:{parent_tax_id}-other"
        if o in precious:
            # print(o + " is precious")
            # Move this node to other
            # then get it's at-rank (or precious) children and move those directly under it
            cur.execute(
                f"""UPDATE statements SET object = 'iedb-taxon:{parent_tax_id}-other'
                WHERE predicate = 'rdfs:subClassOf' AND stanza = '{o}';"""
            )
            other_id = o

        # Create map of parent -> child and ranks
        child_parent = {}
        ranks = {}
        get_descendants_and_ranks(cur, child_parent, ranks, o)

        # Find all at rank level
        at_rank = [x for x, y in ranks.items() if y == "NCBITaxon:" + rank]

        # Also add in any precious
        # (making sure not to remove any important species-subspecies relationships)
        precious_others = set(precious).intersection(set(child_parent.keys()))
        remove_from_po = set()
        for p in precious_others:
            ancestors = set(get_all_ancestors(child_parent, p, o))
            if ancestors.intersection(precious_others) or set(at_rank).intersection(ancestors):
                remove_from_po.add(p)
        precious_others = precious_others - remove_from_po
        at_rank.extend(precious_others)

        # These get bumped up to 'other' then all extra nodes get deleted
        at_rank_str = ", ".join([f"'{x}'" for x in at_rank])
        cur.execute(
            f"""UPDATE statements SET object = '{other_id}'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({at_rank_str})"""
        )

        if not other_id.endswith("other"):
            # Special clean up when the other node is not actually an "other" (o is precious)
            # Find the non-at-ranks and move to other organism
            cur.execute(
                "SELECT stanza FROM statements WHERE predicate = 'rdfs:subClassOf' AND object = ?",
                (other_id,),
            )
            other_organisms = [x[0] for x in cur.fetchall() if x[0] not in at_rank]
            o_str = ", ".join([f"'{x}'" for x in other_organisms])
            cur.execute(
                f"""UPDATE statements SET object = 'iedb-taxon:0100026-other'
                WHERE predicate = 'rdfs:subClassOf' AND subject IN ({o_str})"""
            )

        # Find nodes to move to 'other organism'
        other_organisms = set()
        for s in at_rank:
            if s not in child_parent:
                continue
            ancestors = get_all_ancestors(child_parent, s, "NCBITaxon:" + parent_tax_id)
            if not ancestors:
                continue
            other_organisms.add(ancestors[-1])
        if not other_organisms and o not in precious:
            other_organisms.add(o)
        other_organisms = other_organisms - set(precious)
        if other_organisms:
            o_str = ", ".join([f"'{x}'" for x in other_organisms])
            cur.execute(
                f"""UPDATE statements SET object = 'iedb-taxon:0100026-other'
                WHERE predicate = 'rdfs:subClassOf' AND subject IN ({o_str})"""
            )

    # Finally, move the 'precious' others to the correct level
    if precious_others:
        po_str = ", ".join([f"'{x}'" for x in precious_others])
        cur.execute(
            f"""UPDATE statements SET object = 'iedb-taxon:{parent_tax_id}-other'
            WHERE predicate = 'rdfs:subClassOf' AND subject IN ({po_str})"""
        )


def move_precious_to_other(cur, precious, parent_tax_id, parent_tax_label, others):
    create_other(cur, parent_tax_id, parent_tax_label)
    exclude_from_other_org = set()
    other_id = f"iedb-taxon:{parent_tax_id}-other"
    for o in others:
        if o in precious:
            exclude_from_other_org.add(o)
            continue
        precious_descendants = get_precious_descendants(cur, precious, o)

        # Add all species to those excluded from being totally removed
        exclude_from_other_org.update(precious_descendants)

        # Move species to this-level other node
        pd_str = ", ".join([f"'{x}'" for x in precious_descendants])
        cur.execute(
            f"""UPDATE statements SET object = '{other_id}'
            WHERE predicate = 'rdfs:subClassOf' AND stanza IN ({pd_str})"""
        )
    # Move all others to other organism now that their species are gone
    others = set(others) - exclude_from_other_org
    others_str = ", ".join([f"'{x}'" for x in others])
    cur.execute(
        f"""UPDATE statements set object = 'iedb-taxon:0100026-other'
        WHERE predicate = 'rdfs:subClassOf' AND stanza IN ({others_str})"""
    )
