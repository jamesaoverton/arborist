import sqlite3


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
            cur_new.execute("INSERT INTO statements VALUES " + insert)
            cur_new.execute("CREATE INDEX stanza_idx ON statements (stanza)")
            cur_new.execute("CREATE INDEX subject_idx ON statements (subject)")
            cur_new.execute("CREATE INDEX object_idx ON statements (object)")
            cur_new.execute("ANALYZE")


def get_child_ancestors(child_ancestors, child_parents, child, node):
    p = child_parents.get(node)
    if not p or p == node:
        return
    child_ancestors[child].add(p)
    get_child_ancestors(child_ancestors, child_parents, child, p)


def get_cumulative_counts(cur, count_map, child_ancestors):
    cuml_counts = {}
    for child, ancestors in child_ancestors.items():
        count = count_map.get(child, 0)
        if child in cuml_counts:
            cur_count = cuml_counts[child]
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


def get_descendants(cur, node, limits, descendants):
    # Get the children and maybe iterate
    cur.execute(
        """SELECT DISTINCT subject FROM statements
        WHERE object = ? AND predicate = 'rdfs:subClassOf'""",
        (node,),
    )
    for row in cur.fetchall():
        if row[0] in limits:
            return
        descendants.append(row[0])
        get_descendants(cur, row[0], limits, descendants)


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
