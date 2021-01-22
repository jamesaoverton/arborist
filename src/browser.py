#!/usr/bin/env python

import os
import sqlite3
import urllib.parse

from gizmos import hiccup, tree, search
from jinja2 import Template


# Look for list of database files
# For each db file, create a column in HTML page
# Only search first db for the term (in search bar) - for now, merge results later

# ?dbs=x,y,z&id=foo:bar


def get_rdfa(treename, cur, prefixes, href, stanza, term_id):
    ontology_iri, ontology_title = tree.get_ontology(cur, prefixes)
    if term_id not in top_levels:
        # Get a hierarchy under the entity type
        entity_type = tree.get_entity_type(cur, term_id)
        hierarchy, curies = tree.get_hierarchy(cur, term_id, entity_type)
    else:
        # Get the top-level for this entity type
        entity_type = term_id
        if term_id == "ontology":
            hierarchy = {term_id: {"parents": [], "children": []}}
            curies = set()
            if ontology_iri:
                curies.add(ontology_iri)
        else:
            if term_id == "owl:Individual":
                tls = ", ".join([f"'{x}'" for x in top_levels.keys()])
                cur.execute(
                    f"""SELECT DISTINCT subject FROM statements
                    WHERE subject NOT IN
                        (SELECT subject FROM statements
                         WHERE predicate = 'rdf:type'
                         AND object NOT IN ('owl:Individual', 'owl:NamedIndividual'))
                    AND subject IN
                        (SELECT subject FROM statements
                         WHERE predicate = 'rdf:type' AND object NOT IN ({tls}))"""
                )
            elif term_id == "rdfs:Datatype":
                cur.execute(
                    """SELECT DISTINCT subject FROM statements
                    WHERE predicate = 'rdf:type' AND object = 'rdfs:Datatype'"""
                )
            else:
                pred = "rdfs:subPropertyOf"
                if term_id == "owl:Class":
                    pred = "rdfs:subClassOf"
                # Select all classes without parents and set them as children of owl:Thing
                cur.execute(
                    f"""SELECT DISTINCT subject FROM statements 
                    WHERE subject NOT IN 
                        (SELECT subject FROM statements
                         WHERE predicate = '{pred}'
                         AND object IS NOT 'owl:Thing')
                    AND subject IN 
                        (SELECT subject FROM statements 
                         WHERE predicate = 'rdf:type'
                         AND object = '{term_id}' AND subject NOT LIKE '_:%'
                         AND subject NOT IN ('owl:Thing', 'rdf:type'));"""
                )
            children = [row["subject"] for row in cur.fetchall()]
            hierarchy = {term_id: {"parents": [], "children": children}}
            curies = {term_id}
            for c in children:
                hierarchy[c] = {"parents": [term_id], "children": []}
                curies.add(c)

    # Add all of the other compact URIs in the stanza to the set of compact URIs:
    stanza.sort(key=lambda x: x["predicate"])
    for row in stanza:
        curies.add(row.get("subject"))
        curies.add(row.get("predicate"))
        curies.add(row.get("object"))
    curies.discard("")
    curies.discard(None)

    # Get all the prefixes that are referred to by the compact URIs:
    ps = set()
    for curie in curies:
        if not isinstance(curie, str) or len(curie) == 0 or curie[0] in ("_", "<"):
            continue
        prefix, local = curie.split(":")
        ps.add(prefix)

    # Get all of the rdfs:labels corresponding to all of the compact URIs, in the form of a map
    # from compact URIs to labels:
    labels = {}
    ids = "', '".join(curies)
    cur.execute(
        f"""SELECT subject, value
      FROM statements
      WHERE stanza IN ('{ids}')
        AND predicate = 'rdfs:label'
        AND value IS NOT NULL"""
    )
    for row in cur:
        labels[row["subject"]] = row["value"]
    for t, o_label in top_levels.items():
        labels[t] = o_label
    if ontology_iri and ontology_title:
        labels[ontology_iri] = ontology_title

    obsolete = []
    cur.execute(
        f"""SELECT DISTINCT subject
            FROM statements
            WHERE stanza in ('{ids}')
              AND predicate='owl:deprecated'
              AND value='true'"""
    )
    for row in cur:
        obsolete.append(row["subject"])

    # If the compact URIs in the labels map are also in the tree, then add the label info to the
    # corresponding node in the tree:
    for key in hierarchy.keys():
        if key in labels:
            hierarchy[key]["label"] = labels[key]

    # Initialise a map with one entry for the tree and one for all of the labels corresponding to
    # all of the compact URIs in the stanza:
    data = {"labels": labels, "obsolete": obsolete, treename: hierarchy, "iri": ontology_iri}

    # Determine the label to use for the given term id when generating RDFa (the term might have
    # multiple labels, in which case we will just choose one and show it everywhere). This defaults
    # to the term id itself, unless there is a label for the term in the stanza corresponding to the
    # label for that term in the labels map:
    if term_id in labels:
        selected_label = labels[term_id]
    else:
        selected_label = term_id
    label = term_id
    for row in stanza:
        predicate = row["predicate"]
        value = row["value"]
        if predicate == "rdfs:label" and value == selected_label:
            label = value
            break

    subject = None
    si = None
    subject_label = None
    if term_id == "ontology" and ontology_iri:
        cur.execute(
            f"""SELECT * FROM statements
            WHERE subject = '{ontology_iri}'"""
        )
        stanza = cur.fetchall()
        subject = ontology_iri
        subject_label = data["labels"].get(ontology_iri, ontology_iri)
        si = tree.curie2iri(prefixes, subject)
    elif term_id != "ontology":
        subject = term_id
        si = tree.curie2iri(prefixes, subject)
        subject_label = label

    return tree.term2tree(data, treename, term_id, entity_type, href=href)


def get_tree_html(cur, treename, href, term, search=False):
    # Get prefixes
    cur.execute("SELECT * FROM prefix ORDER BY length(base) DESC")
    all_prefixes = [(x["prefix"], x["base"]) for x in cur.fetchall()]

    if term == "owl:Class":
        stanza = []
    else:
        cur.execute(f"SELECT * FROM statements WHERE stanza = '{term}'")
        stanza = cur.fetchall()
        if not stanza:
            return f"<div><h2>{treename}</h2><p>Term not found</p></div>"

    # Create the prefix element
    pref_strs = []
    for prefix, base in all_prefixes:
        pref_strs.append(f"{prefix}: {base}")
    pref_str = "\n".join(pref_strs)
    
    # get tree rdfa hiccup vector
    body = [get_rdfa(treename, cur, all_prefixes, href, stanza, term)]
    body_wrapper = ["div", {"prefix": pref_str}]
    if search:
        body_wrapper.append(
        [
            "div",
            {"class": "form-row mt-2 mb-2"},
            [
                "input",
                {
                    "id": f"statements-typeahead",
                    "class": "typeahead form-control",
                    "type": "text",
                    "value": "",
                    "placeholder": "Search",
                },
            ],
        ]
    )
    body_wrapper.append(["h2", treename])
    body = body_wrapper + body
    return hiccup.render(all_prefixes, body, href=href)


def clean_value(value):
    if value.startswith("<"):
        value = value.lstrip("<").rstrip(">")
    if value.startswith("http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="):
        value = value.replace("http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=", "")
        value = f'<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={value}">{value}</a>'
    return value


def get_annotations(cur, term_id, href):
    predicate_values = {}
    cur.execute(f"SELECT DISTINCT predicate, value FROM statements WHERE stanza = '{term_id}' AND subject = '{term_id}' AND value IS NOT NULL;")
    for row in cur.fetchall():
        value = clean_value(row["value"])
        predicate_values[row["predicate"]] = value
    
    cur.execute(f"SELECT DISTINCT predicate, object FROM statements WHERE stanza = '{term_id}' AND subject = '{term_id}' AND object IS NOT NULL;")
    for row in cur.fetchall():
        obj_id = row["object"]
        href2 = href.replace("{curie}", obj_id)
        cur.execute(f"SELECT value FROM statements WHERE stanza = '{obj_id}' AND subject = '{obj_id}' AND predicate = 'rdfs:label';")
        label = cur.fetchone()
        if label:
            label = label["value"]
            if label:
                value = f'<a href="{href2}">{label}<a>'
        else:
            value = clean_value(obj_id)
            value = f'<a href="{href2}">{value}</a>'
        predicate_values[row["predicate"]] = value

    predicates = predicate_values.keys()
    predicate_str = ",".join([f"'{x}'" for x in predicates])

    predicate_labels = {}
    cur.execute(f"SELECT DISTINCT subject, value FROM statements WHERE subject IN ({predicate_str}) AND predicate = 'rdfs:label'")
    for row in cur.fetchall():
        predicate_labels[row["subject"]] = row["value"]

    rem_predicates = set(set(predicates) - set(predicate_labels.keys()))
    for rp in rem_predicates:
        predicate_labels[rp] = rp

    return predicate_values, predicate_labels


def build_annotations(annotations, predicate_labels):
    table_names = annotations.keys()
    html = '<table class="table table-striped"><thead><th></th>'
    for tn in table_names:
        html += f"<th>{tn}</th>"
    html += "</thead><tbody>"
    for predicate, label in predicate_labels.items():
        html += f"<tr><td>{label}</td>"
        for tn in table_names:
            value = annotations[tn].get(predicate, "")
            html += f"<td>{value}</td>"
        html += "</tr>"
    html += "</tbody></table>"
    return html


def main():
    if "QUERY_STRING" in os.environ:
        args = dict(urllib.parse.parse_qsl(os.environ["QUERY_STRING"]))
    else:
        print("Content-Type: text/html")
        print("")
        print("Bad page")
        return

    if "dbs" not in args:
        print("Content-Type: text/html")
        print("")
        print("Bad page")
        return

    dbs = args["dbs"].split(",")
    first_db = dbs.pop(0)

    if "format" in args and args["format"] == "json":
        # TODO - maybe we can search both & merge results?
        json = search.search(f"../build/{first_db}.db", args["text"])
        print("Content-Type: application/json")
        print("")
        print(json)
        return

    term = "owl:Class"
    if "id" in args:
        term = args["id"]

    href = "?dbs=" + args['dbs'] + "&id={curie}"

    annotations = {}
    predicate_labels = {}

    with sqlite3.connect(f"../build/{first_db}.db") as conn:
        conn.row_factory = tree.dict_factory
        cur = conn.cursor()
        try:
            first_html = get_tree_html(cur, first_db, href, term, search=True)
        except Exception as e:
            print("Content-Type: text/html")
            print("")
            print("Error when generating HTML for " + db + ":<br>" + str(e))
            return
        if term and term not in top_levels:
            predicate_values, cur_predicate_labels = get_annotations(cur, term, href)
            annotations[first_db] = predicate_values
            for predicate, label in cur_predicate_labels.items():
                if predicate in predicate_labels:
                    if predicate_labels[predicate] == predicate and predicate != label:
                        predicate_labels[predicate] = label
                else:
                    predicate_labels[predicate] = label

    trees = []
    for db in dbs:
        with sqlite3.connect(f"../build/{db}.db") as conn:
            conn.row_factory = tree.dict_factory
            cur = conn.cursor()
            try:
                trees.append(get_tree_html(cur, db, href, term))
                if term and term not in top_levels:
                    predicate_values, cur_predicate_labels = get_annotations(cur, term, href)
                    annotations[db] = predicate_values
                    for predicate, label in cur_predicate_labels.items():
                        if predicate in predicate_labels:
                            if predicate_labels[predicate] == predicate and predicate != label:
                                predicate_labels[predicate] = label
                        else:
                            predicate_labels[predicate] = label
            except Exception as e:
                print("Content-Type: text/html")
                print("")
                print("Error when generating HTML for " + db + ":<br>" + str(e))
                return

    # Load Jinja template with CSS & JS and left & right trees
    with open("index.html.jinja2", "r") as f:
        t = Template(f.read())

    if annotations:
        ann_html = build_annotations(annotations, predicate_labels)
    else:
        ann_html = ""

    html = t.render(first=first_html, trees=trees, title="test", annotations=ann_html)

    # Return with CGI headers
    print("Content-Type: text/html")
    print("")
    print(html)
    return


top_levels = {
    "ontology": "Ontology",
    "owl:Class": "Class",
    "owl:AnnotationProperty": "Annotation Property",
    "owl:DataProperty": "Data Property",
    "owl:ObjectProperty": "Object Property",
    "owl:Individual": "Individual",
    "rdfs:Datatype": "Datatype",
}


if __name__ == '__main__':
    main()
