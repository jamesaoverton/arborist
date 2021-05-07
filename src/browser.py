#!/usr/bin/env python

import json
import os
import re
import sqlite3
import urllib.parse

from collections import defaultdict
from gizmos import hiccup, tree, search
from jinja2 import Template


# Look for list of database files
# For each db file, create a column in HTML page
# Only search first db for the term (in search bar) - for now, merge results later

# ?dbs=x,y,z&id=foo:bar


browsers = {
    "ncbitaxon": {
        "name": "NCBITaxonomy",
        "description": "The official NCBITaxonomy"
    },
    "ncbi-trimmed-plus": {
        "name": "Trimmed NCBI",
        "description": "NCBI with only active IEDB nodes"
    },
    "ncbi-override-plus": {
        "name": "Overriden NCBI",
        "description": "Manual IEDB labels & parents"
    },
    "ncbi-organized-plus": {
        "name": "Organized NCBI",
        "description": "Stable top level organization"
    },
    "ncbi-pruned-plus": {
        "name": "Pruned NCBI",
        "description": "Collapsed intermediate nodes"
    },
    "ncbi-rehomed-plus": {
        "name": "Rehomed NCBI",
        "description": "Terms relocated to 'other' when parents have <1% of epitopes"
    },
    "organism-tree": {
        "name": "Organism Tree",
        "description": "The latest version of the organism tree used for IEDB"
    },
    "subspecies-tree-plus": {
        "name": "Subspecies Tree",
        "description": "The latest version of the subspecies tree provided for IEDB"
    }
}


def get_data(treename, cur, prefixes, term_id, stanza):
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
            children = [row[0] for row in cur.fetchall()]
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
        labels[row[0]] = row[1]
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
        obsolete.append(row[0])

    # If the compact URIs in the labels map are also in the tree, then add the label info to the
    # corresponding node in the tree:
    for key in hierarchy.keys():
        if key in labels:
            hierarchy[key]["label"] = labels[key]

    # Initialise a map with one entry for the tree and one for all of the labels corresponding to
    # all of the compact URIs in the stanza:
    return {
        "labels": labels,
        "obsolete": obsolete,
        treename: hierarchy,
        "iri": ontology_iri,
        "entity_type": entity_type,
    }


def get_rdfa(treename, cur, prefixes, data, href, stanza, term_id):
    ontology_iri = data.get("iri")
    labels = data["labels"]

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
        stanza = tree.create_stanza(cur.fetchall())
        subject = ontology_iri
        subject_label = data["labels"].get(ontology_iri, ontology_iri)
        si = tree.curie2iri(prefixes, subject)
    elif term_id != "ontology":
        subject = term_id
        si = tree.curie2iri(prefixes, subject)
        subject_label = label

    return tree.term2tree(data, treename, term_id, data["entity_type"], href=href)


def get_tree_html(treename, cur, prefixes, data, href, term, stanza, search=False):
    # Create the prefix element
    pref_strs = []
    for prefix, base in prefixes:
        pref_strs.append(f"{prefix}: {base}")
    pref_str = "\n".join(pref_strs)

    # get tree rdfa hiccup vector
    body = [get_rdfa(treename, cur, prefixes, data, href, stanza, term)]
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
    if treename in browsers:
        description = browsers[treename]["description"]
        treename = browsers[treename]["name"]
    else:
        description = ""
    body_wrapper.append(["h2", treename])
    body_wrapper.append(["div", {"style": "height: 30px;"}, ["small", description]])
    body = body_wrapper + body
    body = hiccup.render(prefixes, body, href=href)
    return body


def clean_value(value):
    if value.startswith("<"):
        value = value.lstrip("<").rstrip(">")
    if value.startswith("http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="):
        value = value.replace("http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=", "")
        value = f'<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={value}">{value}</a>'
    return value


def get_annotations(treename, cur, prefixes, data, href, term_id, stanza):
    # The subjects in the stanza that are of type owl:Axiom:
    annotation_bnodes = set()
    for row in stanza:
        if row["predicate"] == "owl:annotatedSource":
            annotation_bnodes.add(row["subject"])

    annotations = defaultdict(dict)
    for row in stanza:
        # subject is the blank node, _:...
        subject = row["subject"]
        if subject not in annotation_bnodes:
            continue

        if subject not in annotations:
            annotations[subject] = {}

        predicate = row["predicate"]
        obj = row["object"]
        value = row["value"]

        if predicate not in [
            "owl:annotatedSource",
            "owl:annotatedTarget",
            "owl:annotatedProperty",
            "rdf:type",
        ]:
            # This is the actual axiom that we care about and contains display value
            annotations[subject]["predicate"] = predicate
            if obj:
                annotations[subject]["object"] = obj
            if value:
                annotations[subject]["value"] = value
            annotations[subject]["annotation"] = row

        if predicate == "owl:annotatedSource":
            annotations[subject]["source"] = obj

        elif predicate == "owl:annotatedProperty":
            annotations[subject]["target_predicate"] = obj

        elif predicate == "owl:annotatedTarget":
            if obj:
                annotations[subject]["target_object"] = obj
            if value:
                annotations[subject]["target_value"] = value

    spv2annotation = {}
    for bnode, details in annotations.items():
        source = details["source"]
        target_predicate = details["target_predicate"]
        target = details.get("target_object", None) or details.get("target_value", None)

        if source in spv2annotation:
            # list of predicate -> values on this target (combo of predicate + value)
            pred2val = spv2annotation[source]
        else:
            pred2val = {}

        if target_predicate in pred2val:
            annotated_values = pred2val[target_predicate]
        else:
            annotated_values = {}

        if target in annotated_values:
            ax_annotations = annotated_values[target]
        else:
            ax_annotations = {}

        # predicate of the annotation
        ann_predicate = details["predicate"]
        if ann_predicate in ax_annotations:
            # values of the annotation
            anns = ax_annotations[ann_predicate]
        else:
            anns = []
        anns.append(details["annotation"])

        ax_annotations[ann_predicate] = anns
        annotated_values[target] = ax_annotations
        pred2val[target_predicate] = annotated_values
        spv2annotation[source] = pred2val

    # s2 maps the predicates of the given term to their corresponding rows (there can be more than
    # one row per predicate):
    s2 = defaultdict(list)
    for row in stanza:
        if row["subject"] == term_id:
            s2[row["predicate"]].append(row)
    pcs = list(s2.keys())

    labels = data["labels"]

    # Loop through the rows of the stanza that correspond to the predicates of the given term:
    predicate_values = {}
    for predicate in pcs:
        if "browser-link" in predicate:
            continue
        # Initialise an empty list of "o"s, i.e., hiccup representations of objects:
        objs = []
        for row in s2[predicate]:
            # Convert the `data` map, that has entries for the tree and for a list of the labels
            # corresponding to all of the curies in the stanza, into a hiccup object `o`:
            o = ["p", tree.row2o(stanza, data, row)]

            # Check for axiom annotations and create nested
            nest = tree.build_nested(
                treename, data, labels, spv2annotation, term_id, row, [], href=href
            )
            if nest:
                o += nest

            # Append the `o` to the list of `os`:
            objs.append(o)
        if objs:
            predicate_values[predicate] = hiccup.render(prefixes, objs.pop(0), href=href)

    predicates = predicate_values.keys()
    predicate_str = ",".join([f"'{x}'" for x in predicates])

    predicate_labels = {}
    cur.execute(
        f"SELECT DISTINCT subject, value FROM statements WHERE subject IN ({predicate_str}) AND predicate = 'rdfs:label'"
    )
    for row in cur.fetchall():
        predicate_labels[row[0]] = row[1]

    rem_predicates = set(set(predicates) - set(predicate_labels.keys()))
    for rp in rem_predicates:
        predicate_labels[rp] = rp

    return predicate_values, predicate_labels


def build_annotations(annotations, predicate_labels):
    table_names = list(annotations.keys())
    html = '<table class="table table-striped"><thead><th></th>'
    for tn in table_names:
        html += f"<th>{tn}</th>"
    html += "</thead><tbody>"
    for predicate, label in predicate_labels.items():
        html += f"<tr><td>{label}</td>"
        first_value = annotations[table_names[0]].get(predicate, "")
        first_iri = None
        if "href=" in first_value:
            first_iri = re.search(r'href="(.+)"', first_value).group(1)
        for tn in table_names:
            value = annotations[tn].get(predicate, "")
            td_class = ""
            if "href=" in value:
                value_iri = re.search(r'href="(.+)"', value).group(1)
                if not first_iri or value_iri != first_iri:
                    td_class = "table-warning"
            elif first_value != value:
                td_class = "table-warning"
            html += f'<td class="{td_class}">{value}</td>'
        html += "</tr>"
    html += "</tbody></table>"
    return html


def main():
    if "QUERY_STRING" in os.environ:
        args = dict(urllib.parse.parse_qsl(os.environ["QUERY_STRING"]))
    else:
        print("Content-Type: text/html")
        print("")
        print("Missing query parameters")
        return

    if "dbs" not in args:
        print("Content-Type: text/html")
        print("")
        print("One or more 'dbs' are required in query parameters")
        return

    dbs = args["dbs"].split(",")

    if "format" in args and args["format"] == "json":
        if not args.get("text"):
            print("Content-Type: application/json")
            print("")
            print(json.dumps([]))
            return
        # Search text in each database
        json_list = []
        search_text = urllib.parse.unquote(args["text"])
        for db in dbs:
            with sqlite3.connect(f"../build/{db}.db") as conn:
                json_list.extend(json.loads(search.search(conn, search_text)))
        # Sort alphabetically by length & name and take the first 20 results
        json_list = sorted(json_list, key=lambda i: i["label"])
        json_list = sorted(json_list, key=lambda i: (-len(i["label"]), i["label"]))[
          :20
        ]
        print("Content-Type: application/json")
        print("")
        print(json.dumps(json_list))
        return

    term = "owl:Class"
    if "id" in args:
        term = args["id"]

    href = "?dbs=" + args["dbs"] + "&id={curie}"

    annotations = {}
    predicate_labels = {}

    trees = []
    for db in dbs:
        with sqlite3.connect(f"../build/{db}.db") as conn:
            cur = conn.cursor()
            # Get prefixes
            cur.execute("SELECT * FROM prefix ORDER BY length(base) DESC")
            all_prefixes = [(x[0], x[1]) for x in cur.fetchall()]
            try:
                if term == "owl:Class":
                    stanza = []
                else:
                    cur.execute(f"SELECT * FROM statements WHERE stanza = '{term}'")
                    stanza = tree.create_stanza(cur.fetchall())

                if term != "owl:Class" and not stanza:
                    trees.append(f"<div><h2>{db}</h2><p>Term not found</p></div>")
                    continue

                data = get_data(db, cur, all_prefixes, term, stanza)
                trees.append(get_tree_html(db, cur, all_prefixes, data, href, term, stanza))
                if term and term not in top_levels:
                    predicate_values, cur_predicate_labels = get_annotations(
                        db, cur, all_prefixes, data, href, term, stanza
                    )
                    annotations[db] = predicate_values
                    for predicate, label in cur_predicate_labels.items():
                        if predicate in predicate_labels:
                            if predicate_labels[predicate] == predicate and predicate != label:
                                predicate_labels[predicate] = label
                        else:
                            predicate_labels[predicate] = label

            except Exception as e:
                raise e
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

    html = t.render(trees=trees, title="test", annotations=ann_html)

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


if __name__ == "__main__":
    main()
