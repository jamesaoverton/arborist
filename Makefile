### Workflow
#
# 1. [Install Requirements](install)
# 2. [Download SOT](refresh_sheets)
# 3. [Build Dependencies](browser_deps)
# 4. [View Browser](./src/browser.py?dbs=ncbi-trimmed-plus,ncbi-organized-plus,ncbi-pruned-plus,ncbi-override-plus,organism-tree&id=NCBITaxon:2759)

build:
	mkdir $@

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	RDFTAB_URL := https://github.com/ontodev/rdftab.rs/releases/download/v0.1.1/rdftab-x86_64-apple-darwin
else
	RDFTAB_URL := https://github.com/ontodev/rdftab.rs/releases/download/v0.1.1/rdftab-x86_64-unknown-linux-musl
endif

build/rdftab: | build
	curl -L -o $@ $(RDFTAB_URL)
	chmod +x $@

build/robot.jar: | build
	curl -Lk -o $@ https://build.obolibrary.io/job/ontodev/job/robot/job/master/lastSuccessfulBuild/artifact/bin/robot.jar

build/ncbitaxon.owl: | build
	curl -Lk http://purl.obolibrary.org/obo/ncbitaxon.owl > $@

build/ncbitaxon.db: src/prefixes.sql build/ncbitaxon.owl | build/rdftab
	rm -f $@
	sqlite3 $@ < src/prefixes.sql
	./build/rdftab $@ < $(word 2,$^)
	sqlite3 $@ "CREATE INDEX idx_stanza ON statements (stanza);"
	sqlite3 $@ "CREATE INDEX idx_subject ON statements (subject);"
	sqlite3 $@ "CREATE INDEX idx_predicate ON statements (predicate);"
	sqlite3 $@ "CREATE INDEX idx_object ON statements (object);"
	sqlite3 $@ "CREATE INDEX idx_value ON statements (value);"
	sqlite3 $@ "ANALYZE;"

build/organism-tree.owl: | build
	# TODO - download from ...

# Manual IEDB Sheets

# Label overrides
build/ncbi_taxa.tsv: | build
	curl -Lk -o $@ "https://docs.google.com/spreadsheets/d/1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4/export?format=tsv&id=1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4&gid=0"

# Parent overrides
build/taxon_parents.tsv: | build
	curl -Lk -o $@ "https://docs.google.com/spreadsheets/d/1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4/export?format=tsv&id=1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4&gid=175863322"

# Manual top-levels
build/top_level.tsv: | build
	curl -Lk -o $@ "https://docs.google.com/spreadsheets/d/1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4/export?format=tsv&id=1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4&gid=74523853"

IEDB_SHEETS := build/ncbi_taxa.tsv build/taxon_parents.tsv build/top_level.tsv
iedb_sheets: $(IEDB_SHEETS)

.PHONY: refresh_sheets
refresh_sheets:
	rm -rf $(IEDB_SHEETS)
	make iedb_sheets

# Tax ID -> epitope count
build/counts.tsv: src/get-counts.sql | build
	$(MIRROR_QUERY) < $< > $@

# Custom IEDB taxa
build/iedb_taxa.tsv: src/get-iedb-taxa.sql | build
	$(MIRROR_QUERY) < $< > $@

# not needed for new tree
build/active-taxa.tsv: build/counts.tsv
	cut -f1 $< > $@

# not needed for new tree
# Active taxa + manually updated taxa (labels & parents)
build/precious.tsv: build/ncbi_taxa.tsv build/iedb_taxa.tsv build/active-taxa.tsv
	tail -n +2 $< | cut -f1 > $@
	tail -n +2 $(word 2,$^) | cut -f1 >> $@
	cat $(word 3,$^) >> $@


### Trees


build/new-subspecies-tree.db: src/prefixes.sql src/run.py build/ncbitaxon.db build/counts.tsv build/iedb_taxa.tsv build/ncbi_taxa.tsv build/taxon_parents.tsv build/top_level.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@ || (rm -rf $@ && exit 1)


### Old Tasks


# IEDB active nodes (including IEDB taxa) + their ancestors (no pruning)
build/ncbi-trimmed.db: src/prefixes.sql src/trim.py build/ncbitaxon.db build/active-taxa.tsv build/iedb_taxa.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@ || (rm -rf $@ && exit 1)

# Get all label overrides based on NCBI synonyms
build/labels.tsv: src/get-labels.py build/ncbi-trimmed.db build/ncbi_taxa.tsv
	python3 $^ > $@

# ncbi-trimmed with manual changes
build/ncbi-override.db: src/prefixes.sql src/override.py build/ncbi-trimmed.db build/labels.tsv build/taxon_parents.tsv build/precious.tsv build/ncbi-trimmed-child-parents.tsv build/counts.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@ || (rm -rf $@ && exit 1)

# ncbi-trimmed organized with stable top levels
build/ncbi-organized.db: src/prefixes.sql src/organize.py build/ncbi-override.db build/top_level.tsv build/precious.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@ || (rm -rf $@ && exit 1)

# ncbi-organized with collapsed nodes based on weights
build/ncbi-pruned.db: src/prefixes.sql src/prune2.py build/ncbi-organized.db build/precious.tsv build/counts.tsv build/ncbi-organized-child-parents.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@ || (rm -rf $@ && exit 1)

# ncbi-pruned with thresholds to move species to "other" (1% of epitopes)
build/ncbi-rehomed.db: src/prefixes.sql src/rehome.py build/ncbi-pruned.db build/precious.tsv build/counts.tsv build/ncbi-pruned-child-parents.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@ || (rm -rf $@ && exit 1)

# Any database with counts
build/%-plus.db: src/prefixes.sql src/add-counts.py build/%.db build/counts.tsv build/%-child-parents.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@ || (rm -rf $@ && exit 1)


### Parent Maps

build/%-child-parents.tsv: src/get-child-parents.py build/%.db
	python3 $^ $@

build/organism-tree.db: src/prefixes.sql build/organism-tree.owl | build/rdftab
	rm -rf $@
	sqlite3 $@ < $<
	./build/rdftab $@ < $(word 2,$^)

build/subspecies-tree.db: src/prefixes.sql build/subspecies-tree.owl | build/rdftab
	rm -rf $@
	sqlite3 $@ < $<
	./build/rdftab $@ < $(word 2,$^)
	sqlite3 $@ "CREATE INDEX idx_stanza ON statements (stanza);"
	sqlite3 $@ "CREATE INDEX idx_subject ON statements (subject);"
	sqlite3 $@ "CREATE INDEX idx_predicate ON statements (predicate);"
	sqlite3 $@ "CREATE INDEX idx_object ON statements (object);"
	sqlite3 $@ "CREATE INDEX idx_value ON statements (value);"
	sqlite3 $@ "ANALYZE;"

.PHONY: install
install: requirements.txt
	python3 -m pip install -r $<

browser_deps: build/new-subspecies-tree-plus.db build/subspecies-tree-plus.db
