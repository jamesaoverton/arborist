### Workflow
#
# 1. [Install Requirements](install)
# 2. [Download SOT](refresh_sheets)
# 3. [Build Dependencies](browser_deps)
# 4. [View Browser](./src/browser.py?dbs=ncbitaxon,ncbi-organized-plus,organism-tree&id=NCBITaxon:2759)

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

build/taxdmp.zip: | build
	curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip

build/ncbitaxon.db: src/ncbitaxon-to-db.py build/taxdmp.zip
	rm -f $@
	sqlite3 $@ < src/prefixes.sql
	python3 $^ $@ || (rm -f $@ && exit 1)
	sqlite3 $@ "CREATE INDEX idx_stanza ON statements (stanza);"
	sqlite3 $@ "CREATE INDEX idx_subject ON statements (subject);"
	sqlite3 $@ "CREATE INDEX idx_predicate ON statements (predicate);"
	sqlite3 $@ "CREATE INDEX idx_object ON statements (object);"
	sqlite3 $@ "CREATE INDEX idx_value ON statements (value);"
	sqlite3 $@ "ANALYZE;"

# Intermediate to organism tree - NCBITaxon + IEDB updates

build/ncbi_taxa.tsv: | build
	curl -Lk -o $@ "https://docs.google.com/spreadsheets/d/1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4/export?format=tsv&id=1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4&gid=0"

build/iedb_taxa.tsv: | build
	curl -Lk -o $@ "https://docs.google.com/spreadsheets/d/1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4/export?format=tsv&id=1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4&gid=1409951076"

build/taxon_parents.tsv: | build
	curl -Lk -o $@ "https://docs.google.com/spreadsheets/d/1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4/export?format=tsv&id=1vza08DSUVEDn1i470RUHNulb_HQYzyczCzcojrH1WU4&gid=1849479413"

IEDB_SHEETS := build/ncbi_taxa.tsv build/iedb_taxa.tsv build/taxon_parents.tsv
iedb_sheets: $(IEDB_SHEETS)

.PHONY: refresh_sheets
refresh_sheets:
	rm -rf $(IEDB_SHEETS)
	make iedb_sheets

build/active-taxa.tsv: build/counts.tsv
	cut -f1 $< > $@

# IEDB active nodes + their ancestors (no pruning)
build/ncbi-trimmed.db: src/prefixes.sql src/trim.py build/ncbitaxon.db build/active-taxa.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 src/trim.py build/ncbitaxon.db build/active-taxa.tsv $@

# active nodes + manual changes
build/ncbi-manual.db: src/prefixes.sql src/update-ncbitaxon.py build/ncbi-trimmed.db $(IEDB_SHEETS)
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@

# Tax ID -> epitope count
build/counts.tsv: | build
	# TODO - query from IEDB

# Tax ID -> parent
build/child-parents.tsv: src/get-child-parents.py build/ncbi-manual.db
	python3 $^ $@

# Active taxa + manually updated taxa
build/precious.tsv: build/ncbi_taxa.tsv build/active-taxa.tsv
	tail -n +2 $< | cut -f1 > $@
	cat $(word 2,$^) >> $@

# Collapse some nodes based on weights
build/ncbi-pruned.db build/pruned-counts.tsv: src/prefixes.sql src/prune.py build/ncbi-manual.db build/precious.tsv build/counts.tsv build/child-parents.tsv
	rm -rf build/ncbi-pruned.db
	sqlite3 build/ncbi-pruned.db < $<
	python3 $(filter-out src/prefixes.sql,$^) build/pruned-counts.tsv build/ncbi-pruned.db

# Pruned with epitope counts
build/ncbi-pruned-plus.db: src/prefixes.sql src/add-counts.py build/ncbi-pruned.db build/pruned-counts.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@ -c

# Organized with stable top levels
build/ncbi-organized.db: src/prefixes.sql src/organize.py build/ncbi-pruned.db build/precious.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 $(filter-out src/prefixes.sql,$^) $@

# The parents have changed, so create a new parent->child map
build/organized-child-parents.tsv: src/get-child-parents.py build/ncbi-organized.db
	python3 $^ $@

# Organized with epitope counts
build/ncbi-organized-plus.db: src/prefixes.sql src/add-counts.py build/ncbi-organized.db build/organized-child-parents.tsv build/counts.tsv
	rm -rf $@
	sqlite3 $@ < $<
	python3 src/add-counts.py build/ncbi-organized.db build/counts.tsv $@ -p build/organized-child-parents.tsv

build/organism-tree.owl: | build
	# TODO - download from ...

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

browser_deps: build/ncbitaxon.db build/ncbi-pruned-plus.db build/ncbi-organized-plus.db build/organism-tree.db
