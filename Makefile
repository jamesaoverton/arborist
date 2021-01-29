### Workflow
#
# 1. [Install Requirements](install)
# 2. [Download SOT](refresh_sheets)
# 3. [Build Dependencies](browser_deps)
# 4. [View Browser](./src/browser.py?dbs=ncbitaxon,ncbi-trimmed,organism-tree&id=NCBITaxon:2759)

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

refresh_sheets:
	rm -rf $(IEDB_SHEETS)
	make iedb_sheets

build/iedb-ncbitaxon.db: build/ncbitaxon.db src/update-ncbitaxon.py $(IEDB_SHEETS)
	cp $< $@
	python3 $(word 2,$^) $@ $(IEDB_SHEETS) || (rm -f $@ && exit 1)

build/active-taxa.tsv: | build
	# TODO - build from IEDB (see current org tree steps to create all-active-taxa)

.PHONY: build/ncbi-trimmed.db
build/ncbi-trimmed.db: build/iedb-ncbitaxon.db src/trim.py build/active-taxa.tsv
	cp $< $@
	python3 $(word 2,$^) $@ $(word 3,$^)

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

.PHONY: browser_deps
browser_deps: build/ncbitaxon.db build/iedb-ncbitaxon.db build/organism-tree.db
