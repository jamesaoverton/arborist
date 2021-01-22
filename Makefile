### Workflow
#
# 1. [Install Requirements](install)
# 2. [Build Dependencies](browser_deps)
# 1. [View Browser](./src/browser.py?dbs=organism-tree,subspecies-tree&id=NCBITaxon:9606)

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
	sqlite3 $@ "CREATE INDEX idx_object ON statements (object);"
	sqlite3 $@ "CREATE INDEX idx_value ON statements (value);"
	sqlite3 $@ "ANALYZE;"

build/ncbitaxon.owl.gz: | build
	curl -Lk http://purl.obolibrary.org/obo/ncbitaxon.owl | gzip > $@

build/ncbitaxon.db: src/prefixes.sql build/ncbitaxon.owl.gz | build/rdftab
	rm -rf $@
	sqlite3 $@ < $<
	zcat $(word 2,$^) | ./build/rdftab $@
	sqlite3 $@ "CREATE INDEX idx_stanza ON statements (stanza);"
	sqlite3 $@ "CREATE INDEX idx_object ON statements (object);"
	sqlite3 $@ "CREATE INDEX idx_value ON statements (value);"
	sqlite3 $@ "ANALYZE;"

.PHONY: install
install: requirements.txt
	python3 -m pip install -r $<

.PHONY: browser_deps
browser_deps: build/ncbitaxon.db build/organism-tree.db build/subspecies-tree.db
