SHELL = /usr/bin/env bash

SNAKEFILE = workflow/Snakefile
PWD = $(shell pwd)

ISCAN_VERSION = 5.73-104.0
MINIFORGE_VERSION = 25.1.1-0

CORES = all
CACHE = $(PWD)/cache # Needs to be an absolute path

SETUP_CACHE = mkdir -p $(CACHE) &&\
              export SNAKEMAKE_OUTPUT_CACHE=$(CACHE)

SNAKEMAKE = $(SETUP_CACHE) &&\
            snakemake --cores $(CORES)\
                      --cache\
                      --printshellcmds

CONFIG = tests/config.yaml
GENOMES = tests/genomes.txt
RESULTS = tests/results

FIG_DIR = graphs
FIG_NAMES = dag filegraph rulegraph
SVGS = $(foreach i,$(FIG_NAMES),$(FIG_DIR)/$(i).svg)
PNGS = $(foreach i,$(FIG_NAMES),$(FIG_DIR)/$(i).png)

ISCAN_SCRIPT = utils/install_iscan.py
ISCAN_DATA = $(PWD)
ISCAN_DRY = --dry

RM_TEST = tests/rm_except_genomes.py

MINIFORGE_INSTALL_DIR = $(shell echo -n $$HOME)/miniforge3
SERVER = https://github.com/conda-forge/miniforge/releases/download/$(MINIFORGE_VERSION)
MINIFORGE = Miniforge3-$(MINIFORGE_VERSION)-Linux-x86_64.sh
LINK_MINIFORGE = $(SERVER)/$(MINIFORGE)
SHA256 = $(MINIFORGE).sha256
LINK_SHA256 = $(SERVER)/$(SHA256)

DEBUG = debug.py
CLEAN = .snakemake $(FIG_DIR) $(RESULTS) $(MINIFORGE) $(SHA256) $(CACHE) $(DEBUG)

R_LIBS_SCRIPT = utils/install_Rlibs.R


.PHONY help:
help:
	@awk -F':' '/^[a-zA-Z0-9_-]+:/ {print $$1}' Makefile


# Debugging print
print-%: ; @echo $* = $($*)


$(MINIFORGE):
	wget '$(LINK_MINIFORGE)'
	wget '$(LINK_SHA256)'
	sha256sum -c '$(SHA256)'


.PHONY install-mamba:
install-mamba: $(MINIFORGE)
	chmod +x $<
	./$< -b -u -p $(MINIFORGE_INSTALL_DIR)


.PHONY install-iscan:
install-iscan: $(ISCAN_SCRIPT)
	$< --reinstall --target $(ISCAN_VERSION) --data $(ISCAN_DATA) $(ISCAN_DRY)


.PHONY install-Rlibs:
install-Rlibs: $(R_LIBS_SCRIPT)
	$<


.PHONY test:
test: $(SNAKEFILE) $(GENOMES) $(CONFIG) $(RM_TEST)
	@printf "Before looking for errors, run:\n"
	@printf "make clean\n\n"
	$(RM_TEST)
	$(SNAKEMAKE) --configfile $(CONFIG)


.PHONY test-dry:
test-dry: $(SNAKEFILE) $(GENOMES) $(CONFIG) $(RM_TEST)
	$(RM_TEST)
	$(SNAKEMAKE) --configfile $(CONFIG) -np


.PHONY debug:
debug: $(SNAKEFILE) $(GENOMES) $(CONFIG)
	$(SNAKEMAKE) --configfile $(CONFIG) -np --print-compilation >| $(DEBUG)
	black $(DEBUG)
	less $(DEBUG)


.PHONY style:
style:
	snakefmt .
	black .
	/usr/bin/Rscript -e 'styler::style_dir("workflow")'
	/usr/bin/Rscript -e 'styler::style_dir("utils")'
	isort --float-to-top -- utils workflow workflow/Snakefile
	isort --float-to-top --ext smk -- utils workflow


$(SVGS): $(SNAKEFILE) $(GENOMES) $(CONFIG)
	mkdir -p $(FIG_DIR)
	$(SNAKEMAKE) --configfile $(CONFIG) --dag       | dot -Tsvg > $(FIG_DIR)/dag.svg
	$(SNAKEMAKE) --configfile $(CONFIG) --filegraph | dot -Tsvg > $(FIG_DIR)/filegraph.svg
	$(SNAKEMAKE) --configfile $(CONFIG) --rulegraph | dot -Tsvg > $(FIG_DIR)/rulegraph.svg


$(PNGS): $(SVGS) $(GENOMES) $(CONFIG)
	parallel convert -background none -size 6000x6000 {} {.}.png ::: $(SVGS)


report.html: $(SNAKEFILE) $(GENOMES) $(CONFIG)
	$(SNAKEMAKE) --configfile $(CONFIG) --report


.PHONY git-config:
git-config:
	git config --global alias.root 'rev-parse --show-toplevel'
	git config push.autoSetupRemote true


.PHONY clean:
clean:
	@rm -rf $(CLEAN)
	git clean -d -n
	@printf "\nTo remove untracked files:\n"
	@printf "    + git clean -d -f\n"
