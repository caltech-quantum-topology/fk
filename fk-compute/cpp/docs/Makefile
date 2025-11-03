# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD  ?= sphinx-build
SOURCEDIR    = source
BUILDDIR     = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Additional targets for development
clean:
	rm -rf $(BUILDDIR)/*

livehtml:
	sphinx-autobuild -b html $(SOURCEDIR) $(BUILDDIR)/html

watch:
	@echo "Starting live reload server..."
	@echo "Open http://localhost:8000 in your browser"
	sphinx-autobuild -b html --host 0.0.0.0 --port 8000 $(SOURCEDIR) $(BUILDDIR)/html

# Check for broken links
linkcheck:
	@$(SPHINXBUILD) -b linkcheck $(SOURCEDIR) $(BUILDDIR)/linkcheck

# Build PDF documentation (requires LaTeX)
latexpdf: latex
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

# Check documentation quality
doctest:
	@$(SPHINXBUILD) -b doctest $(SOURCEDIR) $(BUILDDIR)/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

# Build for deployment (optimized)
deploy: clean html
	@echo "Building optimized documentation for deployment..."
	@echo "Documentation built in $(BUILDDIR)/html"
	@echo "Ready for deployment to web server"

# Development helpers
serve:
	@echo "Serving documentation at http://localhost:8080"
	cd $(BUILDDIR)/html && python3 -m http.server 8080

# Quality checks
spelling:
	@$(SPHINXBUILD) -b spelling $(SOURCEDIR) $(BUILDDIR)/spelling
	@echo "Spell check finished; results in $(BUILDDIR)/spelling."

# Generate API documentation from source code
apidoc:
	@echo "Generating API documentation from source code..."
	# Note: This would require sphinx-apidoc for C++ if available
	@echo "Manual API documentation is maintained in source/api/"