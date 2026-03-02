# fkcompute Documentation

This folder contains the project documentation as standalone Markdown files (similar to `../../../snake_box/docs`).

If you just want to use the tool:

- Start here: `docs/quickstart.md`
- Then: `docs/cli.md`, `docs/config.md`, `docs/python_api.md`, `docs/output.md`

If you are developing/maintaining the project:

- Start here: `docs/overview.md`, `docs/pipeline.md`
- Then: `docs/performance.md`, `docs/inversion.md`, `docs/constraint_system.md`, `docs/ilp.md`, `docs/cpp_backend.md`
- Finally: `docs/development.md`, `docs/troubleshooting.md`

## Build A Single PDF (Optional)

Like `snake_box`, this repo includes a small Pandoc helper script to produce one combined PDF with a TOC.

Requirements:

- `pandoc`
- A TeX engine (`tectonic` recommended; `xelatex`/`lualatex`/`pdflatex` also work)

Build:

```bash
./docs/build_pdf.sh
# or choose the output path
./docs/build_pdf.sh ./docs/fkcompute_docs.pdf
```

Environment overrides (optional):

- `FKCOMPUTE_DOCS_HEADER_TEX` (defaults to `docs/pandoc_pdf_header.tex`)
- `FKCOMPUTE_DOCS_MAINFONT`, `FKCOMPUTE_DOCS_SANSFONT`, `FKCOMPUTE_DOCS_MONOFONT`
