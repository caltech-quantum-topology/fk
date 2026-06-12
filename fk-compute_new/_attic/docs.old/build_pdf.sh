#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

OUT_PDF="${1:-${SCRIPT_DIR}/fkcompute_docs.pdf}"
HEADER_TEX="${FKCOMPUTE_DOCS_HEADER_TEX:-${SCRIPT_DIR}/pandoc_pdf_header.tex}"

if ! command -v pandoc >/dev/null 2>&1; then
  echo "error: pandoc not found on PATH" >&2
  echo "install: https://pandoc.org/installing.html" >&2
  exit 127
fi

PDF_ENGINE=""
if command -v tectonic >/dev/null 2>&1; then
  PDF_ENGINE="tectonic"
elif command -v xelatex >/dev/null 2>&1; then
  PDF_ENGINE="xelatex"
elif command -v lualatex >/dev/null 2>&1; then
  PDF_ENGINE="lualatex"
elif command -v pdflatex >/dev/null 2>&1; then
  PDF_ENGINE="pdflatex"
else
  echo "error: no TeX engine found (tectonic/xelatex/lualatex/pdflatex)" >&2
  echo "hint: install tectonic (recommended) or a LaTeX distribution" >&2
  exit 127
fi

title="fkcompute Docs"
date_str="$(date +%Y-%m-%d)"

# Fonts (only used for unicode-capable engines).
# Override via env vars if desired.
MAIN_FONT="${FKCOMPUTE_DOCS_MAINFONT:-DejaVu Serif}"
SANS_FONT="${FKCOMPUTE_DOCS_SANSFONT:-DejaVu Sans}"
MONO_FONT="${FKCOMPUTE_DOCS_MONOFONT:-DejaVu Sans Mono}"

tmp_md="$(mktemp -t fkcompute_docs_XXXXXX.md)"
trap 'rm -f "${tmp_md}"' EXIT

{
  echo "---"
  echo "title: ${title}"
  echo "date: ${date_str}"
  echo "---"
  echo
} >"${tmp_md}"

preferred=(
  "README.md"
  "overview.md"
  "installation.md"
  "quickstart.md"
  "cli.md"
  "config.md"
  "python_api.md"
  "pipeline.md"
  "performance.md"
  "inversion.md"
  "constraint_system.md"
  "ilp.md"
  "cpp_backend.md"
  "output.md"
  "mathematica.md"
  "development.md"
  "troubleshooting.md"
)

declare -A seen

append_file() {
  local f="$1"
  if [[ ! -f "${f}" ]]; then
    return 0
  fi
  cat "${f}" >>"${tmp_md}"
  printf "\n\n\\newpage\n\n" >>"${tmp_md}"
}

# Add docs in preferred order first.
for base in "${preferred[@]}"; do
  seen["${base}"]=1
  append_file "${SCRIPT_DIR}/${base}"
done

# Append any remaining docs/*.md (alphabetical by basename).
mapfile -t all_basenames < <(cd "${SCRIPT_DIR}" && ls -1 *.md 2>/dev/null | LC_ALL=C sort)
for base in "${all_basenames[@]}"; do
  if [[ -n "${seen["${base}"]+x}" ]]; then
    continue
  fi
  append_file "${SCRIPT_DIR}/${base}"
done

mkdir -p "$(dirname "${OUT_PDF}")"

pandoc_font_args=()
if [[ "${PDF_ENGINE}" != "pdflatex" ]]; then
  pandoc_font_args+=(
    -V "mainfont=${MAIN_FONT}"
    -V "sansfont=${SANS_FONT}"
    -V "monofont=${MONO_FONT}"
  )
fi

pandoc \
  --from=markdown+raw_tex \
  --standalone \
  --toc \
  --toc-depth=2 \
  --number-sections \
  --highlight-style=tango \
  --resource-path="${REPO_ROOT}:${SCRIPT_DIR}" \
  --pdf-engine="${PDF_ENGINE}" \
  -H "${HEADER_TEX}" \
  "${pandoc_font_args[@]}" \
  -V geometry:margin=1in \
  -V colorlinks=true \
  -V linkcolor=blue \
  -V urlcolor=blue \
  -o "${OUT_PDF}" \
  "${tmp_md}"

echo "wrote: ${OUT_PDF}"
