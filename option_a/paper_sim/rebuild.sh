#!/bin/sh
# rebuild.sh -- supplement first (the paper reads its labels), then the paper
set -e
python build_sim_tex.py </dev/null
rm -f localize_sim.bbl localize_sim.aux
xelatex -interaction=nonstopmode -halt-on-error localize_sim_supplement.tex </dev/null >/dev/null
xelatex -interaction=nonstopmode -halt-on-error localize_sim_supplement.tex </dev/null >/dev/null
xelatex -interaction=nonstopmode -halt-on-error localize_sim.tex </dev/null >/dev/null
bibtex localize_sim </dev/null >/dev/null || true
python fix_bbl.py localize_sim.bbl </dev/null
xelatex -interaction=nonstopmode -halt-on-error localize_sim.tex </dev/null >/dev/null
xelatex -interaction=nonstopmode -halt-on-error localize_sim.tex </dev/null >/dev/null
for f in localize_sim localize_sim_supplement; do
  pdfinfo $f.pdf | grep Pages
  grep -c "undefined" $f.log || true
  grep -E "^!|Citation .* undefined|Reference .* undefined" $f.log | head
done
