"""fix_bbl.py -- two repairs to the .bbl that wileyNJD-Chicago.bst writes: it strips spurious \\natexlab suffixes, and it
fills the year of labels the style left empty (it does so for names in braces, such as {le Cessie}), reading the year
from the entry's own text. Refuses to write if the braces do not balance. Run between bibtex and the second xelatex."""
import re, sys

p = sys.argv[1] if len(sys.argv) > 1 else "localize_sim.bbl"
t = open(p, encoding="utf-8").read()
t, k1 = re.subn(r"\{\\natexlab\{[^{}]?\}\}|\{\\natexlab\{\{\}\}", "", t)
items = re.split(r"(?=\\bibitem\[)", t)
k2 = 0
for i, it in enumerate(items):
    m = re.match(r"\\bibitem\[(.*?)\]\{([^}]+)\}", it, re.S)
    if m and re.search(r"\}\{\}$", m.group(1)):
        body = it[m.end():]
        y = re.findall(r"\((\d{4})\)|(?:^|[\s,])((?:19|20)\d{2})[\.,]", body)
        yy = [a or b for a, b in y]
        if not yy:
            sys.exit("no year found for " + m.group(2))
        lab = m.group(1)[:-2] + "{" + yy[-1] + "}"
        items[i] = "\\bibitem[" + lab + "]{" + m.group(2) + "}" + body
        k2 += 1
        print("year", yy[-1], "->", m.group(2))
t = "".join(items)
if t.count("{") != t.count("}"):
    sys.exit("unbalanced braces; not written")
open(p, "w", encoding="utf-8", newline="\n").write(t)
print("natexlab removed", k1, "| years filled", k2)
