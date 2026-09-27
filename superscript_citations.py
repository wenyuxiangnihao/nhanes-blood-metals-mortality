# -*- coding: utf-8 -*-
"""Convert bracketed numeric citations ([12], [4-6,20,27]) in manuscript.docx
to superscript runs, matching Environmental Research's published citation style.
Reference-list entries (paragraphs starting with '<n>.') are left untouched.
Run AFTER md_to_docx.py, BEFORE PDF export."""
import re, sys, os
from docx import Document
from docx.oxml.ns import qn

CITE = re.compile(r"\[(\d+(?:\s*[,\u2013-]\s*\d+)*)\]")

def superserialize_paragraph(p):
    full = "".join(r.text for r in p.runs)
    if not CITE.search(full):
        return 0
    n = 0
    for r in p.runs:
        t = r.text
        if not CITE.search(t):
            continue
        parts = []
        pos = 0
        for m in CITE.finditer(t):
            if m.start() > pos:
                parts.append((t[pos:m.start()], False))
            parts.append((m.group(0), True))
            pos = m.end()
        if pos < len(t):
            parts.append((t[pos:], False))
        # write back: first part stays in this run, rest appended as new runs
        r.text = parts[0][0]
        if parts[0][1]:
            r.font.superscript = True
        prev_el = r._element
        for text, is_sup in parts[1:]:
            new = r._element.makeelement(qn("w:r"), {})
            rPr = r._element.makeelement(qn("w:rPr"), {})
            new.append(rPr)
            tEl = r._element.makeelement(qn("w:t"), {})
            tEl.text = text
            tEl.set(qn("xml:space"), "preserve")
            new.append(tEl)
            prev_el.addnext(new)
            prev_el = new
            # set superscript via python-docx wrapper
            from docx.text.run import Run
            Run(new, p).font.superscript = is_sup
            n += 1
    return n

def main(path):
    doc = Document(path)
    total = 0
    ref_mode = False
    for p in doc.paragraphs:
        text = p.text.strip()
        if text == "References":
            ref_mode = True
            continue
        if ref_mode and re.match(r"^\d+\.\s", text):
            continue  # reference list entry
        total += superserialize_paragraph(p)
    doc.save(path)
    print("superscripted citation runs:", total)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join("submission", "manuscript.docx"))
