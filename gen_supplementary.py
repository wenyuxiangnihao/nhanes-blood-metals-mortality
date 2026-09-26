#!/usr/bin/env python3
import os
ROOT = os.environ.get("NHANES_ROOT") or os.path.dirname(os.path.abspath(__file__))
# -*- coding: utf-8 -*-
"""Build submission/supplementary_tables.md (Table S1-S20) and its DOCX,
plus the STROBE checklist DOCX, reusing the converter in md_to_docx.py."""
import io, os, sys

BASE = ROOT
sys.path.insert(0, os.path.join(BASE, 'scripts'))
import md_to_docx as M

tables = io.open(os.path.join(BASE, 'tables.md'), encoding='utf-8').read()
i = tables.index('**Table S')
supp = "# Supplementary tables\n\n" + tables[i:].rstrip() + "\n"
p_md = os.path.join(BASE, 'submission', 'supplementary_tables.md')
os.makedirs(os.path.dirname(p_md), exist_ok=True)
io.open(p_md, 'w', encoding='utf-8').write(supp)
print('wrote', p_md, len(supp), 'chars;', supp.count('**Table S'), 'supplementary tables')

M.build(p_md, os.path.join(BASE, 'submission', 'supplementary_tables.docx'))
M.build(os.path.join(BASE, 'STROBE_checklist.md'),
        os.path.join(BASE, 'submission', 'STROBE_checklist.docx'),
        title_center=False, double=False)
