#!/usr/bin/env python
# vim: ai

import sql
import sql.aggregate
import logging
import sys

sys.path.append('../pylib')

import glob
import os.path
import subprocess
import itertools
import re

from historytable import Database, HistoryTable, Table

import hashlib

logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(level=logging.INFO)

DB = Database()

slc_like = HistoryTable('slc_like', DB)
slc_like_clusters = Table('slc_like_clusters', DB)
uniprot_proteins = HistoryTable('uniprot_proteins', DB)

q = slc_like.join(
    uniprot_proteins, type_='LEFT',
    condition=(slc_like.accession == uniprot_proteins.accession)
).join(
    slc_like_clusters, type_='LEFT',
    condition=(slc_like.accession == slc_like_clusters.accession)
).select(
    slc_like.accession,
    slc_like_clusters.branch_id,
    slc_like_clusters.family_name,
    uniprot_proteins.tax_id,
    uniprot_proteins.seq_fasta,
    where=(slc_like.status >= 0)
)

cur = DB.query(q, dictionary=True)
all_accessions = cur.fetchall()
cur.close()

#naturalsort = lambda a: tuple((str, int)[x%2](y) for x, y in enumerate(re.split(r'(\d+)', a)))
def naturalsort(a):
    try:
        return tuple((str, int)[x%2](y) for x, y in enumerate(re.split(r'(\d+)', a)))
    except TypeError as e:
        logging.error(f'TypeError with argument: {a}')
        raise e

all_accessions.sort(key=lambda a: naturalsort(a['family_name']))

# rename in tree
tax2tag = {
        9606: 'Homo',
        10090: 'Mus',
        10116: 'Rattus',
        7227: 'Drosophila',
        9031: 'Gallus',
        6239: 'Caenorhabditis',
        7955: 'Danio'
}

# SVG templates

svg_header = """\
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->

<svg
   width="210mm"
   height="297mm"
   viewBox="0 0 210 297"
   version="1.1"
   id="svg11165"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:svg="http://www.w3.org/2000/svg">
  <defs
     id="defs11162" />
  <g
     id="layer1"
     transform="translate({transform_x:.6f}, {transform_y:.6f}) scale({scale:.6f})">

"""

svg_rect = """\
    <rect
       style="fill:{fill};fill-rule:evenodd;stroke:{stroke};stroke-width:0.264583"
       width="{width:0.6f}"
       height="{height:0.6f}"
       x="{x:0.6f}"
       y="{y:0.6f}" />
"""

svg_text = """\
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:{font_weight};font-size:{font_size:0.6f}px;line-height:1.25;font-family:sans-serif;fill:{fill};fill-opacity:1;stroke:none;stroke-width:0.264583"
       x="{x:0.6f}"
       y="{y:0.6f}"
       text-anchor="{text_anchor}"
       ><tspan
         style="stroke-width:0.264583"
         x="{x:0.6f}"
         y="{y:0.6f}">{text}</tspan></text>
"""

svg_text_transformed = """\
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:{font_weight};font-size:{font_size:0.6f}px;line-height:1.25;font-family:sans-serif;fill:{fill};fill-opacity:1;stroke:none;stroke-width:0.264583"
       x="{x:0.6f}"
       y="{y:0.6f}"
       text-anchor="{text_anchor}"
       transform="{transform}"
       ><tspan
         style="stroke-width:0.264583"
         x="{x:0.6f}"
         y="{y:0.6f}">{text}</tspan></text>
"""

svg_footer = """\
  </g>
</svg>
"""

#####

# branch_id IS NOT compatible with family_name
# one family can have multiple branches, sometimes it's not possible to merge
# all branches, e.g. SLC38 or SLC51 families.
pdf_fns = []
svg_lines = []
tab_x = [
    0.0,
] + [60.0 + x*10.0 for x in range(6)]
line_y = 0.0
line_height = 6.0
species_order = [
        'Homo',
        'Mus',
        'Rattus',
        'Danio',
        'Gallus',
        'Drosophila',
        'Caenorhabditis',
]
col_x = 0.0
col_width = 130.0
line_counter = 0
lines_per_col = 152

for family_name, members in itertools.groupby(all_accessions,
                                              key=lambda a: a['family_name']):
    members = sorted(members, key=lambda a: a['accession'])
    # generate key unique to family members
    key = hashlib.sha1(
        ('\n'.join([x['accession'] for x in members])).encode('utf-8')
    ).hexdigest()
    logging.info('{}: ({}) {:d} members'.format(
        members[0]['family_name'],
        ', '.join(['#{:d}'.format(x) for x in set([y['branch_id'] for y in
                                                   members])]),
        len(members)
    ))
    logging.debug('members: '+' '.join([x['accession'] for x in members]))
    rows = []
    if len(members) < 3:
        # single-member cluster, no alignment, no tree, no svg
        # OR
        # two proteins in cluster, no tree, no svg
        orthologs = dict()
        for m in members:
            taxon = tax2tag[m['tax_id']]
            orthologs[taxon] = orthologs.get(taxon, 0) + 1
        rows.append(orthologs)
    else:
        fn = f'branch_{key}/make-tree-svg.err'
        with open(fn, 'rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('INFO:root:Homo'):
                    species = line[10:].strip().split()
                    #print(species)
                    for line in f:
                        if line.startswith('INFO:root:'):
                            orthologs = dict()
                            row = line[10:].split()
                            if len(row) < 7: row.insert(0, '')
                            for s, r in zip(species, row):
                                orthologs[s] = r
                            #print(len(row), row)
                            rows.append(orthologs)
    svg_lines.append(
        svg_text.format(x=tab_x[0] + col_x,
                        y=line_y + line_height*1.9,
                        text=family_name, font_size=7.0, fill='#000000',
                        font_weight='bold',
                        text_anchor='start')
    )
    line_y += line_height * 2.3
    line_counter += 1
    for orthologs in rows:
        if orthologs.get(species_order[0], '') != '':
            label = orthologs[species_order[0]]
            label = '{} ({})'.format(*(label.split('_')[:2]))
            svg_lines.append(
                svg_text.format(x=tab_x[0] + col_x,
                                y=line_y + line_height*0.8,
                                text=label, font_size=5.0, fill='#000000',
                                font_weight='normal',
                                text_anchor='start')
            )

        for i, taxon in enumerate(species_order[1:]):
            label = orthologs.get(taxon, 0)
            if isinstance(label, int) or (isinstance(label, str) and
                                          label.isdigit()):
                label = int(label)
                if label == 1:
                    fill = '#808080'
                elif label >= 2:
                    fill = '#404040'
                else:
                    fill = '#ffffff'
            else:
                if label == '': continue
                elif label == 'M':
                    fill = '#c0c0c0'
                else:
                    fill = '#ffffff'
            svg_lines.append(
                svg_rect.format(x=tab_x[i+1] + col_x,
                                y=line_y,
                                width=10.0,
                                height=line_height,
                                fill=fill,
                                stroke='#000000')
            )
        line_y += line_height
        line_counter += 1
        #if line_counter > lines_per_col:
        if (line_y / line_height) > lines_per_col:
            line_y = 0.0
            line_counter = 0
            col_x += col_width

# labels

x_shift = 0.0
while x_shift <= col_x:
    x = tab_x[0]+3.0 + x_shift
    svg_lines.append(
        svg_text_transformed.format(
            x=x,
            y=-5.0,
            text='human', font_size=5.0, fill='#000000',
            font_weight='normal',
            text_anchor='start',
            transform='rotate(-45 {:.6f} {:.6f})'.format(x, -5.0)
        )
    )
    for i, taxon in enumerate((
        'mouse',
        'rat',
        'zebrafish',
        'chicken',
        'fruitfly',
        'roundworm')):
        x = tab_x[i+1]+3.0 + x_shift
        svg_lines.append(
            svg_text_transformed.format(
                x=x,
                y=-5.0,
                text=taxon, font_size=5.0, fill='#000000',
                font_weight='normal',
                text_anchor='start',
                transform='rotate(-45 {:.6f} {:.6f})'.format(x, -5.0)
            )
        )
    x_shift += col_width

print(svg_header.format(
    transform_x=10.0,
    transform_y=10.0,
    scale=0.25,
))

print(''.join(svg_lines))

print(svg_footer)

