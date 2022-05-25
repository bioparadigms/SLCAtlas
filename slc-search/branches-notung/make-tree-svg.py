#!/usr/bin/env python

from parse_newick import parse_newick
from pprint import pprint
import logging
import re
import sys

logging.basicConfig(level=logging.INFO)

def parse_NHX_features(feature_string):
    '''
    Example: 
    '&&NHX:conf=0.01:name=INTERNAL' -> {'conf':'0.01', 'name': 'INTERNAL'}
    '''
    prefix = '&&NHX:'
    feature_string = feature_string[len(prefix):]
    items = feature_string.split(':')
    key_value_pairs = (item.split('=') for item in items)
    return dict(key_value_pairs)

def tree_builder(label, children, distance, features):
    return {'label': label, 'children': children, 'distance': distance, 'features': features}
    
# branch_*.tree.rooting.*.rearrange.*.reconciled
tree_fn = sys.argv[1]
f = open(tree_fn, 'rt', encoding='utf-8')
buf = f.read()
lines = buf.split('\n')

#print(lines[:-3])

t = parse_newick('\n'.join(lines[:-3]),
                 aggregator=tree_builder,
                 feature_parser=parse_NHX_features)
f.close()

# species tree
prefix = '[&&NOTUNG-SPECIES-TREE'
for line in lines:
    if line.startswith(prefix):
        species_tree_str = line[len(prefix):].strip().rstrip(']') + ';'
species_t = parse_newick(species_tree_str, aggregator=tree_builder)

ancestor_species = []

all_species = []
p = [species_t]
r = 0
# a list for walking the tree from the root downwards, or from the leaves up
while r < len(p):
    if len(p[r]['children']) > 0:
        p.extend(p[r]['children'])
    r += 1
# walk from the leaves up
r = len(p)-1
while r >= 0:
    e = p[r]
    if len(e['children']) > 0:
        e['human_anc'] = any(x['human_anc'] for x in e['children'])
    else:
        e['human_anc'] = (e['label'] == 'Homo')
        all_species.append(e['label'])
    if e['human_anc']: ancestor_species.append(e['label'])
    r -= 1

logging.info(repr(ancestor_species))
ancestor_species = frozenset(ancestor_species)

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
  <text
     xml:space="preserve"
     style="font-style:normal;font-weight:normal;font-size:7.05556px;line-height:1.25;font-family:sans-serif;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.264583"
     x="10"
     y="15"
     id="text26538">{title}</text>
  <path
     style="fill:none;stroke:#000000;stroke-width:0.5;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;stroke-miterlimit:4;stroke-dasharray:none"
     d="M 10,22 H 200"
     id="path44767" />
  {warning}
  {species_tree}
  <g
     id="layer1"
     transform="translate({transform_x:.6f}, {transform_y:.6f}) scale({scale:.6f})">

"""

svg_warning = """\
  <text
     xml:space="preserve"
     style="font-style:normal;font-weight:normal;font-size:3px;line-height:1.25;font-family:sans-serif;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.264583"
     x="10"
     y="28"
     id="text999">
        <tspan x="10" y="28">Warning: this is a large figure that had to be reduced to fit on the page.</tspan>
        <tspan x="10" y="32">Please use the zoom function of your PDF viewer to see the details.</tspan>
        </text>
"""


svg_path = """\
    <path
       style="fill:none;stroke:#000000;stroke-width:0.264583px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       d="M {x1:0.6f},{y1:0.6f} {x2:0.6f},{y2:0.6f}"
        />
"""
svg_branch_path = """\
    <path
       style="fill:none;stroke:#000000;stroke-width:0.264583px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       d="m {x0:0.6f},{y0:0.6f} {h1:0.6f},0.0 0.0,{width:0.6f} {h2:0.6f},0.0"
        />
"""

svg_rect = """\
    <rect
       style="fill:{fill};fill-rule:evenodd;stroke:{stroke};stroke-width:0.264583"
       width="{width:0.6f}"
       height="{width:0.6f}"
       x="{x:0.6f}"
       y="{y:0.6f}" />
"""

svg_text = """\
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:{font_size:0.6f}px;line-height:1.25;font-family:sans-serif;fill:{fill};fill-opacity:1;stroke:none;stroke-width:0.264583"
       x="{x:0.6f}"
       y="{y:0.6f}"
       text-anchor="{text_anchor}"
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

x_scaling = 3.0
y_scaling = 3.0

label_gap = 1.0
duplication_label_gap = 0.3
bootstrap_gap_x = 0.2
bootstrap_gap_y = 0.2
baseline_shift = 0.25


# svg
def branch(x0, y0, e):
    # line from 
	return svg_branch_path.format(x0=x0-e['children'][0]['x'], y0=y0+e['children'][0]['y'], h1=-e['children'][0]['height'], h2=e['children'][1]['height'], width=e['children'][1]['y']-e['children'][0]['y'])

# draw species tree
species_t_svg = []

# SPECIES TREE calculate y (leaves order)
species_t_leaves = [species_t]
r = 0
while r < len(species_t_leaves):
    if len(species_t_leaves[r]['children']) > 0:
        species_t_leaves = species_t_leaves[:r] + species_t_leaves[r]['children'] + species_t_leaves[r+1:]
    else:
        r += 1
for y, e in enumerate(species_t_leaves):
    e['y'] = float(y) * 3.0
    e['y_min'] = float(y) * 3.0
    e['y_max'] = float(y) * 3.0

# SPECIES TREE count n
p = [species_t]
r = 0
# a list for walking the tree from the root downwards, or from the leaves up
while r < len(p):
    if len(p[r]['children']) > 0:
        p.extend(p[r]['children'])
    r += 1
# walk from the leaves up
r = len(p)-1
while r >= 0:
    e = p[r]
    if len(e['children']) > 0:
        e['n'] = sum(x['n'] for x in e['children'])
        # this is for 2 children
        e['children'][0]['height'] = e['children'][1]['n'] * x_scaling
        e['children'][1]['height'] = e['children'][0]['n'] * x_scaling
        e['x'] = e['children'][0]['x'] + e['children'][0]['height']
        e['y_min'] = min([x['y_min'] for x in e['children']])
        e['y_max'] = max([x['y_max'] for x in e['children']])
        e['y'] = (e['y_min'] + e['y_max']) / 2.0
        # SVG
        species_t_svg.append(branch(0.0, 0.0, e))
        # x0-(x of child) - (height of child)
        x = 0.0-e['children'][0]['x'] - e['children'][0]['height']
    else:
        # is a leaf
        e['n'] = 1
        e['x'] = 0
    r -= 1
#pprint(t)

species_t_svg.append(
    svg_path.format(x1=0.0-e['x'], y1=0.0+e['y'],
                    x2=0.0-e['x']-x_scaling, y2=0.0+e['y'])
)

for e in species_t_leaves:
    species_t_svg.append(
        svg_text.format(x=0.0-e['x']+label_gap*x_scaling,
                        y=0.0+e['y']+baseline_shift*y_scaling,
                        #text=e['label']+ ' ({})'.format(e['partition']),
                        text=e['label'],
                        font_size=2.5, 
                        fill='#000000' if 'LOST' not in e['label'] else '#a0a0a0',
                        text_anchor='start')
    )

# walk from the root down
r = 0
while r < len(p):
    e = p[r]
    if len(e['children']) > 0:
        species_t_svg.append(
            svg_text.format(x=0.0-e['x']-bootstrap_gap_x*x_scaling,
                            y=0.0+e['y']-bootstrap_gap_y*y_scaling,
                            text=e['label'], font_size=1.5,
                            fill='#000000',
                            text_anchor='end')
        )
    r += 1






# calculate y (leaves order)
leaves = [t]
r = 0
while r < len(leaves):
    if len(leaves[r]['children']) > 0:
        leaves = leaves[:r] + leaves[r]['children'] + leaves[r+1:]
    else:
        r += 1
for y, e in enumerate(leaves):
    e['y'] = float(y) * y_scaling
    e['y_min'] = float(y) * y_scaling
    e['y_max'] = float(y) * y_scaling

#print([x['label'] for x in leaves])

svg_lines = []

min_x = 1e6

# count n
p = [t]
r = 0
# a list for walking the tree from the root downwards, or from the leaves up
while r < len(p):
    if len(p[r]['children']) > 0:
        p.extend(p[r]['children'])
    r += 1
# walk from the leaves up
r = len(p)-1
while r >= 0:
    e = p[r]
    if len(e['children']) > 0:
        e['n'] = sum(x['n'] for x in e['children'])
        # this is for 2 children
        e['children'][0]['height'] = e['children'][1]['n'] * x_scaling
        e['children'][1]['height'] = e['children'][0]['n'] * x_scaling
        e['x'] = e['children'][0]['x'] + e['children'][0]['height']
        e['y_min'] = min([x['y_min'] for x in e['children']])
        e['y_max'] = max([x['y_max'] for x in e['children']])
        e['y'] = (e['y_min'] + e['y_max']) / 2.0
        # SVG
        svg_lines.append(branch(0.0, 0.0, e))
        # x0-(x of child) - (height of child)
        x = 0.0-e['children'][0]['x'] - e['children'][0]['height']
        if x < min_x: min_x = x
        # count human proteins
        e['n_human'] = sum([x['n_human'] for x in e['children']])
    else:
        # is a leaf
        e['n'] = 1
        e['x'] = 0
        # count human proteins
        e['n_human'] = 1 if e['label'].endswith('_Homo') else 0
    r -= 1
#pprint(t)

svg_lines.append(
    svg_path.format(x1=0.0-e['x'], y1=0.0+e['y'],
                    x2=0.0-e['x']-x_scaling, y2=0.0+e['y'])
)

# walk from the root down
r = 0
while r < len(p):
    e = p[r]
    if ('D' in e['features']) and (e['features']['D'] == 'Y'):
        # draw "duplication" node
        svg_lines.append(
            svg_rect.format(x=0.0-e['x']-0.1*x_scaling,
                            y=0.0+e['y']-0.1*y_scaling,
                            width=0.2*x_scaling,
                            height=0.2*y_scaling,
                            fill='#ff0000',
                            stroke='#ff0000')
        )
        svg_lines.append(
            svg_text.format(x=0.0-e['x']+duplication_label_gap*x_scaling,
                            y=0.0+e['y']+baseline_shift*y_scaling,
                            text='D', font_size=2.0, fill='#ff0000',
                            text_anchor='start')
        )
    if 'B' in e['features']:
        svg_lines.append(
            svg_text.format(x=0.0-e['x']-bootstrap_gap_x*x_scaling,
                            y=0.0+e['y']-bootstrap_gap_y*y_scaling,
                            text=e['features']['B'], font_size=1.5,
                            fill='#00c000',
                            text_anchor='end')
        )
    r += 1

# partition tree
# walk from the root down
r = 0
partitions = ['']
p[0]['partition'] = ''
while r < len(p):
    e = p[r]
    if ('D' in e['features']) and (e['features']['D'] == 'Y'):
        # is duplication node
        if e['features']['S'] in ancestor_species:
            # new partition
            partitions.remove(e['partition'])
            for i, ee in enumerate(e['children']):
                ee['partition'] = e['partition'] + str(i)
                partitions.append(ee['partition'])
        else:
            # no new partition because duplication happened after human
            for i, ee in enumerate(e['children']):
                ee['partition'] = e['partition']
    else:
        # not duplication note
        for i, ee in enumerate(e['children']):
            ee['partition'] = e['partition']
    r += 1
logging.info(repr(partitions))

partition_genes = dict([(x, []) for x in partitions])
multiple = set()
for e in leaves:
    svg_lines.append(
        svg_text.format(x=0.0-e['x']+label_gap*x_scaling,
                        y=0.0+e['y']+baseline_shift*y_scaling,
                        #text=e['label']+ ' ({})'.format(e['partition']),
                        text=e['label'],
                        font_size=2.5, 
                        fill='#000000' if 'LOST' not in e['label'] else '#a0a0a0',
                        text_anchor='start')
    )
    for pp in partitions:
        if pp.startswith(e['partition']):
            partition_genes[pp].append(e['label'])
            if e['partition'] not in partitions:
                multiple.add(e['label'])

for pp, genes in partition_genes.items():
    logging.info(pp)
    logging.info(repr(genes))
logging.info(repr(multiple))

max_y = (len(leaves) + 1) * y_scaling

# scale tree if necessary
transform_x = 135
transform_y = 35
scaling = 1.0
# 290 = 297 - 7 margin
if (max_y + transform_y) > 290:
    scaling = (290-transform_y) / max_y
# 10 = margin
if (min_x*scaling + transform_x) < 10:
    scaling = (10-transform_x) / min_x
logging.info('scaling by {:.6f}'.format(scaling))
warning_message = ''
if scaling < 0.75:
    # scaling too small, we try to make use of more space
    # we reserve 50 for labels
    width = 50 - min_x
    # 20 = margin_left + margin_right
    scaling = (210 - 20) / width
    # 10 = margin_left
    transform_x = -min_x*scaling + 10
    if (max_y*scaling + transform_y) > 290:
        scaling = (290-transform_y) / max_y
        transform_x = -min_x*scaling + 10

if scaling < 0.75:
    warning_message = svg_warning

partition_table = []
#spec_list = ['Homo'] + [x for x in sorted(all_species) if x != 'Homo']
spec_list = ['Homo', 'Rattus', 'Mus', 'Gallus', 'Danio', 'Drosophila',
             'Caenorhabditis']
for pp, genes in partition_genes.items():
    human_gene = [x for x in genes if x.endswith('_Homo')]
    human_label = human_gene[0] if len(human_gene) > 0 else ''
    per_species = dict()
    for label in genes:
        s = label.split('_')[-1]
        if label in multiple:
            # "M" means that this gene (partition) shares an ortholog with
            # other genes in this organism
            per_species[s] = 'M'
        else:
            # a number means this gene has that many orthologs in this organism
            per_species[s] = per_species.get(s, 0) + 1
    line = [human_label]
    for s in spec_list[1:]:
        line.append(str(per_species.get(s, '0')))
    partition_table.append(line)

naturalsort = lambda a: tuple((str, int)[x%2](y) for x, y in
                              enumerate(re.split(r'(\d+)', a)))
partition_table.sort(key=lambda a: naturalsort(a[0]))
template = '{:<30s}' + '{:>15s}'*(len(spec_list)-1)
logging.info(template.format(*spec_list))
for line in partition_table:
    template = '{:<30s}' + '{:>15s}'*(len(line)-1)
    logging.info(template.format(*line))

# SVG generation

species_tree_svg = """\
  <g
     id="species-tree"
     transform="translate(190, 7) scale(0.6)">
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:3px;line-height:1.25;font-family:sans-serif;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.264583"
       x="-25"
       y="-4"
       >Species tree:</text>
  {}
  </g>
""".format(''.join(species_t_svg))

print(svg_header.format(title='@@@TITLE@@@',
                        scale=scaling,
                        transform_x=transform_x,
                        transform_y=transform_y,
                        warning=warning_message,
                        species_tree=species_tree_svg))

print(''.join(svg_lines))

print(svg_footer)


