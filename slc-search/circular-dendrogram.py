#!/usr/bin/env python
# vim: ai

import glob
import sys
import gzip
import logging

import numpy as np
import math
import seaborn as sns
import itertools

import sql
import sql.aggregate

sys.path.append('./pylib')

from historytable import Database, HistoryTable, Table

user = 'circular-dendrogram.py'

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

DB = Database()

sl = Table('pdb_hits', DB)
slc_like = HistoryTable('slc_like', DB)
slc_like_clusters = Table('slc_like_clusters', DB)
slc_like_clusters_linkage = Table('slc_like_clusters_linkage', DB)
uniprot_proteins = HistoryTable('uniprot_proteins', DB)


# human_proteins
# Z
# slc_to_br -> slc_like_clusters accession -> branch_id
# br_to_slc_family -> slc_like_clusters branch_id -> family_name
# branch_of_h
# read_only
read_only = False

# need to prune non-human nodes
# Z matrix:
##### a .. ..          c|              b|
##### b d e           +---+             |
##### c a b          a|  b|      ->     |  
#####                   +---+         +---+
#####                  d|  e|        d|  e|
##### c d e
##### len(c) := len(c) + len(b)
##### num(c) := num(b)
##### --- OR ---
#
# a .. ..            d|              d|
# c a b             +---+           +---+
# d c f            c|  f|          b|  f|
#                 +---+             |
# d b f          a|  b|      ->     |

# we only cluster slc_like where status >= 0
# so after including/excluding proteins, proteins should be re-clustered.
# has to be in the same order as in cluster-slc-like.py
q = slc_like.join(
    slc_like_clusters, type_='LEFT',
    condition=(slc_like_clusters.accession == slc_like.accession)
).join(
    uniprot_proteins, type_='LEFT',
    condition=(uniprot_proteins.accession == slc_like.accession)
).select(
    slc_like.accession,
    uniprot_proteins.symbol,
    uniprot_proteins.tax_id,
    slc_like_clusters.branch_id,
    slc_like_clusters.family_name,
    where=(slc_like.status >= 0),
    order_by=slc_like.record_id.asc
)
cur = DB.query(q)
all_slc_like = []
slc_to_br = dict()
br_to_slc_family = dict()
for row in cur.fetchall():
    accession, gene_symbol, tax_id, branch_id, family_name = row
    slc_to_br[accession] = branch_id
    if branch_id not in br_to_slc_family:
        br_to_slc_family[branch_id] = family_name
    elif br_to_slc_family[branch_id] != family_name:
        logging.warning('multiple definitions for family_name for branch_id '
                        '{:d}'.format(branch_id))
    all_slc_like.append((accession, gene_symbol, tax_id))
cur.close()

q = slc_like_clusters_linkage.select(
    slc_like_clusters_linkage.cluster_id1,
    slc_like_clusters_linkage.cluster_id2,
    slc_like_clusters_linkage.distance,
    slc_like_clusters_linkage.n_elements,
)
cur = DB.query(q)
Z = np.array(cur.fetchall())
cur.close()
print(Z)
print('len(Z) =', len(Z))
print('Z[-1][-1] =', Z[-1][-1])
print('len(all_slc_like) =', len(all_slc_like))
# len(Z)+1 = Z[-1][-1] = len(all_slc_like)

# if we delete a branch that defines a family, we need to provide a replacement
# branch id
# we record these here... if deleting branch "a" then branch "b" will be its
# replacement always (see tree diagram above)
replacement_branches = dict()
def prune(Z, branch_id):
    # property of the linkage matrix
    logging.debug(f'trying to prune branch {branch_id}')
    n_leaves = len(Z)+1
    # we put _id as suffix to indicate it is a branch ID and not index into Z
    a_id = branch_id
    b_id = None
    for c in range(len(Z)):
        br1, br2, l, n = Z[c, :]
        if br1 == a_id: b_id = br2
        if br2 == a_id: b_id = br1
        if b_id is not None:
            logging.debug('found branch: {:d} = {!r}'.format(c+n_leaves, Z[c, :]))
            replacement_branches[c+n_leaves] = int(b_id)
            for d in range(c+1, len(Z)):
                br1, br2, super_l, super_n = Z[d, :]
                if br1 == (c+n_leaves):
                    logging.debug('found branch above: '
                                  '{:d} = {!r}'.format(d+n_leaves, Z[d, :]))
                    Z[d, 0] = b_id
                    if a_id < n_leaves: Z[d, 3] -= 1
                    else: Z[d, 3] -= Z[a_id-n_leaves, 3]
                    logging.debug('new branch: '
                                  '{:d} = {!r}'.format(d+n_leaves, Z[d, :]))
                    break
                if br2 == (c+n_leaves):
                    logging.debug('found branch above: '
                                  '{:d} = {!r}'.format(d+n_leaves, Z[d, :]))
                    Z[d, 1] = b_id
                    if a_id < n_leaves: Z[d, 3] -= 1
                    else: Z[d, 3] -= Z[a_id-n_leaves, 3]
                    logging.debug('new branch: '
                                  '{:d} = {!r}'.format(d+n_leaves, Z[d, :]))
                    break
            Z[c, 0] = -1
            Z[c, 1] = -1
            break

# prune non-human leaves
human_proteins = []
human_proteins_names = []
human_proteins_idx = []
branch_lookup = dict()
for i, (acc, gene_symbol, tax_id) in enumerate(all_slc_like):
    if tax_id == 9606:
        human_proteins.append(acc)
        if gene_symbol is None:
            human_proteins_names.append(acc)
        else:
            human_proteins_names.append(gene_symbol)
        branch_lookup[i] = len(human_proteins_idx)
        human_proteins_idx.append(i)
    else:
        prune(Z, i)

# compact the Z-matrix: remove unused branch_ids and renumber remaining ones

new_Z = []
for i, l in enumerate(Z):
    if l[0] < 0:
        try:
            assert(i+len(Z)+1 in replacement_branches)
        except AssertionError as e:
            print('i =', i)
            print('l =', l)
            raise(e)
        continue
    if l[1] < 0:
        try:
            assert(i+len(Z)+1 in replacement_branches)
        except AssertionError as e:
            print('i =', i)
            print('l =', l)
            raise(e)
        continue
    branch_lookup[i+len(all_slc_like)] = len(new_Z)+len(human_proteins)
    new_Z.append([float(branch_lookup[int(l[0])]),
                  float(branch_lookup[int(l[1])]),
                  l[2],
                  l[3]])
print('len(new_Z) =', len(new_Z))
print('len(branch_lookup) =', len(branch_lookup))

# adjust family br_ids, as well, because maybe those were deleted above
for acc, br_id in slc_to_br.items():
    while br_id in replacement_branches:
        br_id = replacement_branches[br_id]
    try:
        # the -1 is generated when br_id is a deleted leaf
        br_id = branch_lookup.get(br_id, -1)
    except KeyError as e:
        print('br_id =', br_id)
        print('Z[br_id, :] =', Z[br_id-len(all_slc_like), :])
        raise(e)
    slc_to_br[acc] = br_id

new_br_to_slc_family = dict()
for br_id, family_name in list(br_to_slc_family.items()):
    while br_id in replacement_branches:
        br_id = replacement_branches[br_id]
    try:
        # the -1 is generated when br_id is a deleted leaf
        br_id = branch_lookup.get(br_id, -1)
    except KeyError as e:
        print('br_id =', br_id)
        print('Z[br_id, :] =', Z[br_id-len(all_slc_like), :])
        raise(e)
    if br_id < 0: continue
    new_br_to_slc_family[br_id] = family_name
br_to_slc_family = new_br_to_slc_family

print('***** slc_to_br:')
print(slc_to_br)
print('***** br_to_slc_family:')
print(br_to_slc_family)

# finalize new Z-matrix
Z = np.array(new_Z)


svg_header_base = """<?xml version="1.0" encoding="utf-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="{0:d}pt" height="{1:d}pt" viewBox="0 0 {0:d} {1:d}" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
 <defs>
  <style type="text/css">
*{{stroke-linecap:butt;stroke-linejoin:round;stroke-miterlimit:100000;}}
  </style>
 </defs>
 <g id="figure_1">
"""

svg_footer = """ </g>
</svg>
"""

def p2c(cx, cy, r, phi):
    phi_rad = (phi-90)*math.pi/180.0
    return (cx + r*math.cos(phi_rad), cy + r*math.sin(phi_rad))

def arc(x, y, r, phi0, phi1):
    return '<path d="M {0[0]:f} {0[1]:f} A {2:f} {2:f} 0 {3:d} 1 {1[0]:f} {1[1]:f}" style="fill:none;stroke:#000000;"/>'.format(
		p2c(x, y, r, phi0), p2c(x, y, r, phi1), r, (phi1-phi0) > 180)

def sector(cx, cy, r0, phi0, r1, phi1, color='#000000', opacity=1.0):
    # cx, cy = centre of circle
    # r0, phi0 = top-left in polar coordinates
    # r1, phi1 = bottom-right in polar coordinates
    return '''<path d="M {p0[0]:f} {p0[1]:f}
    L {p1[0]:f} {p1[1]:f} 
    A {irad:f} {irad:f} 0 {large:d} 1 {p2[0]:f} {p2[1]:f} 
    L {p3[0]:f} {p3[1]:f}
    A {orad:f} {orad:f} 0 {large:d} 0 {p0[0]:f} {p0[1]:f}" style="fill: {color}; stroke: none; fill-opacity: {opacity:f}"/>
    '''.format(
        p0=p2c(cx, cy, r0, phi0),
        p1=p2c(cx, cy, r1, phi0),
        irad=r1,
        large=(phi1-phi0) > 180,
        p2=p2c(cx, cy, r1, phi1),
        p3=p2c(cx, cy, r0, phi1),
        orad=r0,
        color=color,
        opacity=opacity)

def branch(cx, cy, ix, dx, color='#000000'):
    # cx, cy = centre of circle
    # ix = 4x float, element of icoords = phi
    # dx = 4x float, corresponding element of dcoords = radius
    return '''<path d="M {start[0]:f} {start[1]:f}
    L {p1[0]:f} {p1[1]:f} 
    A {irad:f} {irad:f} 0 {large:d} 1 {p2[0]:f} {p2[1]:f} 
    L {end[0]:f} {end[1]:f}" style="fill: none; stroke: {color};"/>
    '''.format(
        start=p2c(cx, cy, dx[0], ix[0]),
        p1=p2c(cx, cy, dx[1], ix[1]),
        irad=dx[1],
        large=(ix[2]-ix[1]) > 180,
        p2=p2c(cx, cy, dx[2], ix[2]),
        end=p2c(cx, cy, dx[3], ix[3]),
        color=color)

def label(c, phi, text, color='#000000', font_size='4pt', dy='0.35'):
    # c = centre (x, y)
    return '''<text style="font-family: Bitstream Vera Sans; font-size: {font_size}; font-style: normal; fill: {color};
    text-anchor: {anchor};"
    transform="rotate({angle:f}, {start[0]:f}, {start[1]:f}) translate(0 {dy:f})" x="{start[0]:f}" y="{start[1]:f}">{text:s}</text>\n'''.format(
        start=c, angle=phi-90.0 if phi < 180.0 else 90.0+phi, text=text, anchor=("start", "end")[int(phi >= 180.0)],
        color=color, font_size=font_size, dy=dy)

color_codes = 'bgrcmykw'
color_values = ['#000000'] + \
['#{:02x}{:02x}{:02x}'.format(int(x[0]*255), int(x[1]*255), int(x[2]*255)) for x in sns.color_palette('colorblind')]

color_map = dict(zip(color_codes, color_values[:len(color_codes)]))
print(color_map)

# we expand here the whole tree to come up with the order the proteins should
# be presented on the X axis
leaves = [len(human_proteins)+len(Z)-1]
# the last cluster uniting everything
#leaves = [len(all_slc_like)+len(Z)-1]
i = 0
while i < len(leaves):
    if leaves[i] < len(human_proteins): i += 1
    #if leaves[i] < len(all_slc_like): i += 1
    else:
        #leaves = leaves[:i] + list(map(int, Z[leaves[i]-len(all_slc_like)][:2])) + leaves[i+1:]
        leaves = leaves[:i] + list(map(int, Z[leaves[i]-len(human_proteins)][:2])) + leaves[i+1:]
print('len(leaves) =', len(leaves))
print(leaves)
print('len(human_proteins) =', len(human_proteins))
x_coords = [float(leaves.index(x)) for x in range(len(human_proteins))]
y_coords = [float(0.0)] * len(human_proteins)

# however, it seems final branch coordinates are generated in icoord and dcoord
# the x_coords and y_coords are just used locally here
icoord = []
dcoord = []
for l in Z:
    x1 = x_coords[int(l[0])]
    x2 = x_coords[int(l[1])]
    y1 = y_coords[int(l[0])]
    y2 = y_coords[int(l[1])]
    icoord.append((x1, x1, x2, x2))
    dcoord.append((y1, l[2], l[2], y2))
    x_coords.append( (x1+x2)/2 )
    y_coords.append( l[2] )
print('len(icoord) =', len(icoord))
print('len(dcoord) =', len(dcoord))

# we scale the edges to make branches near the edge and the centre of the circle more visible
# here we use an arcsin-based function, but anything else could be used
bound = lambda x: 1.0 if x > 1.0 else (0.0 if x < 0.0 else x)
scale = lambda x: math.asin(bound(x)*2-1)/math.pi+0.5
#scale = lambda x: (np.sign(bound(x)*2-1)*math.pow(abs(bound(x)*2-1), 1.0/5.0)+1)*0.5
#scale = lambda x: (math.pow(bound(x)*2-1, 5.0)+1)*0.5
#scale0 = lambda x: (np.sign(x)*math.pow(abs(x), 1.0/3.0)+1)*0.5
#print scale0(-1.0), scale0(-0.5), scale0(0.0), scale0(0.5), scale0(1.0)
print(scale(0.0), scale(0.25), scale(0.5), scale(0.75), scale(1.0))

#scaled_dcoord = [[scale(x) for x in y] for y in dnd['dcoord']]
scaled_dcoord = [[scale(x) for x in y] for y in dcoord]

# SVG generation
svg_header = svg_header_base.format(1000, 1000)

# we build the colors from scratch
color_map = []
assigned_colors = dict()
cl_list = []
hp = []
for i in leaves:
    acc = human_proteins[i]
    hp.append(acc)
    # here we use the name cl_id because it doesn't index our branches any more!
    cl_id = slc_to_br[acc]
    #cl_name = br_to_slc_family.get(cl_id, str(cl_id))
    cl_list.append(cl_id)

# count consecutive family/cluster labels
cm = []
color_cycle = itertools.cycle(color_values[1:])
running_x = 0
discs = []
for cl_id, group in itertools.groupby(cl_list):
    count = sum(1 for x in group)
    cl_name = br_to_slc_family.get(cl_id, str(cl_id))
    if count == 1:
        color = color_values[0]
    elif cl_name in assigned_colors: color = assigned_colors[cl_name]
    else:
        color = next(color_cycle)
        if (len(cm) > 0) and (color == cm[-1]): color = next(color_cycle)
    if count > 1:
        # record here the data of the discs we will draw
        top = scaled_dcoord[cl_id-len(human_proteins)][1]
        discs.append([float(running_x)-0.5, -0.01, float(running_x+count-1)+0.5, top+0.01, color, cl_name])
    cm.extend([color]*count)
    running_x += count
    
human_protein_colors = dict(zip(hp, cm))

# we need this in original order!
color_map = [human_protein_colors[x] for x in human_proteins]

# now color branches
for l in Z:
    if (l[0] < 0) or (l[1] < 0): continue
    c1, c2 = map(lambda x: color_map[int(x)], l[:2])
    if c1 == c2:
        color_map.append(c1)
    else:
        color_map.append(color_values[0])

# circular dendrogram
cx, cy = 500.0, 500.0
radius = 400.0
outer_margin = 5.0

polar_icoord = np.array(icoord)/len(human_proteins) * 360.0
#icoord = ((np.array(icoord)-5.0)/Z[-1][-1]/10.0)*360.0
#icoord = ((np.array(dnd['icoord'])-5.0)/Z[-1][-1]/10.0)*360.0
#dcoord = 100.0-(np.array(dnd['dcoord'])*0.7*100.0)
# there is another scaling by 0.7 here to make a hole in the centre of the circle
polar_dcoord = 100.0-(np.array(scaled_dcoord)*0.7*100.0)

# transform
polar_dcoord = polar_dcoord*(radius-outer_margin)/100.0
# draw backwards to avoid overlap
svg_body = [branch(cx, cy, polar_icoord[x], polar_dcoord[x], color=color_map[len(human_proteins)+x]) 
            for x in range(len(human_proteins)-2, 0, -1)]
## the final line into the centre
#svg_body.append('<path d="M {start[0]:f} {start[1]:f} L {end[0]:f} {end[1]:f}" style="fill:none;stroke:#000000;"/>\n'.format(
#                start=p2c(cx, cy, dcoord[-1][1], (icoord[-1][1]+icoord[-1][2])/2.0),
#                end=(cx, cy))
#               )

# labels
label_margin = 5

svg_labels = [label(p2c(cx, cy, radius+label_margin,
                        i*360.0/len(human_proteins)),
                    i*360.0/len(human_proteins), 
                    human_proteins_names[x].replace("_Homo", ""), color=color_map[x], dy=2.0)
              for i, x in enumerate(leaves)]

# families (discs)

disc_label_margin = 40

svg_discs = []
svg_disc_labels = []
for disc in discs:
    phi = np.array([disc[0], disc[2]])
    r = np.array([disc[1], disc[3]])
    r = 100.0-(r*0.7*100.0)
    r = r*(radius-outer_margin)/100.0
    phi = phi/len(human_proteins) * 360.0
    svg_discs.append(sector(cx, cy, r[0], phi[0], r[1], phi[1], color=disc[4], opacity=0.2))
    if disc[5].startswith('SLC'):
        svg_disc_labels.append(label(p2c(cx, cy, r[0]+disc_label_margin, (phi[0]+phi[1])/2.0), (phi[0]+phi[1])/2.0,
                                     disc[5], color=disc[4], font_size='12pt', dy=7.0))

# write to file
if not read_only:
    with open('circular-dendrogram.svg', 'wt') as f:
        f.write(svg_header)
        f.write('<g id="families">')
        f.writelines(svg_discs)
        f.write('</g>')
        f.write('<g id="tree">')
        f.writelines(svg_body)
        f.write('</g>')
        f.write('<g id="slc_labels">\n')
        f.writelines(svg_labels)
        f.write('</g>\n')
        f.write('<g id="family_labels">\n')
        f.writelines(svg_disc_labels)
        f.write('</g>\n')
        f.write(svg_footer)


