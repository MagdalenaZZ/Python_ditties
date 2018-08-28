#!/usr/bin/env python3

import argparse
import os
import fastaq


parser = argparse.ArgumentParser(
    description = 'Makes plots of contig layout from nucmer files',
    usage = '%(prog)s [options] <ref name> <ref.fasta.fai> <nucmer info file> <read depth file> <outfile>')
parser.add_argument('refname')
parser.add_argument('ref_fai')
parser.add_argument('nucmer_fofn')
parser.add_argument('read_depth_file')
parser.add_argument('outfile')
options = parser.parse_args()




class Coords:
    def __init__(self, a, b):
        self.start = a
        self.end = b

    def __lt__(self, other):
        return self.start < other.start or (self.start == other.start and self.end > other.end)

    def __str__(self):
        return str(self.start) + ',' + str(self.end)



def nucmer_to_coords(refname, filename):
    f = fastaq.utils.open_file_read(filename)
    lines = f.readlines()
    fastaq.utils.close(f)
    i = 0
    while not lines[i].startswith('[S1]'):
        i += 1
    lines = [x.rstrip() for x in lines[i+1:]]
    coords = []

    for line in lines:
        a = line.split('\t')
        if a[11] == refname:
            s = min(int(a[0]), int(a[1]))
            e = max(int(a[0]), int(a[1]))
            coords.append(Coords(s, e))
    coords.sort()
    return coords


def coords_to_plot_positions(coords):
    pos = []

    for c in coords:
        print('loop start', c)
        for l in pos:
            print('   ', ' .. '.join([str(x) for x in l]))

        if len(pos) == 0:
            pos = [[c]]
            continue

        possible = [i for i in range(len(pos)) if c.start > pos[i][-1].end + 1]
        print('possible:', possible)
        if len(possible) == 0:
            print('   append', c)
            pos.append([c])
        else:
            first_i = possible[0]
            for  i in possible:
                if pos[i][-1].end < pos[first_i][-1].end:
                    first_i = i

            pos[first_i].append(c)

    return pos


def load_cov_plot(filename, refname):
    lines = fastaq.utils.syscall_get_stdout('tabix ' + filename + ' ' + refname)
    read_depth = []
    snps = []
    for line in lines:
        a = line.split()
        read_depth.append(int(a[2]))
        snps.append(int(a[3]))

    return read_depth, snps

def list_to_plot_coords(l, xmin, xmax, ymin, ymax, window=1):
    max_y = max(l)
    min_y = min(l)
    print(min_y, max_y)
    coords = []
    for i in range(0, len(l), window):
        y = sum(l[i:i+window]) / window
        y = 1 - (l[i] - min_y) / (max_y - min_y) # scaled between 0 and 1
        y = ymin + (y * (ymax - ymin)) # scaled to plot area
        x = i / len(l) # scaled between 0 and 1
        x = xmin + (x * (xmax - xmin))
        coords.append((x, y))

    return coords




# coords  = list of tuples [(x1, y1), (x2, y2) ...]
def svg_polygon(coords, fill_colour, border_colour, border_width = 1, opacity=-1):
    return_string = '<polygon points="' + ' '.join([str(x[0])+','+str(x[1]) for x in coords]) + '" ' \
            + 'fill="' + fill_colour + '" '

    if opacity != -1:
        return_string += 'fill-opacity="' + str(opacity) + '" '

    return_string += 'stroke="' + border_colour + '" ' \
            + 'stroke-width="' + str(border_width) + '" ' \
            + '/>'

    return return_string

def svg_polyline(coords, colour, thickness = 1):
    return '<polyline points="' + ' '.join([str(x[0])+','+str(x[1]) for x in coords]) + '" ' \
            + 'stroke="' + colour + '" ' \
            + 'stroke-width="' + str(thickness) + '" ' \
            + 'fill="none"' \
            + '/>'





coords = {}

f = fastaq.utils.open_file_read(options.nucmer_fofn)
for line in f:
    (assembler, filename) = line.rstrip().split()
    coords[assembler] = nucmer_to_coords(options.refname, filename)
fastaq.utils.close(f)


ref_lengths = {}
fastaq.tasks.lengths_from_fai(options.ref_fai, ref_lengths)
ref_length = ref_lengths[options.refname]
contig_height = 2
contig_space = 2
assembly_space = 4
contig_rows = sum([len(x) for x in coords.values()])
total_contig_height = contig_rows * (contig_height + contig_space) + len(coords) * assembly_space
plot_height = 50
plot_space = 10
height = total_contig_height + plot_height * 2 + plot_space * 2


f = fastaq.utils.open_file_write(options.outfile)
width = 500
print (r'''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC " -//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg width="''' + str(width) + '" height="' + str(height) + '">', file=f)

colours = {
    'PRICE': 'red',
    'Vicuna': 'blue',
    'IVA': 'green'
}

h = height
for assembler in ['PRICE', 'Vicuna', 'IVA']:
    positions = coords_to_plot_positions(coords[assembler])
    for l in positions:
        print(assembler, h)
        top = h - contig_height
        bottom = h
        for p in l:
            left = width * p.start / ref_length
            right = width * p.end / ref_length
            corners = [(left, bottom), (left, top), (right, top), (right, bottom)]
            print(svg_polygon(corners, colours[assembler], colours[assembler], border_width=0), file=f)
        h -= (contig_height + contig_space)
    h -= assembly_space


top = h - contig_height / 2
bottom = h
print(svg_polygon([(0, bottom), (0, top), (width, top), (width, bottom)], 'black', 'black'), file=f)
h -= plot_space


read_depth, snps = load_cov_plot(options.read_depth_file, options.refname)
read_depth_coords = list_to_plot_coords(read_depth, 0, width, h - plot_height, h)
print(svg_polyline(read_depth_coords, 'black', thickness=0.25), file=f)
print(svg_polyline([(0,h), (width, h)], 'gray', thickness=0.1), file=f)

h -= plot_space + plot_height
snps_coords = list_to_plot_coords(snps, 0, width, h - plot_height, h)
print(svg_polyline(snps_coords, 'black', thickness=0.25), file=f)
print(svg_polyline([(0,h), (width, h)], 'gray', thickness=0.1), file=f)


h -= plot_space + plot_height
snps_per_read = [snps[i] / read_depth[i] if read_depth[i] > 0 else 0 for i in range(len(snps))]
snps_per_read_coords =  list_to_plot_coords(snps_per_read, 0, width, h - plot_height, h)
print(svg_polyline(snps_per_read_coords, 'black', thickness=0.25), file=f)
print(svg_polyline([(0,h), (width, h)], 'gray', thickness=0.1), file=f)


print('</svg>', file=f)
fastaq.utils.close(f)




