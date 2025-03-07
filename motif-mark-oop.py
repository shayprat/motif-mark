#!/usr/bin/env python

import argparse
import cairo
import re


def get_args():
    parser = argparse.ArgumentParser(description="A script to identiyfy and visualize motifs from a fasta file")
    parser.add_argument("-f", help="absolute filepath for input fasta file", required= True)
    parser.add_argument("-m", help="absolute filepath input motifs_file", required = True)
    
    return parser.parse_args()

args = get_args()
fasta_file = args.f 
motif_file = args.m


#list of colors using RGB codes
color_list = [
    (67/255, 114/255, 196/255), #blue
    (255/255, 192/255, 3/255), #yellow
    (237/255, 124/255, 49/255), #orange
    (220/255, 30/255, 200/255) #magenta
]

#pycairo image scaling
FEATURE_HEIGHT = 15
ALIGN = 40
IMAGE_HEIGHT = 70

#Classes

class Motif:
    def __init__(self, start, stop, gene_count, color):
       
        self.start = start
        self.stop = stop
        self.color = color 
        self.gene_count = gene_count

    def draw(self, ctx):
        ctx.set_source_rgb(self.color[0], self.color[1], self.color[2]) 
        y = (IMAGE_HEIGHT * self.gene_count) + ALIGN/2
        ctx.rectangle(self.start + ALIGN, y + ALIGN*1.3, self.stop - self.start, FEATURE_HEIGHT)        #(x0,y0,x1,y1) (moves left to right, moves up and down, width of box, height of box)
        ctx.fill()


class Gene:
    def __init__(self, start, stop, gene_name, gene_count):
    
        self.start = start
        self.stop = stop
        self.gene_name = gene_name
        self.gene_count = gene_count

    def draw(self, ctx):
        ctx.set_source_rgb(0, 0, 0) #black line
        ctx.set_line_width(1)
        y_axis = (IMAGE_HEIGHT * self.gene_count) + ALIGN*2
        ctx.move_to(self.start + ALIGN, y_axis)        
        ctx.line_to(self.stop + ALIGN, y_axis)      
        ctx.stroke()

        ctx.select_font_face("Arial", cairo.FONT_WEIGHT_BOLD)
        ctx.set_font_size(12)
        ctx.move_to(IMAGE_HEIGHT/5, y_axis - ALIGN/2)
        ctx.show_text(self.gene_name)



class Exon:

    EXON_HEIGHT = 10
    
    def __init__(self, start, stop, gene_count):
        
        self.start = start
        self.stop = stop
        self.gene_count = gene_count

    def draw(self, ctx):
        ctx.set_source_rgb(0, 0, 0) #black rectangle
        y_axis = (IMAGE_HEIGHT * self.gene_count) + ALIGN/2
        ctx.rectangle(self.start + ALIGN, y_axis + ALIGN*1.3, self.stop - self.start, FEATURE_HEIGHT) 
        ctx.fill()

#find the length of the exon, this will be needed to draw rectangles!
def find_exon(sequence: str):
    for i in range(len(sequence)):
        if sequence[i].isupper():
            start = int(i)
            break
    
    for i in range(start, len(sequence)): 
        if sequence[i].islower():
            stop = int(i) 
            break

    return (start, stop)

#making a dictionary linking the degenerate/ambiguous base with it's potential DNA base
ambiguous_nt_dict = {
    "W":"[A|T]",
    "S":"[C|G]",
    "M":"[A|C]",
    "K":"[G|T]",
    "R":"[A|G]",
    "Y":"[C|T]",
    "B":"[C|G|T]",
    "D":"[A|G|T]",
    "H":"[A|C|T]",
    "V":"[A|C|G]",
    "N":"[A|C|G|T]",
}

def non_ambig_motifs(motif: str):
    '''This function generates '''
    motif = motif.upper()
    for degen_base, corresponding_base in ambiguous_nt_dict.items(): 
        motif = motif.replace(degen_base, corresponding_base) 

    return motif

motif_dict = {}
with open(motif_file, 'r') as motif_txt:
    for line in motif_txt:
        line = line.strip()
        if line == "":
            break
        motif_dict[line] = non_ambig_motifs(line)
#print(motif_dict)


fasta_dict = {} #{value = header : key = sequence}
with open(fasta_file, "r") as input_fasta:
    while True:
        line = input_fasta.readline().strip()
        if line == "":
            break
        if line.startswith(">"):
            header = line
            fasta_dict[header] = ""
        else: 
            fasta_dict[header] += line
#print(fasta_dict)

#width of image surface
width = 800 + ALIGN

#height of image surface 
height = (IMAGE_HEIGHT * len(fasta_dict)) +100

#create image surface 
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
ctx = cairo.Context(surface)
ctx.set_source_rgb(1, 1, 1) 
ctx.paint()


for i, header in enumerate(fasta_dict): 
    gene = Gene(0, len(fasta_dict[header]), header, i)
    gene.draw(ctx)

    #find start and end of exon 
    (start, stop) = find_exon(fasta_dict[header])

    exon = Exon(start, stop, i)
    exon.draw(ctx)

    for j, motif in enumerate(motif_dict):
        matches = re.finditer(motif_dict[motif], fasta_dict[header].upper())
        for match in matches:
            start = match.start()
            stop = match.start() + len(match.group())
            motif = Motif(start, stop, i, color_list[j])
            motif.draw(ctx)


#context.set_source_rgba(0,0,0) #black
#context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
                  
#legend
x_start = width - 825
ctx.set_font_size(10)
ctx.move_to(x_start, 40)

#motif legend
for i, motif in enumerate(motif_dict):
    #set font & color key
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL)
    ctx.set_source_rgb(color_list[i][0], color_list[i][1], color_list[i][2])
    y_start = 320
    ctx.rectangle(x_start, y_start +25, 9, 15)        
    ctx.fill()
    
    #add motif seq 
    ctx.set_font_size(10)
    ctx.move_to(x_start + 15, y_start + 35)
    ctx.show_text(motif)

    x_start += 60

# Exon legend
ctx.set_font_size(12)
ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL)
ctx.set_source_rgb(0, 0, 0)
ctx.rectangle(x_start + 35, y_start +25, 30, 15)
ctx.fill()
ctx.move_to(x_start + 70  , y_start + 36)
ctx.show_text("Exon")

#Gene legend
ctx.set_font_size(12)
ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL)
ctx.set_source_rgb(0, 0, 0)
ctx.move_to(x_start + 105, y_start + 33)
ctx.line_to(x_start + 125, y_start +33)
ctx.stroke()
ctx.move_to(x_start + 128 , y_start + 36)
ctx.show_text("Gene")


prefix = fasta_file.split(".")
surface.write_to_png(f"{prefix[0]}.png")