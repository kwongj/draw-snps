#!/usr/bin/env python
# Script by Jason Kwong
# Show SNP positions along reference

# Use modern print function from python 3.x
from __future__ import print_function

# Import modules
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import svgwrite
import os
import sys
#import csv
#import subprocess
#from subprocess import Popen


# Usage
parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Show SNP positions along reference',
	usage='\n  %(prog)s --ref REF <alignment>')
parser.add_argument('aln', metavar='ALN', nargs=1, help='alignment of gene sequences in FASTA format (required)')
parser.add_argument('--ref', metavar='REF', nargs=1, required=True, help='reference sequence')
parser.add_argument('--out', metavar='FILE', default='snps.svg', help='output file')
parser.add_argument('--svgsize', metavar='WIDExHIGH', default='800x600', help='specify width and height of SVG in pixels (default="800x600")')
#parser.add_argument('--svgorder', metavar='FILE', help='specify file containing list of taxa (1 per line) in desired order')
parser.add_argument('--version', action='version', version=
	'%(prog)s v0.1\n'
	'Updated 11-Sep-2016 by Jason Kwong\n'
	'Dependencies: Python 2.x')
args = parser.parse_args()

ref = args.ref[0]


# Functions
# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);
	
# Check file exists
def check_file(f):
	if os.path.isfile(f) == False:
		err('ERROR: Cannot find "{}". Check file exists in the specified directory.'.format(f))

# SVG functions
def label(n,p):                 # Sequence ID labels
	dwg.add(dwg.text(n, insert=(width+15, (p*h)-(0.2*h)), fill=colour, style=font))

#def ticks(q):                   # Draw tick marks
#	count = interval
#	while count < seqlen:
#		dwg.add(dwg.line((count*s,((q)*h)), (count*s,((q)*h)-3), stroke=main_colour))
#		dwg.add(dwg.text((count/1000000.0), insert=((count*s)-6, ((q)*h)+13), fill=main_colour, style=font))
#		count = count + interval

#def box(p, h, width):			# Draw box
#	dwg.add(dwg.line((0,p*h), (width,p*h), stroke=colour))
#	dwg.add(dwg.line((width,p*h), (width,(p-0.8)*h), stroke=colour))
#	dwg.add(dwg.line((0,(p-0.8)*h), (width,(p-0.8)*h), stroke=colour))
#	dwg.add(dwg.line((0,(p-0.8)*h), (0,p*h), stroke=colour))
	
#def box(width, h):               # Box with 5px margin
#	dwg.add(dwg.line((0,0), ((width),0), stroke=colour))
#	dwg.add(dwg.line(((width),0), ((width),(numseqs*h)), stroke=colour))
#	dwg.add(dwg.line(((width),(numseqs*h)), (0,(numseqs*h)), stroke=colour))
#	dwg.add(dwg.line((0,(numseqs*h)), (0,0), stroke=colour))

def rect(x,p,w,bgcolour):
	dwg.add(dwg.rect(insert=(x, (p-0.8)*h), size=(w, h*0.8), fill=bgcolour))

for record in SeqIO.parse(ref, 'fasta'):
	refSEQ = record.seq
	refLEN = len(record.seq)

aln = open(args.aln[0], 'rU')
seqALN = SeqIO.parse(aln, 'fasta')
numseqs = len(list(seqALN))
#numseqs = len(list(SeqIO.parse(aln, 'fasta')))

# SVG
dwg = svgwrite.Drawing(args.out)

svgsize = args.svgsize.split('x',1)
width = int(svgsize[0])
height = int(svgsize[1])
s = width/float(int(refLEN))
h = height/float(numseqs)
fsize = int(h*0.9)
font = 'font-size:{}px; font-family:Arial'.format(fsize)
colour = 'black'
bgcolour = 'lightgrey'
blankcolour = 'white'
interval = 100		               # Tick interval in bp
p = 1

aln = open(args.aln[0], 'rU')
seqALN = SeqIO.parse(aln, 'fasta')
for record in seqALN:
	seqNAME = record.id
	print('Drawing {} ...'.format(seqNAME))
	rect(0,p,width,bgcolour)
	label(seqNAME, p)
	if len(record.seq) != refLEN:
		err('ERROR: all sequences in alignment must be the same length as the reference')
	for n in range(0,refLEN-1):
		if record.seq[n] == '-':
			dwg.add(dwg.line((n*s,p*h), (n*s,(p-0.8)*h), stroke=blankcolour))
		elif record.seq[n] != refSEQ[n]:
			dwg.add(dwg.line((n*s,p*h), (n*s,(p-0.8)*h), stroke=colour))
	p = p+1
#box(width, h)
dwg.save()
