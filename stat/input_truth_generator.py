import sys
import numpy as np
from sets import Set
from heapq import *

name2seq = {}
names = []
relation = []
used = Set()
remap = {}

with open(sys.argv[1], "r") as ins:
	for line in ins:
		text = line.split()
		relation.append((text[0], text[2]))
		seq = text[1].split("_")[-1]
		name2seq[text[2]] = seq
		if text[1][0] == "p" and seq not in used:
			names.append((text[2], seq))
			used.add(seq)

relation = relation[::-1]
name2seq[relation[0][0]] = sys.argv[2]
for e in relation:
	if name2seq[e[0]] == name2seq[e[1]]:
		remap[e[0]] = e[1]
relation = [e for e in relation if name2seq[e[0]] != name2seq[e[1]]]
relation = [(remap[e[0]], e[1]) if e[0] in remap else e for e in relation]
relation = [(e[0], remap[e[1]]) if e[1] in remap else e for e in relation]

print ">root"
print sys.argv[2]

for e in names:
	print ">" + e[0]
	print e[1]

print >>sys.stderr, len(names)
print >>sys.stderr, len(relation)
print >>sys.stderr, hash(relation[0][0])
for e in names:
	print >>sys.stderr, e[0], hash(e[0])
for e in relation:
	print >>sys.stderr, hash(e[0]), hash(e[1])