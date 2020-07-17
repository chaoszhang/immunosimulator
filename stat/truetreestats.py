import sys
import numpy as np

children = {}

name2seq = {}
names = []
relation = []
used = set()
remap = {}

def dist(s1, s2):
	return len([i for i in range(len(s1)) if s1[i] != s2[i]])
	
def f(n, depth):
	if n in children:
		sd = (depth if name2seq[n] in used else 0)
		md = 0
		sbl = 0
		nlb = 0
		for c in children[n]:
			d = dist(name2seq[n], name2seq[c])
			if d != 1:
				nlb += 1
			rt = f(c, depth + d)
			sd += rt[0]
			md = max(md, rt[1])
			sbl += d + rt[2]
			nlb += rt[3]
		return (sd, md, sbl, nlb)
	else:
		return (depth, depth, 0, 0)
		
def fphy(n, depth):
	if n in children:
		sd = (depth if name2seq[n] in used else 0)
		md = 0
		sbl = 0
		nlb = 0
		for c in children[n]:
			d = dist(name2seq[n], name2seq[c])
			if d != 1:
				nlb += 1
			d = -np.log(1 - 4 / 3 * d / len(name2seq[n]))
			rt = fphy(c, depth + d)
			sd += rt[0]
			md = max(md, rt[1])
			sbl += d + rt[2]
			nlb += rt[3]
		return (sd, md, sbl, nlb)
	else:
		return (depth, depth, 0, 0)
	
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

root = relation[0][0]
for e in relation:
	if e[0] in children:
		children[e[0]].append(e[1])
	else:
		children[e[0]] = [e[1]]

sumDepth, maxDepth, sumBranchLength, nLongBranch = f(root, 0)
print(sumDepth, maxDepth, sumBranchLength, nLongBranch)
sumDepth, maxDepth, sumBranchLength, nLongBranch = fphy(root, 0)
print(sumDepth, maxDepth, sumBranchLength, nLongBranch)