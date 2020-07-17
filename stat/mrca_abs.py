import sys
import numpy as np

def update(p, children, dist):
	depth = [[(p, 0)]]
	if p in children:
		depth += [update(c, children, dist) for c in children[p]]	
	for s1 in depth:
		for s2 in depth:
			if s1 != s2:
				for e1 in s1:
					for e2 in s2:
						dist[(e1[0], e2[0])] = (e1[1], e2[1])
	depth = [(e[0], e[1] + 1) for s in depth for e in s]
	return depth

f1 = open(sys.argv[1], "r")
line = f1.readline()
arr = line.split()
m1 = int(arr[0])
line = f1.readline()
arr = line.split()
n1 = int(arr[0])
line = f1.readline()
arr = line.split()
root1 = int(arr[0])
name1 = {}
for i in xrange(m1):
	line = f1.readline()
	arr = line.split()
	if arr[0] == "root":
		continue
	name1[arr[0]] = int(arr[1])
children1 = {}
for i in xrange(n1):
	line = f1.readline()
	arr = line.split()
	p = int(arr[0])
	c = int(arr[1])
	if p not in children1:
		children1[p] = [c]
	else:
		children1[p].append(c)
dist1 = {}
update(root1, children1, dist1)

f2 = open(sys.argv[2], "r")
line = f2.readline()
arr = line.split()
m2 = int(arr[0])
line = f2.readline()
arr = line.split()
n2 = int(arr[0])
line = f2.readline()
arr = line.split()
root2 = int(arr[0])
name2 = {}
for i in xrange(m2):
	line = f2.readline()
	arr = line.split()
	if arr[0] == "root":
		continue
	name2[arr[0]] = int(arr[1])
children2 = {}
for i in xrange(n2):
	line = f2.readline()
	arr = line.split()
	p = int(arr[0])
	c = int(arr[1])
	if p not in children2:
		children2[p] = [c]
	else:
		children2[p].append(c)
dist2 = {}
update(root2, children2, dist2)

total = 0
for a in name1:
	for b in name1:
		if a < b:
			d1 = dist1[(name1[a], name1[b])]
			d2 = dist2[(name2[a], name2[b])]
			total += abs(d1[0] - d2[0]) + abs(d1[1] - d2[1])

print total