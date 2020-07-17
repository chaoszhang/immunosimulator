import sys
import numpy as np
from sets import Set

def update(p, children, cluster, rname):
	if p in rname:
		h = hash(rname[p])
	else:
		h = 0
	if p in children:
		for c in children[p]:
			h += update(c, children, cluster, rname)	
	cluster.add(h)
	return h

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
rname1 = {}
for i in xrange(m1):
	line = f1.readline()
	arr = line.split()
	if arr[0] == "root":
		continue
	name1[arr[0]] = int(arr[1])
	rname1[int(arr[1])] = arr[0]
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
cluster1 = Set()
update(root1, children1, cluster1, rname1)

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
rname2 = {}
for i in xrange(m2):
	line = f2.readline()
	arr = line.split()
	if arr[0] == "root":
		continue
	name2[arr[0]] = int(arr[1])
	rname2[int(arr[1])] = arr[0]
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
cluster2 = Set()
update(root2, children2, cluster2, rname2)

both = len([x for x in cluster1 if x in cluster2])

print both, len(cluster1) - both, len(cluster2) - both