import sys
import numpy as np

children = {}
name = {}
hasname = {}

def fnNodes(p):
	if p in children:
		return 1 + sum([fnNodes(c) for c in children[p]])	
	return 1

def fnLeaves(p):
	if p in children:
		return sum([fnLeaves(c) for c in children[p]])	
	return 1

def fnInternalLabels():
	return len([n for n in name if name[n] in children])

def fnCherry(p):
	if p in children:
		n = len([c for c in children[p] if c not in children])
		return (n * (n - 1) / (len(children[p]) - 1) if len(children[p]) > 1 else n) / 2 + sum([fnCherry(c) for c in children[p]])	
	return 0

def fmaxHeight(p):
	if p in children:
		return 1 + max([fmaxHeight(c) for c in children[p]])	
	return 0
	
def fsumHeight(p):
	if p in children:
		lst = [fsumHeight(c) for c in children[p]]
		s = sum([e[0] + e[1] for e in lst])
		cnt = sum([e[1] for e in lst])
		if p in hasname:
			cnt += 1
		return (s, cnt)
	return (0, 1)
	
f = open(sys.argv[1], "r")
line = f.readline()
arr = line.split()
m = int(arr[0])
line = f.readline()
arr = line.split()
n = int(arr[0])
line = f.readline()
arr = line.split()
root = int(arr[0])
for i in range(m):
	line = f.readline()
	arr = line.split()
	if arr[0] == "root":
		continue
	name[arr[0]] = int(arr[1])
	hasname[int(arr[1])] = arr[0]
for i in range(n):
	line = f.readline()
	arr = line.split()
	p = int(arr[0])
	c = int(arr[1])
	if p not in children:
		children[p] = [c]
	else:
		children[p].append(c)

print(fnNodes(root), fnLeaves(root), fnInternalLabels(), fnCherry(root), fmaxHeight(root), fsumHeight(root)[0])
