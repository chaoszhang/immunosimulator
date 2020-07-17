import sys
import numpy as np
from sets import Set

rotate_table = {}
format_table = {}
formated_set = Set()
def is_proper(t):
	for i in range(10):
		if (t >> (3 * i)) & 7 not in [0, 1, 2, 4]:
			return False
	return True
	
def rotate_subtree(s):
	n1 = s & 7
	n2 = (s >> 3) & 7
	n3 = (s >> 6) & 7
	return n1 | (max(n2, n3) << 3) | (min(n2, n3) << 6)

def rotate(t):
	if t in rotate_table:
		return rotate_table[t]
	r = t & 7
	s1 = rotate_subtree((t >> 3) & 511)
	s2 = rotate_subtree((t >> 12) & 511)
	s3 = rotate_subtree((t >> 21) & 511)
	
	ns1 = max(s1, s2, s3)
	ns3 = min(s1, s2, s3)
	ns2 = s1 + s2 + s3 - ns1 - ns3
	
	rotate_table[t] = r | (ns1 << 3) | (ns2 << 12) | (ns3 << 21)
	return rotate_table[t]
	
def format_subtree(s):
	s = rotate_subtree(s)
	n1 = s & 7
	n2 = (s >> 3) & 7
	n3 = (s >> 6) & 7
	
	if n1 == 0 and n3 == 0:
		return n2
	
	return s

def format(t):
	if t in format_table:
		return format_table[t]
	r = t & 7
	s1 = format_subtree((t >> 3) & 511)
	s2 = format_subtree((t >> 12) & 511)
	s3 = format_subtree((t >> 21) & 511)
	
	ns1 = max(s1, s2, s3)
	ns3 = min(s1, s2, s3)
	ns2 = s1 + s2 + s3 - ns1 - ns3
	
	if r == 0 and ns2 == 0:
		format_table[t] = (ns1 & 7) | (((ns1 >> 3) & 7) << 3) | (((ns1 >> 6) & 7) << 12)
	else:
		format_table[t] = r | (ns1 << 3) | (ns2 << 12) | (ns3 << 21)
	formated_set.add(format_table[t])
	return format_table[t]

dist_list = {}
def	adjacency(t, df, b):
	if b == 0:
		if df == 0:
			return
		t0 = format(t)
		t1 = format(t ^ df)
		dist_list[(t0, t1)] = 1
		dist_list[(t1, t0)] = 1
		dist_list[(t0, t0)] = 0
		dist_list[(t1, t1)] = 0
	else:
		for p in range(10):
			s = b << (3 * p)
			adjacency(t | s, df, b >> 1)
			if df == 0 and p in [1, 4, 7]:
				adjacency(t | s, b | s, b >> 1)
				adjacency(t | s, s * 9, b >> 1)
				adjacency(t | s, s * 65, b >> 1)

def update(p, children, mrcad, depth):
	groups = [[p]]
	mrcad[(p, p)] = depth
	if p in children:
		groups += [update(c, children, mrcad, depth + 1) for c in children[p]]	
	for s1 in groups:
		for s2 in groups:
			if s1 != s2:
				for e1 in s1:
					for e2 in s2:
						mrcad[(e1, e2)] = depth
	return [e for s in groups for e in s]

def depth2triplet(ab, bc, ac, aa, bb, cc):
	v = 0
	if ab == bc and bc == ac:
		v += (1 if aa == ab else (1 << 3))
		v += (2 if bb == ab else (2 << 12))
		v += (4 if cc == ab else (4 << 21))
	if ab > bc:
		v += ((1 << 3) if aa == ab else (1 << 6))
		v += ((2 << 3) if bb == ab else (2 << 9))
		v += (4 if cc == bc else (4 << 12))
	if bc > ab:
		v += (1 if aa == ab else (1 << 12))
		v += ((2 << 3) if bb == bc else (2 << 6))
		v += ((4 << 3) if cc == bc else (4 << 9))
	if ac > ab:
		v += ((1 << 3) if aa == ac else (1 << 6))
		v += (2 if bb == ab else (2 << 12))
		v += ((4 << 3) if cc == ac else (4 << 9))
	return format(v)
	
adjacency(0, 0, 4)
for k in formated_set:
	for i in formated_set:
		if (i, k) in dist_list:
			for j in formated_set:
				if (k, j) in dist_list:
					if (i, j) not in dist_list or dist_list[(i, j)] > dist_list[(i, k)] + dist_list[(k, j)]:
						dist_list[(i, j)] = dist_list[(i, k)] + dist_list[(k, j)]
#print len(dist_list)
#print len([e[0] for e in dist_list if is_proper(e[0]) and is_proper(e[1])])
#print Set([e[0] for e in dist_list if is_proper(e[0]) and is_proper(e[1])])
#print dist_list[(273, 84)]
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
mrcad1 = {}
update(root1, children1, mrcad1, 0)

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
mrcad2 = {}
update(root2, children2, mrcad2, 0)

names = [s for s in name1]
n = len(names)
total = 0
zerocnt = 0
for i in range(n):
	for j in range(i):
		for k in range(j):
			a1 = name1[names[i]]
			b1 = name1[names[j]]
			c1 = name1[names[k]]
			a2 = name2[names[i]]
			b2 = name2[names[j]]
			c2 = name2[names[k]]
			t = dist_list[(depth2triplet(mrcad1[(a1, b1)], mrcad1[(b1, c1)], mrcad1[(a1, c1)], mrcad1[(a1, a1)], mrcad1[(b1, b1)], mrcad1[(c1, c1)]),
				depth2triplet(mrcad2[(a2, b2)], mrcad2[(b2, c2)], mrcad2[(a2, c2)], mrcad2[(a2, a2)], mrcad2[(b2, b2)], mrcad2[(c2, c2)]))]
			total += t
			if t == 0:
				zerocnt += 1
print total, zerocnt