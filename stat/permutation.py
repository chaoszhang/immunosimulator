import sys
import numpy as np
f = open(sys.argv[1], "r")
line = f.readline()
arr = line.split()
m = int(arr[0])
print m
line = f.readline()
arr = line.split()
n = int(arr[0])
print n
line = f.readline()
arr = line.split()
root = int(arr[0])
print root
names = []
ids = []
for i in xrange(m):
	line = f.readline()
	arr = line.split()
	names.append(arr[0])
	ids.append(int(arr[1]))
np.random.shuffle(ids)
for i in xrange(m):
	print names[i], ids[i]
for i in xrange(n):
	line = f.readline()
	arr = line.split()
	p = int(arr[0])
	c = int(arr[1])
	print p, c