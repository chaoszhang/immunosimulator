import sys
import numpy as np

names = []
present = []
with open(sys.argv[2], "r") as ins:
	isTag = True
	for line in ins:
		text = line.split()[0]
		if isTag:
			names.append(text[1:])
		isTag = not isTag

f = open(sys.argv[1], "r")
line = f.readline()
arr = line.split()
m = int(arr[0])
print len(names)
line = f.readline()
arr = line.split()
n = int(arr[0])
print n + (len(names) - m)
line = f.readline()
arr = line.split()
root = int(arr[0])
print root
ids = []
for i in xrange(m):
	line = f.readline()
	arr = line.split()
	present.append(arr[0])
	print arr[0], arr[1]
for e in [x for x in names if x not in present]:
	print e, hash(e)
for i in xrange(n):
	line = f.readline()
	arr = line.split()
	p = int(arr[0])
	c = int(arr[1])
	print p, c
for e in [x for x in names if x not in present]:
	print root, hash(e)