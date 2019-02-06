import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fn = sys.argv[1]
ofn = sys.argv[2]

xs = []
ys = []
lines = []
with open(fn) as f:
    lines = f.readlines()
    for line in lines[3:]:
        a, b = line.split()
        xs.append(int(a))
        ys.append(int(b))

maxx = xs[len(xs)-1]
plt.hist(xs, bins = [i for i in range(0, 100)], weights=ys)
plt.title(lines[0])
plt.xlabel(lines[1])
plt.ylabel(lines[2])
plt.savefig(ofn)
