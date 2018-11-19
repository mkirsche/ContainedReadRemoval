import matplotlib.pyplot as plt
import sys
import numpy as np
import os.path

if len(sys.argv) < 2:
    print('Usage: python ' + sys.argv[0] + ' [reads filename] [non-contained reads filename]')
    sys.exit(0)
    
fn1 = sys.argv[1]
fn2 = sys.argv[2]

with open(fn1) as f1, open(fn2) as f2:
    lines1 = f1.readlines()
    lines2 = f2.readlines()
    set1 = []
    set2 = []
    for i in range(1, len(lines1), 2):
        set1.append(len(lines1[i]))
    for i in range(1, len(lines2), 2):
        set2.append(len(lines2[i]))
    minLength = min(min(set1), min(set2))
    maxLength = max(max(set1), max(set2))
    bins = np.arange(8000, 20000, 1000)
    print bins
    plt.figure()
    plt.hist(set1, bins, stacked=True, label = 'Full data')
    plt.hist(set2, bins, stacked=True, label = 'Subset')
    plt.legend(loc='upper right')
    plt.show()
