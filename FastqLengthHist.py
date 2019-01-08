# Generates a histogram of read lengths from a file in FASTQ format
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print 'Usage: python FastqLengthHist.py <readsfile> <outputfile>'
    sys.exit()
    
fn = sys.argv[1]
ofn = sys.argv[2]

lens = []

with open(fn) as f:
    idx = 0
    for line in f:
        idx += 1
        if idx%4 == 2:
            lens.append(len(line))
            
plt.hist(lens)
plt.savefig(ofn)
            
