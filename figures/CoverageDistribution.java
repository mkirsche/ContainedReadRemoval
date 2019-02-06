/*
 * Takes in Minimap2 output from reads aligned to a reference and 
 * outputs a coverage distribution.  Each line is of the form:
 * <coverage> <number of reference bases>
 * This can then be passed to a python program for plotting.
 */
import java.util.*;
import java.io.*;
public class CoverageDistribution {
public static void main(String[] args) throws IOException
{
	String fn = args[0];
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	TreeMap<Long, Integer> delta = new TreeMap<Long, Integer>();
	TreeMap<String, Long> offsets = new TreeMap<String, Long>();
	long totOffset = 0;
	input = new Scanner(new FileInputStream(new File(fn)));
	while(input.hasNext())
	{
		PafEntry cur = new PafEntry(input.nextLine());
		if(!offsets.containsKey(cur.tname))
		{
			offsets.put(cur.tname, totOffset);
			totOffset += cur.tlength;
		}
		long curOffset = offsets.get(cur.tname);
		long a = cur.tstart + curOffset, b = cur.tend + curOffset;
		delta.put(a, delta.containsKey(a) ? (1 + delta.get(a)) : 1);
		delta.put(b+1, delta.containsKey(b+1) ? (-1 + delta.get(b+1)) : -1);
	}
	long last = 0;
	int coverage = 0;
	TreeMap<Integer, Long> covFreq = new TreeMap<Integer, Long>();
	for(long x : delta.keySet())
	{
		long numBP = x - last;
		covFreq.put(coverage, covFreq.containsKey(coverage) ? (numBP + covFreq.get(coverage)) : numBP);
		last = x;
		coverage += delta.get(x);
	}
	System.out.println("Coverage for " + fn);
	System.out.println("Coverage");
	System.out.println("Frequency (bp)");
	for(int x : covFreq.keySet())
	{
		System.out.println(x+" "+covFreq.get(x));
	}
	double mean = 0;
	long count = 0;
	for(int x : covFreq.keySet())
	{
		mean += 1.0 * x * covFreq.get(x);
		count  += covFreq.get(x);
	}
}
static class PafEntry
{
	String qname, tname;
	int qlength, qstart, qend, tlength, tstart, tend;
	char strand;
	int matches, blockLength, qual;
	PafEntry(String line)
	{
		String[] tokens = line.split("\t");
		qname = tokens[0];
		qlength = Integer.parseInt(tokens[1]);
		qstart = Integer.parseInt(tokens[2]);
		qend = Integer.parseInt(tokens[3]);
		strand = tokens[4].charAt(0);
		tname = tokens[5];
		tlength = Integer.parseInt(tokens[6]);
		tstart = Integer.parseInt(tokens[7]);
		tend = Integer.parseInt(tokens[8]);
		matches = Integer.parseInt(tokens[9]);
		blockLength = Integer.parseInt(tokens[10]);
		qual = Integer.parseInt(tokens[11]);
	}
}
}
