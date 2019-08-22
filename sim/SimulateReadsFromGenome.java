/*
 * Simulates reads from a genome and computes containment scores for each of them.
 * 
 * For contained reads, the score is the maximum harmonic mean of extra bases on each end of any containing read
 * For uncontained reads, it's the number of bases in this read which are not contained in the read r 
 *   such that r contains the highest proportion of this read.  It is also multiplied by -1.
 */

import java.util.*;
import java.io.*;
import org.apache.commons.math3.distribution.NormalDistribution;

public class SimulateReadsFromGenome {
	static long maxLength = 100000000;
	static int coverage = 30;
	static double mean = 15000;
	static double stdev = 1500;
	static double errorRate = 0.12;
	static NormalDistribution nd;
	static String genomeFn, readOfn, scoreOfn;
	static String sampleFn;
	static Random r;
	static String outRefFn;
	static  Distribution dist;
	static String distType;
	static String intervalOfn;
public static void main(String[] args)  throws IOException
{
	genomeFn = "/home/mkirsche/references/genome.fa";
	readOfn = "simulatedreads.fa";
	scoreOfn = "simulatedscores.txt";
	outRefFn = "simulatedgenome.txt";
	intervalOfn = "simulatedintervals.txt";
	sampleFn = "";
	
	int argParseErrorKey = parseArgs(args);
	if(argParseErrorKey != 0)
	{
		return;
	}
	
	if(sampleFn.length() > 0)
	{
		dist = getDistribution(sampleFn);
	}
	
	Scanner input = new Scanner(new FileInputStream(new File(genomeFn)));
	PrintWriter readOut = new PrintWriter(new File(readOfn));
	PrintWriter scoreOut = new PrintWriter(new File(scoreOfn));
	PrintWriter genomeOut = new PrintWriter(new File(outRefFn));
	PrintWriter intervalOut = new PrintWriter(new File(intervalOfn));
	
	nd = new NormalDistribution(mean, stdev);
	r = new Random(333);
	
	ArrayList<Pair> coords = simulatePositions();
	System.err.println("Simulated " + coords.size() + " reads");
	
	System.err.println("Printing reads");
	printReads(input, coords, readOut, genomeOut);
	
	System.err.println("Printing intervals");
	for(int i = 0; i<coords.size(); i++)
	{
		intervalOut.println(getName(coords.get(i).i)+" "+coords.get(i).a+" "+coords.get(i).b);
	}
	
	System.err.println("Computing scores");
	printContainmentScores(coords, scoreOut);
	
	genomeOut.close();
	scoreOut.close();
	readOut.close();
	intervalOut.close();
	
	PrintWriter sortedOut = new PrintWriter(new File("sorted.txt"));
	Collections.sort(coords, new Comparator<Pair>() {

		@Override
		public int compare(Pair a, Pair b) {
			long alength = a.b - a.a;
			long blength = b.b - b.a;
			return Long.compare(blength, alength);
		}
	});
	
	int test = 100000;
	boolean[] done = new boolean[coords.size()];
	int[] counts = new int[Math.min(coords.size(), test)];
	long[] cumulength = new long[counts.length];
	int[] tots = new int[counts.length];
	int tot = 0;
	IntervalUnion iu = new IntervalUnion();
	for(int i = 0; i<test && i < coords.size(); i++)
	{
		int count = 0;
		for(int j = i+1; j<coords.size(); j++)
		{
			if(coords.get(j).a >= coords.get(i).a && coords.get(j).b <= coords.get(i).b)
			{
				count++;
				if(!done[j])
				{
					tot++;
				}
				done[j] = true;
			}
		}
		iu.add(coords.get(i).a, coords.get(i).b);
		counts[i] = count;
		tots[i] = tot;
		cumulength[i] = iu.length;
	}
	for(int i = 0; i<coords.size(); i++) iu.add(coords.get(i).a, coords.get(i).b);
	System.out.println(iu.length);
	for(int i = 0; i<test && i < coords.size(); i++)
	{
		sortedOut.println(getName(coords.get(i).i)+" "+coords.get(i).a+" "+coords.get(i).b+" "+counts[i]+" "+tots[i]+" "+cumulength[i]);
	}
	sortedOut.close();
}
static class IntervalUnion
{
	TreeSet<Interval> set;
	long length;
	
	IntervalUnion()
	{
		set = new TreeSet<>();
		length = 0;
	}
	
	class Interval implements Comparable<Interval>
	{
		long a, b;
		Interval(long aa, long bb)
		{
			a = aa; b = bb;
		}
		@Override
		public int compareTo(Interval o) {
			// TODO Auto-generated method stub
			if(a != o.a) return (int) (a - o.a);
			return (int) (b - o.b);
		}
	}
	
	void add(long a, long b)
	{
		long removedLength = 0;
		Interval toAdd = new Interval(a, b);
		Interval floor = set.floor(toAdd);
		while(floor != null)
		{
			if(floor.b >= toAdd.a)
			{
				toAdd.a = Math.min(toAdd.a, floor.a);
				toAdd.b = Math.max(toAdd.b, floor.b);
				removedLength += floor.b - floor.a + 1;
				set.remove(floor);
				floor = set.floor(toAdd);
			}
			else
			{
				break;
			}
		}
		Interval ceiling = set.ceiling(toAdd);
		while(ceiling != null)
		{
			if(ceiling.a <= toAdd.b)
			{
				toAdd.b = Math.max(toAdd.b, ceiling.b);
				toAdd.a = Math.min(toAdd.a, ceiling.a);
				removedLength += ceiling.b - ceiling.a + 1;
				set.remove(ceiling);
				ceiling = set.ceiling(toAdd);
			}
			else
			{
				break;
			}
		}
		length += toAdd.b - toAdd.a + 1 - removedLength;
		set.add(toAdd);
	}
}
static String meanKey = "mean";
static String stdevKey = "stdev";
static String genomeKey = "genomefile";
static String readKey = "readfile";
static String scoreKey = "scorefile";
static String intervalKey = "intervalfile";
static String coverageKey = "coverage";
static String maxLengthKey = "maxlength";
static String errorKey = "error";
static String sampleKey = "sample";
static String distKey = "dist";
static int parseArgs(String[] args)
{
	if(args.length == 0)
	{
		return 0;
	}
	String firstArg = args[0].replaceAll("-", "");
	if(firstArg.equals("help") || firstArg.equals("h"))
	{
		help();
		return 1;
	}
	for(String arg : args)
	{
		arg = arg.replace("-", "");
		if(arg.indexOf('=') == -1)
		{
			continue;
		}
		
		String key = arg.substring(0, arg.indexOf('=')).toLowerCase();
		String value = arg.substring(1 + arg.indexOf('='));
		
		if(key.equals(genomeKey))
		{
			genomeFn = value;
		}
		else if(key.equals(readKey))
		{
			readOfn = value;
		}
		else if(key.equals(scoreKey))
		{
			scoreOfn = value;
		}
		else if(key.equals(intervalKey))
		{
			intervalOfn = value;
		}
		else if(key.equals(sampleKey))
		{
			sampleFn = value;
		}
		else if(key.equals(meanKey))
		{
			mean = Double.parseDouble(value);
		}
		else if(key.equals(stdevKey))
		{
			stdev = Double.parseDouble(value);
		}
		else if(key.equals(maxLengthKey))
		{
			maxLength = Long.parseLong(value);
		}
		else if(key.equals(coverageKey))
		{
			coverage = Integer.parseInt(value);
		}
		else if(key.equals(errorKey))
		{
			errorRate = Double.parseDouble(value);
			if(errorRate > 1)
			{
				errorRate /= 100;
			}
		}
		else if(key.equals("distKey"))
		{
			distType = value;
		}
	}
	
	return 0;
}
static void help()
{
	System.out.println("Usage: java SimulateReadsFromGenome [optional parameters]");
	System.out.println("\nParameters:");
	System.out.println("  " + genomeKey + ": fasta file with genome to simulate from");
	System.out.println("  " + readKey + ": file name to output reads to");
	System.out.println("  " + scoreKey + ": file name to output containment scores to");
	System.out.println("  " + intervalKey + ": file name to output read intervals to to");
	System.out.println("  " + meanKey + ": average read length");
	System.out.println("  " + stdevKey + ": standard deviation of read length");
	System.out.println("  " + coverageKey + ": coverage to simulate");
	System.out.println("  " + maxLengthKey + ": maximum prefix of genome to simulate from");
	System.out.println("  " + errorRate + ": error rate of reads relative to the reference");
	System.out.println("  " + sampleKey + ": file to draw sample reads from to use a real distribution");
	System.out.println("  " + distKey + ": the distribution to draw reads from if no sample specified (normal)");
}
static ArrayList<Pair> simulatePositions()
{
	ArrayList<Pair> res = new ArrayList<>();
	long totalLength = 0;
	long neededLength = 1L * maxLength * coverage;
	int readCount = 0;
	while(totalLength < neededLength)
	{
		Pair cur = sampleRead(readCount);
		readCount++;
		totalLength += cur.b - cur.a;
		res.add(cur);
	}
	return res;
}
static int sampleLength()
{
	if(dist != null)
	{
		return dist.sample();
	}
	return (int)(nd.sample() + 0.5);
}
static Pair sampleRead(int i)
{
	int length = sampleLength();
	long start = r.nextLong() % (maxLength - length);
	if(start < 0) start += maxLength - length;
	long end = start + length;
	return new Pair(start, end, i);
}
static void printContainmentScores(ArrayList<Pair> coords, PrintWriter out)
{
	Collections.sort(coords);
	int n = coords.size();
	
	int radius = n;
	
	for(int i = 0; i<n; i++)
	{
		int maxContainmentScore = -1;
		int minOverhang = (int)(coords.get(i).b - coords.get(i).a);
		int bestContaining = i;
		for(int j = Math.max(0, i - radius); j <= i + radius && j < n; j++)
		{
			if(j == i)
			{
				continue;
			}
			if(coords.get(j).a <= coords.get(i).a && coords.get(j).b >= coords.get(i).b)
			{
				// Contained so update score
				int curScore = score(coords.get(i).a, coords.get(i).b, coords.get(j).a, coords.get(j).b);
				if(curScore > maxContainmentScore)
				{
					bestContaining = j;
				}
				maxContainmentScore = Math.max(maxContainmentScore, curScore);
			}
			else
			{
				// Not contained so update max overhang
				int curScore = (int)(Math.max(0, coords.get(j).a - coords.get(i).a) + Math.max(0, coords.get(i).b - coords.get(j).b));
				if(maxContainmentScore < 0 && curScore < minOverhang)
				{
					bestContaining = j;
				}
				minOverhang = Math.min(minOverhang, curScore);
			}
		}
		out.print(getName(coords.get(i).i) + " ");
		int score = maxContainmentScore >= 0 ? maxContainmentScore : -minOverhang;
		coords.get(i).score = score;
		out.print(score+" ");
		out.print(coords.get(i).a+" ");
		out.print(coords.get(i).b+" ");
		out.println(getName(coords.get(bestContaining).i));
	}
}
static int score(long innera, long innerb, long outera, long outerb)
{
	long extraa = innera - outera;
	long extrab = outerb - innerb;
	return (int)(.5 + harmonicMean(extraa, extrab));
}
static double harmonicMean(double x, double y)
{
	return 2.0 * x * y / (x+y);
}
static void printReads(Scanner input, ArrayList<Pair> coords, PrintWriter out, PrintWriter genomeOut)
{
	int bufferLength = 500000;
	Collections.sort(coords, new Comparator<Pair>() {
		public int compare(Pair a, Pair b) {
			
			if(a.b != b.b) return Long.compare(a.b, b.b);
			return Long.compare(a.a, b.a);
		}
	});
	
	int characterIndex = 0;
	int coordsIndex = 0;
	
	ArrayDeque<Character> buffer = new ArrayDeque<>();
	genomeOut.println(">genome");
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.startsWith(">"))
		{
			continue;
		}
		
		for(int i = 0; i < line.length(); i++)
		{
			char c = line.charAt(i);
			if(!isBasePair(c)) continue;
			
			buffer.addLast(c);
			
			genomeOut.print(c);
			
			if(buffer.size() > bufferLength)
			{
				buffer.pollFirst();
			}
			
			characterIndex++;
			
			while(coordsIndex < coords.size() && characterIndex == coords.get(coordsIndex).b)
			{
				String simulatedRead = getSuffix(buffer, (int)(coords.get(coordsIndex).b - coords.get(coordsIndex).a));
				out.println(">" + getName(coords.get(coordsIndex).i));
				out.println(simulatedRead);
				coordsIndex++;
			}
			
			if(characterIndex == maxLength)
			{
				break;
			}
		}
		
		if(characterIndex == maxLength)
		{
			break;
		}
	}
	genomeOut.println();
}
static String getName(int index)
{
	return "read" + (index+1);
}
static boolean isBasePair(char c)
{
	return c == 'a' || c == 'A' || c == 'c' || c == 'C' || c == 'g' || c == 'G' || c == 't' || c == 'T';
}
static String getSuffix(ArrayDeque<Character> buffer, int length)
{
	StringBuilder sb = new StringBuilder();
	Iterator<Character> it = buffer.descendingIterator();
	for(int i = 0; i<length; i++)
	{
		char c = it.next();
		if(r.nextDouble() < errorRate)
		{
			sb.append(mutate(c));
		}
		else
		{
			sb.append(c);
		}
	}
	return new String(sb.reverse());
}
static char mutate(char c)
{
	char[] bases = new char[] {'a', 'A', 'c', 'C', 'g', 'G', 't', 'T'};
	for(int i = 0; i<8; i++)
	{
		if(c == bases[i])
		{
			return bases[(i + (r.nextInt(3) * 2) + 2)%8];
		}
	}
	return c;
}
static class Pair implements Comparable<Pair>
{
	long a, b;
	int i;
	int score;
	Pair(long aa, long bb, int ii)
	{
		a = aa;
		b = bb;
		i = ii;
	}
	public int compareTo(Pair o)
	{
		if(a != o.a) return Long.compare(a, o.a);
		return Long.compare(b, o.b);
	}
}
/*
 * Allows sampling from a discrete distribution given the counts of all elements in it
 */
static class Distribution
{
	int total;
	TreeMap<Integer, Integer> invFrequency;
	Random r;
	Distribution(ReadUtils.OrderedFrequencyMap<Integer> ofm, int seed)
	{
		invFrequency = new TreeMap<Integer, Integer>();
		int cFreq = 0;
		for(int x : ofm.freq.keySet())
		{
			cFreq += ofm.count(x);
			invFrequency.put(cFreq-1, x);
		}
		r = new Random(seed);
		total = cFreq;
	}
	int sample()
	{
		int cur = r.nextInt(total);
		return invFrequency.get(invFrequency.ceilingKey(cur));
	}
}
static Distribution getDistribution(String fn) throws IOException
{
	ReadUtils.OrderedFrequencyMap<Integer> ofm = ReadUtils.getLengths(fn);
	int seed = 171;
	Distribution d = new Distribution(ofm, seed);
	return d;
}

public static class ReadUtils {
	static enum FileType {
		EMPTY, FASTA, FASTQ;
	};
static void testGetLengths() throws IOException
{
	FileType[] fts = new FileType[] {FileType.FASTA, FileType.FASTQ, FileType.FASTA, FileType.FASTQ, FileType.FASTA};
	int[][] lengths = new int[][] {
		{10, 20, 30, 40, 50},
		{10, 20, 30, 40, 50},
		{4, 4, 4, 4, 3, 3, 3, 2, 2, 1},
		{4, 4, 4, 4, 3, 3, 3, 2, 2, 1},
		{100}
	};
	int numTests = fts.length;
	for(int i = 0; i<numTests; i++)
	{
		boolean result = test(fts[i], lengths[i]);
		System.out.println("TEST " + (i+1) + ": "+fts[i] + " " 
				+ Arrays.toString(lengths[i]) + " " + (result ? "PASSED" : "FAILED"));
	}
}
static boolean test(FileType ft, int[] lengths) throws IOException
{
	String fn = "test.txt";
	File f = new File(fn);
	PrintWriter out = new PrintWriter(f);
	for(int i = 0; i<lengths.length; i++)
	{
		if(ft == FileType.FASTA)
		{
			int lineLength = 80;
			out.println(">read"+(i+1));
			String read = randomRead(lengths[i], i);
			for(int j = 0; j+lineLength<lengths[i]; j+=lineLength)
			{
				out.println(read.substring(j, Math.min(lengths[i], j+lineLength)));
			}
		}
		else
		{
			out.println(">read"+(i+1));
			String read = randomRead(lengths[i], i);
			out.println(read);
			out.println("+");
			for(int j = 0; j<lengths[i]; j++) out.print("*");
			out.println();
		}
	}
	out.close();
	
	OrderedFrequencyMap<Integer> myLengths = getLengths(fn);
	OrderedFrequencyMap<Integer> trueLengths = getLengths(fn);
	
	f.delete();
	
	return myLengths.toString().equals(trueLengths.toString());
}
/*
 * Generaets a random read of a fixed length for testing purposes
 */
static String randomRead(int length, int seed)
{
	char[] bases = new char[] {'A', 'C', 'G', 'T'};
	Random r = new Random(seed);
	char[] res = new char[length];
	for(int i = 0; i<length; i++)
	{
		res[i] = bases[r.nextInt(bases.length)];
	}
	return new String(res);
}
/*
 * Gets the multi-set of read lengths from a read file
 */
static OrderedFrequencyMap<Integer> getLengths(String fn) throws IOException
{
	OrderedFrequencyMap<Integer> res = new OrderedFrequencyMap<Integer>();
	PeekableScanner  input = new PeekableScanner(new FileInputStream(new File(fn)));
	FileType ft = getFileType(input);
	if(ft == FileType.EMPTY)
	{
		return res;
	}
	while(true)
	{
		String curRead = getUnlabelledRead(input, ft);
		if(curRead == null || curRead.length() == 0)
		{
			break;
		}
		if(curRead.length() > maxLength)
		{
			continue;
		}
		res.add(curRead.length());
	}
	return res;
}
public static String getName(PeekableScanner input, FileType ft)
{
	if(!input.hasNext())
	{
		return null;
	}
	if(ft == FileType.FASTQ)
	{
		String res = "";
		for(int i = 0; i<4; i++)
		{
			if(!input.hasNext())
			{
				return null;
			}
			if(i == 0)
			{
				res = input.nextLine();
				if(res.length() != 0)
				{
					res = res.substring(1);
				}
				else
				{
					res = null;
				}
			}
			else
			{
				input.nextLine();
			}
		}
		return res;
	}
	else if(ft == FileType.FASTA)
	{
		if(!input.hasNext())
		{
			return null;
		}
		String res = input.nextLine();
		
		if(res.length() != 0)
		{
			res =  res.substring(1);
		}
		else
		{
			res = null;
		}
		
		while(input.hasNext())
		{
			String curLine = input.peekLine();
			if(curLine.length() == 0 || curLine.startsWith(">"))
			{
				break;
			}
			input.nextLine();
		}
		
		return res;
	}
	else
	{
		return null;
	}
}
/*
 * Scans the next read and returns its name and sequence
 */
static String[] getLabelledRead(PeekableScanner input, FileType ft)
{
	String[] res = new String[2];
	if(!input.hasNext())
	{
		return null;
	}
	if(ft == FileType.FASTQ)
	{
		for(int i = 0; i<4; i++)
		{
			if(!input.hasNext())
			{
				return null;
			}
			if(i == 0) res[0] = input.nextLine();
			else if(i == 1) res[1] = input.nextLine();
			else input.nextLine();
		}
		if(res[0].length() == 0)
		{
			return null;
		}
		else
		{
			res[0] = res[0].substring(1);
		}
		return res;
	}
	else if(ft == FileType.FASTA)
	{
		if(!input.hasNext())
		{
			return null;
		}
		res[0] = input.nextLine();
		
		StringBuilder sb = new StringBuilder("");
		
		while(input.hasNext())
		{
			String curLine = input.peekLine();
			if(curLine.length() == 0 || curLine.startsWith(">"))
			{
				break;
			}
			sb.append(input.nextLine());
		}
		
		res[1] = sb.toString();
		
		if(res[0].length() == 0)
		{
			return null;
		}
		else
		{
			res[0] = res[0].substring(1);
		}
		
		return res;
	}
	else
	{
		return null;
	}
}
/*
 * Scans the next read and returns the sequence associated with it
 */
static String getUnlabelledRead(PeekableScanner input, FileType ft)
{
	if(!input.hasNext())
	{
		return null;
	}
	if(ft == FileType.FASTQ)
	{
		String res = "";
		for(int i = 0; i<4; i++)
		{
			if(!input.hasNext())
			{
				return null;
			}
			if(i == 1) res = input.nextLine();
			else input.nextLine();
		}
		return res;
	}
	else if(ft == FileType.FASTA)
	{
		if(!input.hasNext())
		{
			return null;
		}
		input.nextLine();
		
		StringBuilder res = new StringBuilder("");
		
		while(input.hasNext())
		{
			String curLine = input.peekLine();
			if(curLine.length() == 0 || curLine.startsWith(">"))
			{
				break;
			}
			res.append(input.nextLine());
		}
		
		return res.toString();
	}
	else
	{
		return null;
	}
}
/*
 * Gets the filetype of a read file: fasta, fastq, or other/empty
 */
static FileType getFileType(PeekableScanner input) throws IOException
{
	if(!input.hasNext())
	{
		return FileType.EMPTY;
	}
	String line = input.peekLine();
	if(line.startsWith(">"))
	{
		return FileType.FASTA;
	}
	else if(line.startsWith("@"))
	{
		return FileType.FASTQ;
	}
	else
	{
		return FileType.EMPTY;
	}
}
/*
 * Similar to a TreeSet but keeps a map from element to frequency instead
 */
static class OrderedFrequencyMap<T>
{
	TreeMap<T, Integer> freq;
	OrderedFrequencyMap()
	{
		freq = new TreeMap<T, Integer>();
	}
	OrderedFrequencyMap(T[] data)
	{
		freq = new TreeMap<T, Integer>();
		for(T x : data)
		{
			add(x);
		}
	}
	void add(T x)
	{
		freq.put(x, freq.containsKey(x) ? (1 + freq.get(x)) : 1);
	}
	int count(T x)
	{
		return freq.containsKey(x) ? freq.get(x) : 0;
	}
	public String toString()
	{
		return freq.toString();
	}
}
/*
 * File scanner with the additional ability to peek at the next line
 */
static class PeekableScanner
{
	String lastLine;
	boolean hasLastLine;
	Scanner sc;
	PeekableScanner(InputStream is)
	{
		sc = new Scanner(is);
		hasLastLine = false;
	}
	String peekLine()
	{
		if(hasLastLine)
		{
			return lastLine;
		}
		else
		{
			lastLine = sc.nextLine();
			hasLastLine = true;
			return lastLine;
		}
	}
	String nextLine()
	{
		if(hasLastLine)
		{
			hasLastLine = false;
			String res = lastLine;
			lastLine = null;
			return res;
		}
		else
		{
			return sc.nextLine();
		}
	}
	boolean hasNext()
	{
		return hasLastLine || sc.hasNext();
	}
}
}
}
