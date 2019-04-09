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
	static Random r;
	static String outRefFn;
public static void main(String[] args)  throws IOException
{
	genomeFn = "/home/mkirsche/2018_08_Crossdel/genome.fa";
	readOfn = "simulatedreads.fa";
	scoreOfn = "simulatedscores.txt";
	outRefFn = "simulatedgenome.txt";
	
	int argParseErrorKey = parseArgs(args);
	if(argParseErrorKey != 0)
	{
		return;
	}
	
	Scanner input = new Scanner(new FileInputStream(new File(genomeFn)));
	PrintWriter readOut = new PrintWriter(new File(readOfn));
	//readOut = new PrintWriter(System.out);
	PrintWriter scoreOut = new PrintWriter(new File(scoreOfn));
	PrintWriter genomeOut = new PrintWriter(new File(outRefFn));
	
	nd = new NormalDistribution(mean, stdev);
	r = new Random();
	
	ArrayList<Pair> coords = simulatePositions();
	System.err.println("Simulated " + coords.size() + " reads");
	
	System.err.println("Printing reads");
	printReads(input, coords, readOut, genomeOut);
	
	System.err.println("Computing scores");
	printContainmentScores(coords, scoreOut);
	
	genomeOut.close();
	scoreOut.close();
	readOut.close();
}
static String meanKey = "mean";
static String stdevKey = "stdev";
static String genomeKey = "genomefile";
static String readKey = "readfile";
static String scoreKey = "scorefile";
static String coverageKey = "coverage";
static String maxLengthKey = "maxlength";
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
	System.out.println("  " + meanKey + ": average read length");
	System.out.println("  " + stdevKey + ": standard deviation of read length");
	System.out.println("  " + coverageKey + ": coverage to simulate");
	System.out.println("  " + maxLengthKey + ": maximum prefix of genome to simulate from");
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
static Pair sampleRead(int i)
{
	int length = (int)(nd.sample() + 0.5);
	long start = r.nextLong() % (maxLength - length);
	if(start < 0) start += maxLength - length;
	long end = start + length;
	return new Pair(start, end, i);
}
static void printContainmentScores(ArrayList<Pair> coords, PrintWriter out)
{
	Collections.sort(coords);
	int n = coords.size();
	
	int radius = 100;
	
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
	int bufferLength = 50000;
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
}
