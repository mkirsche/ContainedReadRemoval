import java.util.*;
import java.io.*;
public class HashContainment {
	static int FREQ_MINIMIZERS = 8; // Frequency will be about 1/(2^x)
	static int K = 20;
	static double CONTAINMENT_THRESHOLD = 0.85;
	static int samples = 5;
public static void main(String[] args) throws IOException
{
	String fn = "/home/mkirsche/ccs/chr22.fastq";
	if(args.length > 0)
	{
		if(args.length == 1)
		{
			System.out.println("readfilename freqminimizers k containmentthreshold seed=");
			return;
		}
		else
		{
			fn = args[0];
			FREQ_MINIMIZERS = Integer.parseInt(args[1]);
			K = Integer.parseInt(args[2]);
			CONTAINMENT_THRESHOLD = Integer.parseInt(args[3]) * 1. / 100;
			for(String s : args)
			{
				if(s.startsWith("seed="))
				{
					samples = Integer.parseInt(s.substring("seed=".length()));
					break;
				}
			}
		}
	}
	HashMap<Long, HashSet<Integer>> map = new HashMap<>();
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	
	ArrayList<Read> rs = new ArrayList<Read>();
	int count = 0;
	boolean fastq = !fn.endsWith("fasta") && !fn.endsWith("fa");
	while(input.hasNext())
	{
		count++;
		if(count%1000 == 0) System.err.println("Input " + count + " reads");
		rs.add(new Read(input.nextLine(), input.nextLine()));
		if(fastq)
		{
			input.nextLine();
			input.nextLine();
		}
	}
	Collections.sort(rs);
	int n = rs.size();
	boolean[] contained = new boolean[n];
	System.err.println("Total reads: " + n);
	Random r = new Random(50);
	for(int i = 0; i<n; i++)
	{
		if(i % 1000 == 999) System.err.println("Processed " + (i+1) + " reads");
		HashSet<Integer> check = new HashSet<Integer>();
		
		int sz = rs.get(i).ms.length;
		if(sz == 0)
		{
			contained[i] = true;
			continue;
		}
		if(samples > 0)
		{
			for(int ss = 0; ss<samples; ss++)
			{
				int idx = r.nextInt(sz);
				if(!map.containsKey(rs.get(i).ms[idx])) continue;
				for(int x : map.get(rs.get(i).ms[idx])) check.add(x);
			}
			//System.out.println(check.size());
			for(int j : check)
			{
				if(i == j) continue;
				if(rs.get(j).contains(rs.get(i)))
				{
					contained[i] = true;
					break;
				}
			}
		}
		else
		{
			for(int j = 0; j<i; j++)
			{
				if(rs.get(j).contains(rs.get(i)))
				{
					contained[i] = true;
					break;
				}
			}
		}
		if(samples > 0 && !contained[i])
		{
			for(long x : rs.get(i).ms)
			{
				if(!map.containsKey(x)) map.put(x, new HashSet<Integer>());
				map.get(x).add(i);
			}
		}
	}
	int countContained = 0;
	for(int i = 0; i<n; i++) if(contained[i]) countContained++;
	System.err.println("Number contained: " + countContained);
	String ofn = fn + ".uncontained_hash" + "." + FREQ_MINIMIZERS + "_" + K + "_" 
			+ String.format("%.2f", CONTAINMENT_THRESHOLD) + "_" + samples;
	System.out.println(ofn);
	PrintWriter out = new PrintWriter(new File(ofn));
	for(int i = 0; i<n; i++)
		if(!contained[i])
			out.println(rs.get(i).name);
	out.close();
}
static long hash(long val, long m)
{
	long x = (~val + (val << 21));
	x = x ^ (x >> 24);
	x = (x + (x<<3) + (x<<8));
	x = x ^ (x >> 14);
	x = (x + (x<<2) + (x<<4));
	x = x ^ (x >> 28);
	x = (x + (x << 31)) & ((1L<<m)-1);
	return x;
}
static int map(char c)
{
	if(c == 'A') return 0;
	else if(c == 'C') return 1;
	else if(c == 'G') return 2;
	return 3;
}
static long revComp(long x)
{
	long res = 0;
	for(int i = 0; i<K; i++)
	{
		res<<=2;
		res |= ((x&3)^3);
		x>>=2;
	}
	return res;
}
static long[] getModimizers(String s)
{
	TreeSet<Long> kmers = new TreeSet<Long>();
	int n = s.length();
	long kmer = 0;
	for(int i = 0; i<K; i++) kmer = (kmer << 2) | map(s.charAt(i));
	long kmer2 = revComp(kmer);
	long hash2 = hash(kmer2, FREQ_MINIMIZERS);
	long hash = hash(kmer, FREQ_MINIMIZERS);
	if(hash == 0) kmers.add(kmer);
	if(hash2 == 0) kmers.add(kmer2);
	for(int i = K; i<n; i++)
	{
		kmer = kmer & ((1L << (2*K - 2)) - 1);
		kmer <<= 2;
		kmer |= map(s.charAt(i));
		kmer2 = revComp(kmer);
		hash2 = hash(kmer2, FREQ_MINIMIZERS);
		hash = hash(kmer, FREQ_MINIMIZERS);
		//System.out.println(kmer+" "+hash);
		if(hash == 0) kmers.add(kmer);
		if(hash2 == 0) kmers.add(kmer2);
	}
	//System.out.println(kmers);
	long[] res = new long[kmers.size()];
	int idx = 0;
	for(long km : kmers) res[idx++] = km;
	//System.out.println(n+" "+res.length);
	return res;
}
static class Read implements Comparable<Read>
{
	String name;
	long[] ms;
	int len;
	Read(String l1, String l2)
	{
		name = l1.substring(1);
		ms = getModimizers(l2.toUpperCase());
		len = l2.length();
	}
	boolean contains(Read r)
	{
		int common = 0, n = ms.length, m = r.ms.length;
		int i = 0, j = 0;
		while(i < n && j < m)
		{
			long a = ms[i], b = r.ms[j];
			if(a < b) i++;
			else if(a > b) j++;
			else
			{
				i++;
				j++;
				common++;
			}
		}
		if(common > CONTAINMENT_THRESHOLD * m - 1e-9)
			return true;
		return false;
	}
	public int compareTo(Read o) {
		return o.len - len;
	}
}
}
