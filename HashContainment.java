import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.io.*;
public class HashContainment {
	static int FREQ_MINIMIZERS = 8; // Frequency will be about 1/(2^x)
	static int K = 20;
	static double CONTAINMENT_THRESHOLD = 0.85;
	static int samples = -3;
	static int LIMIT = 0;
	static boolean fnOnly = false; // Whether or not to just output the file name and exit
	
	static Random r;
	static boolean[] contained;
	static ArrayList<Read> rs;
	static ConcurrentHashMap<Long, ConcurrentLinkedDeque<Integer>> map;
	static int NUM_THREADS = 1;
	static int PREPROCESS = 5000;
	static AtomicInteger processed;
public static void main(String[] args) throws Exception
{
	long startTime = System.currentTimeMillis();
	processed = new AtomicInteger();
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
			for(String s : args)
			{
				if(s.startsWith("limit="))
				{
					LIMIT = Integer.parseInt(s.substring("limit=".length()));
					break;
				}
			}
			for(String s : args)
			{
				if(s.equals("--fnOnly"))
				{
					fnOnly = true;
					break;
				}
			}
			for(String s : args)
			{
				if(s.startsWith("threads="))
				{
					NUM_THREADS = Integer.parseInt(s.substring("threads=".length()));
					break;
				}
			}
			for(String s : args)
			{
				if(s.startsWith("preprocess="))
				{
					PREPROCESS = Integer.parseInt(s.substring("preprocess=".length()));
					break;
				}
			}
		}
	}
	String ofn = fn + ".uncontained_hash" + "." + FREQ_MINIMIZERS + "_" + K + "_" 
			+ String.format("%.2f", CONTAINMENT_THRESHOLD) + "_" + samples;
	if(fnOnly)
	{
		System.out.println(ofn);
		return;
	}
	map = new ConcurrentHashMap<>();
	//Scanner input = new Scanner(new FileInputStream(new File(fn)));
	BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
	int countLines = 0;
	while(true)
	{
		try {
			String ss = input.readLine();
			if(ss == null) break;
			countLines ++;
		} catch(Exception e) {
			break;
		}
	}
	input = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
	rs = new ArrayList<Read>();
	int count = 0;
	boolean fastq = !fn.endsWith("fasta") && !fn.endsWith("fa");
	countLines /= (fastq ? 4 : 2);
	System.err.println("Reads: " + countLines);
	contained = new boolean[countLines];
	//while(input.hasNext())
	ArrayList<MyThread> ts = new ArrayList<MyThread>();
	int iter = 5000;
	r = new Random(50);
	int thread = 0;
	int lastEnd = -1;
	while(true)
	{
		try {
		count++;
		if(count == PREPROCESS || count > PREPROCESS && (count - PREPROCESS)%iter == 0)
		{
			int start = lastEnd + 1;
			int end = count - 1;
			lastEnd = end;
			if(ts.size() < NUM_THREADS)
			{
				ts.add(new MyThread(start, end, 0));
				ts.get(ts.size() - 1).start();
			}
			else
			{
				int idx = thread%ts.size();
				ts.get(idx).join();
				ts.set(idx, new MyThread(start, end, 0));
				ts.get(idx).start();
			}
			thread++;
			System.err.println("Input " + count + " reads (threads = " + ts.size() + ")");
		}
		rs.add(new Read(input.readLine(), input.readLine()));
		if(fastq)
		{	
			input.readLine();
			input.readLine();
		}
		} catch(Exception e) {
			break;
		}
	}
	int start = lastEnd + 1;
	int end = countLines - 1;
	if(ts.size() < NUM_THREADS)
	{
		ts.add(new MyThread(start, end, 0));
		ts.get(ts.size() - 1).start();
	}
	else
	{
		int idx = thread%ts.size();
		ts.get(idx).join();
		ts.set(idx, new MyThread(start, end, 0));
		ts.get(idx).start();
	}
	for (MyThread th : ts) {
	    th.join();
	}
	Collections.sort(rs);
	int n = rs.size();
	System.err.println("Total reads: " + n);
	for(int i = 0; i<PREPROCESS; i++) process(i);
	int[] starts = new int[NUM_THREADS], ends = new int[NUM_THREADS];
	starts[0] = PREPROCESS;
	int per = Math.max(0, (n - starts[0]) / NUM_THREADS);
	for(int i = 1; i<NUM_THREADS; i++) starts[i] = starts[i-1] + per;
	for(int i = 0; i<NUM_THREADS; i++) ends[i] = i == (NUM_THREADS - 1) ? n-1 : (starts[i+1] - 1);
	for(int i = 0; i<NUM_THREADS; i++)
	{
		ts.set(i, new MyThread(starts[i], ends[i], 1));
		ts.get(i).start();
	}
	for (MyThread th : ts) {
	    th.join();
	}
//	for(int i = 0; i<n; i++)
//	{
//		process(i);
//	}
	int countContained = 0;
	for(int i = 0; i<n; i++) if(contained[i]) countContained++;
	System.err.println("Number contained: " + countContained);
	System.out.println(ofn);
	PrintWriter out = new PrintWriter(new File(ofn));
	for(int i = 0; i<n; i++)
		if(!contained[i])
			out.println(rs.get(i).name);
	out.close();
	long endTime = System.currentTimeMillis();
	System.err.println("Time (ms): " + (endTime - startTime));
}
static void process(int i) throws Exception
{
	int done = processed.incrementAndGet();
	if(done%1000 == 0) System.err.println("Processed " + done + " reads");
	//rs.get(i).init();
	HashSet<Integer> check = new HashSet<Integer>();
	
	int sz = rs.get(i).ms.length;
	if(sz == 0)
	{
		contained[i] = true;
		return;
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
	else if(samples < 0)
	{
		HashMap<Integer, Integer> checkmap = new HashMap<Integer, Integer>();
		for(int ss = 0; ss<-samples; ss++)
		{
			int idx = r.nextInt(sz);
			if(!map.containsKey(rs.get(i).ms[idx])) continue;
			for(int x : map.get(rs.get(i).ms[idx]))
			{
				if(checkmap.containsKey(x)) checkmap.put(x, checkmap.get(x)+1);
				else checkmap.put(x, 1);
			}
		}
		for(int x : checkmap.keySet()) if(checkmap.get(x) > 1) check.add(x);
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
	if(samples != 0 && !contained[i])
	{
		for(long x : rs.get(i).ms)
		{
			if(!map.containsKey(x)) map.put(x, new ConcurrentLinkedDeque<Integer>());
			if(LIMIT == 0 || map.get(x).size() < LIMIT) map.get(x).add(i);
		}
	}
}
static class MyThread extends Thread
{
	int x, y, p;
	public MyThread(int i, int j, int phase)
	{
		x = i; y = j; p = phase;
	}
	public void run() {
		 try {

		      for(int i = x; i<=y; i++){
		    	  if(p == 0) rs.get(i).init();
		    	  if(p == 1) process(i);
		      }

		    } catch(Exception e) {

		      System.err.println("Error: " + e.getMessage());      
		    }
	}
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
	HashSet<Long> kmers = new HashSet<Long>();
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
	Arrays.sort(res);
	//System.out.println(n+" "+res.length);
	return res;
}
static class Read implements Comparable<Read>
{
	String name;
	long[] ms;
	int len;
	String line;
	Read(String l1, String l2)
	{
		name = l1.substring(1);
		line = l2;
		//ms = getModimizers(l2.toUpperCase());
		len = l2.length();
	}
	void init()
	{
		ms = getModimizers(line.toUpperCase());
		line = "";
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
