/*
 * Using hashing to determine which reads in a dataset are contained within
 * other reads and remove those reads
 */
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.io.*;
public class HashContainment {
	/*
	 * Parameters:
	 * 
	 * FREQ_MINIMIZERS: The sampling frequency when creating a sketch of each read.
	 *   If this value is x, approximately 1 in 2^x kmers will be stored for mod hashing,
	 *   and the window size will be 2^x for min hashing.
	 *   
	 * K: The number of base pairs in each kmer
	 * 
	 * CONTAINMENT_THRESHOLD: A modified Jacard distance similarity metric is used to determine
	 *   if a read r is contained in another read r': If S is the sketch of r, and S' is the
	 *   sketch of r', then similarity(r in r') = |S intersect S'| / |S|, and r is considered
	 *   contained in r' if similarity(r in r') is at least this threshold  
	 *   
	 * SAMPLES: Checking all pairs of reads is often too slow.  In order to avoid this,
	 *   given a read r, in order to find a read that contains r, |SAMPLES| kmers are sampled
	 *   from the sketch of r, and the remaining reads are filtered according to a scheme which
	 *   depends on the value of SAMPLES
	 *     SAMPLES > 0: Only reads which contain at least one kmer in the sample
	 *     SAMPLES < 0: Only reads which contain at least two kmers in the sample
	 *     SAMPLES = 0: No filter
	 * 
	 * LIMIT: In the sampling schemes above, further filter the results by keeping only a set
	 *   of at most LIMIT reads for each kmer, and rather than just containing a kmer, reads
	 *   must be in that kmer's set
	 *   
	 * NUM_THREADS: The maximum number of threads to execute at any given time
	 * 
	 * PREPROCESS:  When running on multiple threads, there is no guarantee that a containing
	 * read has been processed because processing on multiple threads causes the reads to be
	 * considered out of order.  This parameter specifies a number of reads to process before
	 * the program forks into multiple threads to guarantee that at least the longest reads will
	 * have already been processed
	 * 
	 * HASH_TYPE: Which hashing scheme to use - 0 is mod-hash and 1 is min-hash
	 * 
	 * READ_TYPE: Auto-tune parameters to a specific read type.  Options are
	 *   ccs, pacbio, and nanopore
	 *   
	 * fn: The filename of the input file containing the reads to be analyzed
	 */
	static int FREQ_MINIMIZERS = 8;
	static int K = 20;
	static double CONTAINMENT_THRESHOLD = 0.85;
	static int SAMPLES = -3;
	static int LIMIT = 0;
	static int NUM_THREADS = 8;
	static int PREPROCESS = 5000;
	static int HASH_TYPE = 0;
	static String READ_TYPE = "";
	static String fn = "/home/mkirsche/ccs/chr22.fastq";
	
	// Set this flag to output only the name of the output file produced given the parameters
	static boolean fnOnly = false;
	
	// The filename of the output file
	static String ofn;
	
	static Random r;
	
	// A list of all reads
	static ArrayList<Read> rs;
	
	// Whether or not each read is contained in some other read
	static boolean[] contained;
	
	//Keeps a set of reads containing each kmer
	static ConcurrentHashMap<Long, ConcurrentLinkedDeque<Integer>> map;
	
	// A counter keeping track of how many reads have been processed so far
	static AtomicInteger processed;
	
	// The number of reads to assign each thread at a time
	static 	int iter = 5000;
	
public static void main(String[] args) throws Exception
{
	long startTime = System.currentTimeMillis();
	processed = new AtomicInteger();
	
	parseArgs(args);
	
	generateOutputFilename();
	
	if(fnOnly)
	{
		System.out.println(ofn);
		return;
	}
	
	// Initialize data structures for holding read information
	map = new ConcurrentHashMap<>();
	rs = new ArrayList<Read>();
	
	// Check which format the reads are in - currently based on file name
	boolean fastq = !fn.endsWith("fasta") && !fn.endsWith("fa");
	
	// Initialize threads and random number generator
	ArrayList<MyThread> ts = new ArrayList<MyThread>();
	r = new Random(50);
	int thread = 0;
	int lastEnd = -1;
	
	// Scan through reads and produce a sketch for each read
	int countInput = 0;
	BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
	while(true)
	{
		try {
		countInput++;
		if(countInput%iter == 0)
		{
			int start = lastEnd + 1;
			int end = countInput - 1;
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
			System.err.println("Input " + countInput + " reads (threads = " + ts.size() + ")");
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
	int end = rs.size() - 1;
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
	contained = new boolean[n];
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
	
	// Wait for all threads to finish before finalizing results
	for (MyThread th : ts) {
	    th.join();
	}
	
	// Count how many reads were contained
	int countContained = 0;
	for(int i = 0; i<n; i++) if(contained[i]) countContained++;
	System.err.println("Number contained: " + countContained);
	System.out.println(ofn);
	
	// Output non-contained readnames to a file
	PrintWriter out = new PrintWriter(new File(ofn));
	for(int i = 0; i<n; i++)
		if(!contained[i])
			out.println(rs.get(i).name);
	out.close();
	
	// Output the total runtime of the program
	long endTime = System.currentTimeMillis();
	System.err.println("Time (ms): " + (endTime - startTime));
}
static void generateOutputFilename()
{
	if(READ_TYPE.length() > 0) ofn = fn + ".uncontained_hash" + "." + READ_TYPE;
	else ofn = fn + ".uncontained_hash" + "." + FREQ_MINIMIZERS + "_" 
	+ K + "_" + String.format("%.2f", CONTAINMENT_THRESHOLD) + "_" 
			+ SAMPLES + "_" + LIMIT + "_" + HASH_TYPE;
}
/*
 * Parse command line arguments
 */
static void parseArgs(String[] args)
{
	if(args.length > 0)
	{
		if(args.length == 1)
		{
			// Invalid number of arguments - output usage message
			System.out.println("Usage:\n"
					+ "java HashContainment readfilename freqminimizers k containmentthreshold");
			System.out.println("Optional parameters:\n"
					+ "seed= limit= --fnOnly threads= preprocess= hashtype= readtype=");
			return;
		}
		fn = args[0];
		for(String s : args)
		{
			if(s.startsWith("readtype="))
			{
				READ_TYPE = s.substring("readtype=".length());
				break;
			}
		}
		// TODO - determine best parameters for each type of reads
		if(READ_TYPE.equals("ccs"))
		{
			
		}
		else if(READ_TYPE.equals("pacbio"))
		{
			
		}
		else if(READ_TYPE.equals("nanopore"))
		{
			
		}
		else
		{
			FREQ_MINIMIZERS = Integer.parseInt(args[1]);
			K = Integer.parseInt(args[2]);
			CONTAINMENT_THRESHOLD = Integer.parseInt(args[3]) * 1. / 100;
			for(String s : args)
			{
				if(s.startsWith("seed="))
				{
					SAMPLES = Integer.parseInt(s.substring("seed=".length()));
					break;
				}
				if(s.startsWith("limit="))
				{
					LIMIT = Integer.parseInt(s.substring("limit=".length()));
					break;
				}
				if(s.equals("--fnOnly"))
				{
					fnOnly = true;
				}
				if(s.startsWith("threads="))
				{
					NUM_THREADS = Integer.parseInt(s.substring("threads=".length()));
				}
				if(s.startsWith("preprocess="))
				{
					PREPROCESS = Integer.parseInt(s.substring("preprocess=".length()));
				}
				if(s.startsWith("hashtype="))
				{
					HASH_TYPE = Integer.parseInt(s.substring("hashtype=".length()));
				}
			}
		}
	}
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
	if(SAMPLES > 0)
	{
		for(int ss = 0; ss<SAMPLES; ss++)
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
	else if(SAMPLES < 0)
	{
		HashMap<Integer, Integer> checkmap = new HashMap<Integer, Integer>();
		for(int ss = 0; ss<-SAMPLES; ss++)
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
	if(SAMPLES != 0 && !contained[i])
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
	long[] res = new long[kmers.size()];
	int idx = 0;
	for(long km : kmers) res[idx++] = km;
	Arrays.sort(res);
	return res;
}
static long[] getMinimizers(String s)
{
	HashSet<Long> kmers = new HashSet<Long>();
	int n = s.length();
	long kmer = 0;
	for(int i = 0; i<K; i++) kmer = (kmer << 2) | map(s.charAt(i));
	long kmer2 = revComp(kmer);
	long mod = (1L<<50) - 1;
	long hash2 = hash(kmer2, mod);
	long hash = hash(kmer, mod);
	int window = 1 << FREQ_MINIMIZERS;
	MinQueue mq = new MinQueue();
	mq.add(hash, kmer);
	mq.add(hash2,  kmer2);
	for(int i = K; i<n; i++)
	{
		kmer = kmer & ((1L << (2*K - 2)) - 1);
		kmer <<= 2;
		kmer |= map(s.charAt(i));
		kmer2 = revComp(kmer);
		hash2 = hash(kmer2, FREQ_MINIMIZERS);
		hash = hash(kmer, FREQ_MINIMIZERS);
		mq.add(hash, kmer);
		mq.add(hash2, kmer2);
		if(mq.size > 2*window)
		{
			mq.remove();
			mq.remove();
			long minHash = mq.minIndex();
			kmers.add(minHash);
		}
	}
	long[] res = new long[kmers.size()];
	int idx = 0;
	for(long km : kmers) res[idx++] = km;
	Arrays.sort(res);
	return res;
}
/*
 * Data structure which supports the following operations
 * Add - append a value to the end of the queue
 * Remove - remove the first value from the queue
 * Min - return the smallest value in the queue
 * Also, each value can have an index associated with it
 */
static class MinQueue
{
	int size;
	ArrayDeque<Long> d;
	ArrayDeque<Long> maxIndex;
	Queue<Long> q;
	public MinQueue()
	{
		d = new ArrayDeque<Long>();
		maxIndex = new ArrayDeque<Long>();
		q = new LinkedList<Long>();
		size = 0;
	}
	long min()
	{
		return d.getFirst();
	}
	// Returns index associated with min value
	long minIndex()
	{
		return maxIndex.getFirst();
	}
	void add(long x, long idx)
	{
		size++;
		q.add(x);
		while(!d.isEmpty() && d.getLast() > x)
		{
			d.removeLast();
			maxIndex.removeLast();
		}
		d.add(x);
		maxIndex.add(idx);
	}
	void remove()
	{
		if(q.isEmpty()) return;
		size--;
		q.poll();
		if(!d.isEmpty() && d.getFirst() == q.peek())
		{
			d.poll();
			maxIndex.poll();
		}
	}
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
		len = l2.length();
	}
	void init()
	{
		line = line.toUpperCase();
		ms = HASH_TYPE == 0 ? getModimizers(line) : getMinimizers(line);
		line = null;
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
