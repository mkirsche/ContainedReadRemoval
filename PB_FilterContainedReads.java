/*
 * Filter contained reads from a PacBio dataset
 * Algorithm tuned to handle error rates of 10-15%
 */
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.io.*;
public class PB_FilterContainedReads {
	
	/*
	 * The kmer length and window size for sketching the reads
	 */
	static int K1 = 8, W1 = 5;
	
	/*
	 * The proportion of kmers in a read which must be contained 
	 * in another read to imply the read being contained
	 */
	static double CONTAINMENT_THRESHOLD = 0.1;
	
	/*
	 * The sketch parmeter for making an intial pass for speeding up the algorithm
	 * A read A is thought to possibly contain another read B only if they share at least one (K2, W2) minimizer
	 */
	static int K2 = 14, W2 = 50;
	
	/*
	 * Whether or not the read file is in fastq format as opposed to fasta
	 */
	static boolean fastq = true;
	
	/*
	 * The number of threads to use
	 */
	static int NUM_THREADS = 24;
	
	/*
	 * The maximum number of reads to consider for containment
	 */
	static int LIMIT = 50;
	
	/*
	 * Input and output filenames
	 */
	static String fn, ofn;
	
	/*
	 * List of reads
	 */
	static ArrayList<Read> rs;
	
	/*
	 * A counter keeping track of how many reads have been processed so far
	 */
	static AtomicInteger processed;
		
	/*
	 * The number of reads to assign each thread at a time
	 */
	static 	int iter = 5000;
	
	/*
	 * The number of reads to process on a single thread before procesing the rest in parallel
	 */
	static int PREPROCESS = 5000;
	
	/*
	 * A set of reads containing each kmer
	 */
	static ConcurrentHashMap<Long, ConcurrentLinkedDeque<Integer>> map;
	
	/*
	 * Whether or not each read is contained
	 */
	static boolean[] contained;
	
@SuppressWarnings("resource")
public static void main(String[] args) throws Exception
{
	// Get system time for timing the program
	long startTime = System.currentTimeMillis();
	
	// Initialize filenames
	fn = "ERR2173373.fastq";
	generateOutputFilename();
	
	// Initialize input reader and global variables
	BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
	map = new ConcurrentHashMap<>();
	rs = new ArrayList<Read>();
	processed = new AtomicInteger();
	
	// List of threads
	ArrayList<MyThread> ts = new ArrayList<MyThread>();
	
	// Keep track of which thread we are adding to and what range of reads to send it
	int thread = 0;
	int lastEnd = -1;
	
	// Counter of how many reads have been input
	int countInput = 0;
	
	// Scan through reads and produce a sketch for each read
	while(true)
	{
		try 
		{
			countInput++;
			if(countInput%iter == 0)
			{
				// End of batch so send to a thread to compute minimizers
				int start = lastEnd + 1;
				int end = countInput - 1;
				lastEnd = end;
				
				// Check if max number of threads has been reached; if so, replace one with this one
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
			
			// Add read to read list
			String name = input.readLine(), read = input.readLine();
			rs.add(new Read(name, read));
			if(fastq)
			{	
				input.readLine();
				input.readLine();
			}
		} 
		catch(Exception e)
		{
			break;
		}
	}
	
	// Finish partial batch at the end
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
	
	// Wait for all threads to finish
	for (MyThread th : ts) {
	    th.join();
	}
	
	// Sort reads by descending length
	System.err.println("Sorting reads");
	Collections.sort(rs);
	
	// Initialize all reads as non-contained
	int n = rs.size();
	contained = new boolean[n];
	
	// Process a batch initially to check all future reads against
	for(int i = 0; i<PREPROCESS; i++) process(i);
	
	// Determine range of reads each thread will process
	int[] starts = new int[NUM_THREADS], ends = new int[NUM_THREADS];
	starts[0] = PREPROCESS;
	int per = Math.max(0, (n - starts[0]) / NUM_THREADS);
	for(int i = 1; i<NUM_THREADS; i++) starts[i] = starts[i-1] + per;
	for(int i = 0; i<NUM_THREADS; i++) ends[i] = i == (NUM_THREADS - 1) ? n-1 : (starts[i+1] - 1);
	
	// Start all threads
	for(int i = 0; i<NUM_THREADS; i++)
	{
		ts.set(i, new MyThread(starts[i], ends[i], 1));
		ts.get(i).start();
	}
	
	// Wait for all threads to finish before finalizing results
	for (MyThread th : ts) {
	    th.join();
	}
	
	// Output names of non-contained reads
	int count = 0;
	PrintWriter out = new PrintWriter(new File(ofn));
	for(int i = 0; i<n; i++)
		if(!contained[i])
			out.println(rs.get(i).name);
		else count++;
	
	// Output number contained and runtime
	System.out.println(count + " contained out of " + n);
	out.close();
	
	long endTime = System.currentTimeMillis();
	System.out.println("Time (ms): " + (endTime - startTime));
}
/*
 * Process read with index i and check if it is contained
 */
static void process(int i) throws Exception
{
	int done = processed.incrementAndGet();
	if(done%1000 == 0) System.err.println("Processed " + done + " reads");
	
	// Set of other reads to check
	HashSet<Integer> check = new HashSet<Integer>();
	
	if(rs.get(i) == null || rs.get(i).ms2 == null)
	{
		contained[i] = true;
		return;
	}
	int sz = rs.get(i).ms2.length;
	if(sz == 0)
	{
		contained[i] = true;
		return;
	}
	
	// Find reads with shared (k2, w2) minimizer and add them to set of candidates
	for(long x : rs.get(i).ms2)
	{
		if(!map.containsKey(x)) continue;
		for(int y : map.get(x)) check.add(y);
	}
	
	// Check if any of the candidates contain this read
	for(int j : check)
	{
		if(i == j) continue;
		if(rs.get(j).contains(rs.get(i)))
		{
			contained[i] = true;
			break;
		}
	}
	
	// If this read is not contained, let it be a candidate for future reads
	if(!contained[i])
	{
		for(long x : rs.get(i).ms2)
		{
			if(!map.containsKey(x)) map.put(x, new ConcurrentLinkedDeque<Integer>());
			if(map.get(x).size() < LIMIT) map.get(x).add(i);
		}
	}
	rs.get(i).ms2 = null;
}

/*
 * Thread for computing minimizers or checking containment
 */
static class MyThread extends Thread
{
	// Performs this phase's function (based on p) on reads in range [x, y]
	int x, y, p;
	
	// phase 0 = compute minimizers; phase 1 = check containment
	public MyThread(int i, int j, int phase)
	{
		x = i; y = j; p = phase;
	}
	public void run() {
		 try {

		      for(int i = x; i<=y; i++){
		    	  if(p == 0) rs.get(i).init();
		    	  else process(i);
		      }

		    } catch(Exception e) {

		      System.err.println("Error: " + e.getMessage());     
		    }
	}
}

/*
 * Generate output filename based on input filename and parameters
 */
static void generateOutputFilename()
{
	ofn = fn + "." + K1 + "." + K2 + "." + W1 + "." + W2 + String.format("%.2f", CONTAINMENT_THRESHOLD); 
}

/*
 * Computes the (K, W) minimizers of a string s
 */
static long[] getWindowMinimizers(String s, int K, int W)
{
	// String too short for any windows - return empty list
	if(s.length() <= K + W) return  new long[] {};
	
	// Set to hold minimizers
	HashSet<Long> kmers = new HashSet<Long>();
	
	// Current values of kmer and reverse complement
	long[] ks = new long[2];
	
	// Current kmers and hashes in the window
	long[] kWindow = new long[W<<1];
	long[] vals = new long[W<<1];
	
	// Fill initial window
	int idx = 0;
	for(int i = 0; i<K+W; i++)
	{
		ks = updateKmers(K, ks, s.charAt(i));
		if(i>=K)
		{
			for(int j = 0; j<2; j++)
			{
				kWindow[idx] = ks[j];
				vals[idx] = hash(ks[j]);
				idx++;
				if(idx == vals.length) idx = 0;
			}
		}
	}
	
	// Add minimizer from first window
	int bestIdx = 0;
	for(int i = 1; i<vals.length; i++) if(vals[i] < vals[bestIdx]) bestIdx = i;
	for(int i = 0; i<vals.length; i++) if(vals[i] == vals[bestIdx]) kmers.add(kWindow[i]);
	
	// Compute subsequent windows
	for(int c = K+W; c<s.length(); c++)
	{
		// Update window
		ks = updateKmers(K, ks, s.charAt(c));
		for(int j = 0; j<2; j++)
		{
			kWindow[idx] = ks[j];
			vals[idx] = hash(ks[j]);
			idx++;
			if(idx == vals.length) idx = 0;
		}
		
		// Add minimizer
		bestIdx = 0;
		for(int i = 1; i<vals.length; i++) if(vals[i] < vals[bestIdx]) bestIdx = i;
		for(int i = 0; i<vals.length; i++) if(vals[i] == vals[bestIdx]) kmers.add(kWindow[i]);
	}
	
	// Convert set to array
	long[] res = new long[kmers.size()];
	idx = 0;
	for(long x : kmers) res[idx++] = x;
	return res;
}

/*
 * Updates a sliding kmer and its reverse complement based on the next character
 */
static long[] updateKmers(int K, long[] kmers, char c)
{
	long first = kmers[0] & ((1L << (2*K - 2)) - 1);
	first = (first << 2) | map(c);
	long second = kmers[1] >> 2;
	second = second | ((map(c) ^ 3) << (2*K-2));
	return new long[] {first, second};
}

/*
 * Maps a nucleotide into a value in [0, 3]
 */
static int map(char c)
{
	if(c == 'A') return 0;
	else if(c == 'C') return 1;
	else if(c == 'G') return 2;
	else return 3;
}

/*
 * Hashes a values into a pseudorandom 32-bit integer for min hash
 */
static long hash(long val)
{
	long m = 32;
	long x = (~val + (val << 21));
	x = x ^ (x >> 24);
	x = (x + (x<<3) + (x<<8));
	x = x ^ (x >> 14);
	x = (x + (x<<2) + (x<<4));
	x = x ^ (x >> 28);
	x = (x + (x << 31)) & ((1L<<m)-1);
	return x;
}
/*
 * A class representing a genomic read and its sketch(es)
 */
static class Read implements Comparable<Read>
{
	// The name of the read
	String name;
	
	// The primary sketch
	long[] ms;
	
	// The secondary sketch used to find candidates
	long[] ms2;
	
	// The length of the read (bp)
	int len;
	
	// The content of the read, which is removed once sketches are built
	String line;
	
	// Takes as input two lines of a fasta/fastq file
	Read(String l1, String l2)
	{
		name = l1.substring(1);
		line = l2;
		len = l2.length();
	}
	
	// Initialize read by building sketches
	void init()
	{
		line = line.toUpperCase();
		ms = getWindowMinimizers(line, K1, W1);
		ms2 = getWindowMinimizers(line, K2, W2);
		line = null;
	}
	
	// Check if this read contains another read based on shared proportion of primary sketch
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
	
	// When sorting, sort by descending read length
	public int compareTo(Read o) {
		return o.len - len;
	}
}
}
