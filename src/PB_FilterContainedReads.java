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
	 * For lengths shorter than this threshold, the containment threshold is cut in half,
	 * prioritizing removing shorter reads with the filter.
	 */
	static int LENGTH_FILTER = 12000;
	
	/*
	 * The kmer length and window size for sketching the reads
	 */
	static int K1 = 12, W1 = 5;
	
	/*
	 * The proportion of kmers in a read which must be contained 
	 * in another read to imply the read being contained
	 */
	static double CONTAINMENT_THRESHOLD = 0.25;
	static double CONTAINMENT_THRESHOLD2 = 0.25;
	static boolean setCt2 = false;
	static double REPEAT_THRESHOLD = 0.25;
	
	/*
	 * The sketch parameter for making an initial pass for speeding up the algorithm
	 * A read A is thought to possibly contain another read B only if they share at least one (K2, W2) minimizer
	 */
	static int K2 = 18, W2 = 50;
	
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
	 * The number of reads to process on a single thread before processing the rest in parallel
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
	
	/*
	 * If true, only output the output filename and exit
	 */
	static boolean fnOnly = false;
	
	/*
	 * Output for debugging and analysis
	 */
	static String debugFn = "containment_debug.txt";
	static String[] debugLines;
	
	/*
	 * Possible filtering methods
	 */
	static enum FilterMethod {
			THRESHOLD, ADAPTIVE_THRESHOLD, RECTANGLE
	};
	
	static FilterMethod method = FilterMethod.ADAPTIVE_THRESHOLD;
	
@SuppressWarnings("resource")
public static void main(String[] args) throws Exception
{
	// Get system time for timing the program
	long startTime = System.currentTimeMillis();
	
	// Initialize filenames
	fn = "ERR2173373.fastq";
	
	ofn = "";
	
	parseArgs(args);
	
	if(ofn.length() == 0)
	{
		generateOutputFilename();
	}
	
	if(fnOnly)
	{
		System.out.println(ofn);
		return;
	}
	
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
			if(name.startsWith(">"))
			{
				fastq = false;
			}	
			rs.add(new Read(name.split(" ")[0], read));
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
	debugLines = new String[n];
	
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
	System.err.println(count + " contained out of " + n);
	out.close();
	
	long endTime = System.currentTimeMillis();
	System.err.println("Time (ms): " + (endTime - startTime));

	PrintWriter debugOut = new PrintWriter(new File("debug.txt"));
	for(boolean b : contained) debugOut.println(b);
	debugOut.close();
	
	debugOut = new PrintWriter(new File(debugFn));
	for(String s : debugLines)
	{
		debugOut.println(s);
	}
	debugOut.close();
}
/*
 * Parse command line arguments
 */
static void parseArgs(String[] args)
{
	if(args.length == 0) return;
	if(args[0].equals("help"))
	{
		System.out.println("Usage:\n");
		System.out.println("  java PB_FilterContainedReads <readfile>\n");
		System.out.println("Optional arguments:");
		System.out.println("  nt=[num_threads (int)]");
		System.out.println("  k1=[k1 (int)]");
		System.out.println("  w1=[w1 (int)]");
		System.out.println("  k2=[k2 (int)]");
		System.out.println("  w2=[w2 (int)]");
		System.out.println("  ct=[containment_threshold (float)]");
		System.out.println("  ct2=[containment_threshold2 (float)]");
		System.out.println("  lf=[length_filter (int)]");
		System.out.println("  ofn=[output filename (string)]");
		System.out.println("  dfn=[debug filename (string)]");
		System.out.println("  rt=[repeat threshold (float)]");
		System.out.println("  method=[filtering method (THRESHOLD, ADAPTIVE_THRESHOLD, or RECTANGLE)]");
		System.out.println("  fnonly");
	}
	fn = args[0];
	for(int i = 1; i<args.length; i++)
	{
		if(args[i].startsWith("nt="))
		{
			NUM_THREADS = Integer.parseInt(args[i].substring(1 + args[i].indexOf('=')));
		}
		else if(args[i].startsWith("k1="))
		{
			K1 = Integer.parseInt(args[i].substring(1 + args[i].indexOf('=')));
		}
		else if(args[i].startsWith("w1="))
		{
			W1 = Integer.parseInt(args[i].substring(1 + args[i].indexOf('=')));
		}
		else if(args[i].startsWith("k2="))
		{
			K2 = Integer.parseInt(args[i].substring(1 + args[i].indexOf('=')));
		}
		else if(args[i].startsWith("w2="))
		{
			W2 = Integer.parseInt(args[i].substring(1 + args[i].indexOf('=')));
		}
		else if(args[i].startsWith("ct="))
		{
			CONTAINMENT_THRESHOLD = Double.parseDouble(args[i].substring(1 + args[i].indexOf('=')));
			if(CONTAINMENT_THRESHOLD > 1) CONTAINMENT_THRESHOLD *= 0.01; // handle ct given as percent
		}
		else if(args[i].startsWith("ct2="))
		{
			setCt2 = true;
			CONTAINMENT_THRESHOLD2 = Double.parseDouble(args[i].substring(1 + args[i].indexOf('=')));
			if(CONTAINMENT_THRESHOLD2 > 1) CONTAINMENT_THRESHOLD2 *= 0.01; // handle ct given as percent
		}
		else if(args[i].startsWith("rt="))
		{
			REPEAT_THRESHOLD = Double.parseDouble(args[i].substring(1 + args[i].indexOf('=')));
			if(REPEAT_THRESHOLD > 1) REPEAT_THRESHOLD *= 0.01; // handle ct given as percent
		}
		else if(args[i].startsWith("ofn="))
		{
			ofn = args[i].substring(1 + args[i].indexOf('='));
		}
		else if(args[i].startsWith("dfn="))
		{
			debugFn = args[i].substring(1 + args[i].indexOf('='));
		}
		else if(args[i].equals("fnonly"))
		{
			fnOnly = true;
		}
		else if(args[i].startsWith("method="))
		{
			String type = args[i].substring(1+args[i].indexOf('=')).toUpperCase();
			if(type.equals("THRESHOLD"))
			{
				method = FilterMethod.THRESHOLD;
			}
			else if(type.equals("ADAPTIVE_THRESHOLD"))
			{
				method = FilterMethod.ADAPTIVE_THRESHOLD;
			}
			else if(type.equals("RECTANGLE"))
			{
				method = FilterMethod.RECTANGLE;
			}
		}
		else if(args[i].startsWith("lf="))
		{
			LENGTH_FILTER = Integer.parseInt(args[i].substring(1 + args[i].indexOf('=')));
		}
	}
	if(!setCt2)
	{
		CONTAINMENT_THRESHOLD2 = CONTAINMENT_THRESHOLD;
	}
}
/*
 * Process read with index i and check if it is contained
 */
static void process(int i) throws Exception
{
	int done = processed.incrementAndGet();
	if(done%1000 == 0) System.err.println("Processed " + done + " reads");
	
	debugLines[i] = rs.get(i).name.split(" ")[0] + " 1.0 " + " 1.0 " + " 1.0 " + rs.get(i).len;
	
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
	
	if(rs.get(i).repeats && rs.get(i).len > LENGTH_FILTER)
	{
		debugLines[i] = rs.get(i).name + " 0.0 " + " 0.0 " + " 0.0 " + rs.get(i).len;
		contained[i] = false;
		return;
	}
	
	// Find reads with shared (k2, w2) minimizer and add them to set of candidates
	for(long x : rs.get(i).ms2)
	{
		if(rs.get(i).repeats) break;
		if(!map.containsKey(x)) continue;
		if(map.get(x).size() == LIMIT)
		{
			// Kmer seen too many times
			//continue;
		}
		for(int y : map.get(x))
		{
			check.add(y);
		}
	}
	
	double maxContainment = 0.0;
	double[] maxEdges = new double[] {0, 0};
	int bestIdx = -1;
	boolean foundContained = false;
	
	// Check if any of the candidates contain this read
	for(int j : check)
	{
		if(i == j) continue;
		double[] score =  rs.get(j).containmentScore(rs.get(i));
		if(score[0] > maxContainment)
		{
			if(!foundContained)
			{
				bestIdx = j;
				maxContainment = score[0];
				maxEdges[0] = score[1];
				maxEdges[1] = score[2];
			}
		}
		if(rs.get(j).contains(rs.get(i)))
		{
			contained[i] = true;
			if(!foundContained || score[0] > maxContainment)
			{
				foundContained = true;
				bestIdx = j;
				maxContainment = score[0];
				maxEdges[0] = score[1];
				maxEdges[1] = score[2];
			}
			//break;
		}
	}
	
	debugLines[i] = rs.get(i).name + " " + maxContainment + " " + maxEdges[0] + " " + maxEdges[1] + " " + rs.get(i).len;
	if(bestIdx != -1)
	{
		debugLines[i] += " " + rs.get(bestIdx).name;
	}
	
	// If this read is not contained, let it be a candidate for future reads
	//if(!contained[i])
	{
		for(long x : rs.get(i).ms2)
		{
			if(!map.containsKey(x)) map.put(x, new ConcurrentLinkedDeque<Integer>());
			if(map.get(x).size() < LIMIT)
			{
				map.get(x).add(i);
			}
		}
	}
	rs.get(i).ms2 = null;
}

static <T> void increment(TreeMap<T, Integer> map, T key)
{
	map.put(key, 1 + (map.containsKey(key) ? map.get(key) : 0));
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
		      e.printStackTrace();
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
static long[][] getWindowMinimizers(String s, int K, int W, int endThreshold)
{
	// String too short for any windows - return empty list
	if(s.length() <= K + W) return  new long[][] {{},{}, {}};
	
	// Set to hold minimizers
	TreeMap<Long, Integer> kmers = new TreeMap<Long, Integer>();
	TreeMap<Long, Integer> startKmers = new TreeMap<Long, Integer>();
	TreeMap<Long, Integer> endKmers = new TreeMap<Long, Integer>();
		
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
	for(int i = 0; i<vals.length; i++)
	{
		if(vals[i] == vals[bestIdx])
		{
			increment(kmers, kWindow[i]);
			if(endThreshold > 0)
			{
				increment(startKmers, kWindow[i]);
			}
		}
	}
	
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
		for(int i = 0; i<vals.length; i++) if(vals[i] == vals[bestIdx])
		{
			if(c - (K + W - 1) < endThreshold)
			{
				increment(startKmers, kWindow[i]);
			}
			if(s.length() - c - 1 < endThreshold)
			{
				increment(endKmers, kWindow[i]);
			}
			increment(kmers, kWindow[i]);
		}
	}
		
	// Convert set to array
	long[] res = new long[kmers.size()*2];
	idx = 0;
	for(long x : kmers.keySet())
	{
		res[idx] = x;
		res[idx+1] = kmers.get(x);
		idx += 2;
	}
		
	long[] start = new long[startKmers.size()*2];
	long[] end = new long[endKmers.size()*2];
	idx = 0;
	for(long x : startKmers.keySet())
	{
		start[idx] = x;
		start[idx+1] = startKmers.get(x);
		idx += 2;
	}
	idx = 0;
	for(long x : endKmers.keySet())
	{
		end[idx] = x;
		end[idx+1] = endKmers.get(x);
		idx += 2;
	}
	
	return new long[][] {res, start, end};
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

static long[] everyOther(long[] a)
{
	int n = a.length/2;
	long[] res = new long[n];
	for(int i = 0; i<n; i++) res[i] = a[2*i];
	return res;
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
	
	// The sketches around the ends of the read
	long[] starts;
	long[] ends;
	long[] starts2;
	long[] ends2;
	
	// The length of the read (bp)
	int len;
	
	boolean repeats;
	
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
		long[][] k1w1Mini = getWindowMinimizers(line, K1, W1, 500);	
		ms = k1w1Mini[0];
		int repeatCount = 0, totCount = 0;
		for(int i = 1; i<ms.length; i+=2)
		{
			totCount += ms[i];
			if(ms[i] > 10)
			{
				repeatCount += ms[i];
				ms[i] = 1;
			}
		}
		if(repeatCount > REPEAT_THRESHOLD * totCount)
		{
			repeats = true;
		}
		if(repeatCount > REPEAT_THRESHOLD * totCount)
		{
			System.out.println("Repeat read: " +name+" with repetitiveness "+1.0 * repeatCount / totCount + " (needed " + REPEAT_THRESHOLD + ")");
		}
		long[][] k1Mini = getWindowMinimizers(line, K1, 1, 100);
		starts2 = k1Mini[1];
		ends2 = k1Mini[2];
		starts = k1w1Mini[1];
		ends = k1w1Mini[2];
		long[][] k2w2Mini = getWindowMinimizers(line, K2, W2, 500);
		ms2 = everyOther(k2w2Mini[0]);
		line = null;
	}
	
	// Check if this read contains another read based on shared proportion of primary sketch
	boolean contains(Read r)
	{
		double[] containmentScores = this.containmentScore(r);
		double[] thresholds = new double[] {CONTAINMENT_THRESHOLD, CONTAINMENT_THRESHOLD2, CONTAINMENT_THRESHOLD2};
		if(method == FilterMethod.THRESHOLD)
		{
			for(int i = 0; i<containmentScores.length; i++)
			{
				if(containmentScores[i] < thresholds[i] - 1e-9)
				{
					return false;
				}
			}
			return true;
		}
		else if(method == FilterMethod.ADAPTIVE_THRESHOLD)
		{
			for(int i = 0; i<containmentScores.length; i++)
			{
				if(containmentScores[i] < thresholds[i] - 1e-9 && r.len >= LENGTH_FILTER)
				{
					return false;
				}
				else if(containmentScores[i] < thresholds[i] * .5 - 1e-9)
				{
					return false;
				}
			}
			return true;
		}
		else if(method == FilterMethod.RECTANGLE)
		{
			if(r.len < LENGTH_FILTER)
			{
				return true;
			}
			boolean lessThanDouble = true;
			for(int i = 0; i<3; i++)
			{
				if(containmentScores[i] < thresholds[i] - 1e-9)
				{
					return false;
				}
				else
				{
					lessThanDouble &= containmentScores[i] < 2 * thresholds[i] - 1e-9;
				}
			}
			if(lessThanDouble)
			{
				// Check for off-by-one in ends
				for(int i = 1; i<3; i++)
				{
					long[] kmers = i == 1 ? r.starts2 : r.ends2;
					double cur = compareErrorKmers(kmers, ms, K1);
					System.out.println(cur);
					if(cur < thresholds[i] * 2 - 1e-9)
					{
						return false;
					}
				}
				return true;
			}
			else
			{
				return true;
			}
		}
		else
		{
			return false;
		}
	}
	
	static double compareErrorKmers(long[] error, long[] database, int k)
	{
		ArrayList<IndexedKmer> errored = new ArrayList<IndexedKmer>();
		for(int index = 0; index < error.length; index+=2)
		{
			long e = error[index];
			int count = (int)error[index+1];
			for(int i = 0; i<k; i++)
			{
				for(int newBase = 0; newBase < 4; newBase++)
				{
					long newHash = e;
					long sub = e & (3L << (2*i));
					long add = newBase << (2*i);
					if(sub == add)
					{
						continue;
					}
					errored.add(new IndexedKmer(newHash+add-sub, count, index/2));
				}
			}
		}
		Collections.sort(errored);
		int n = errored.size(), m = database.length, i = 0, j = 0;
		long[] count = new long[error.length/2];
		int common = 0, totCount = 0;
		while(i < n && j < m)
		{
			long a = errored.get(i).kmer, b = database[j];
			if(a < b) i++;
			else if(a > b) j+=2;
			else
			{
				count[errored.get(i).index] += Math.min(errored.get(i).count, database[j+1]);
				i++;
			}
		}
		
		for(i = 1; i<error.length; i+=2)
		{
			totCount += error[i];
		}
		for(i = 0; i<count.length; i++)
		{
			common += count[i];
		}
		return 1.0 * common / totCount;
	}
	
	static class IndexedKmer implements Comparable<IndexedKmer>
	{
		long kmer;
		int index;
		int count;
		IndexedKmer(long kk, int cc, int ii)
		{
			kmer = kk; index = ii; count = cc;
		}
		@Override
		public int compareTo(IndexedKmer o) {
			// TODO Auto-generated method stub
			if(kmer != o.kmer) return Long.compare(kmer, o.kmer);
			return index - o.index;
		}
	}
	
	void print(long val)
	{
		char[] cs = new char[] {'A', 'C', 'G', 'T'};
		char[] res = new char[K1], res2 = new char[K1];
		for(int i = 0; i<K1; i++)
		{
			res[res.length-1-i] = cs[(int)(val%4)];
			res2[i] = cs[(int)(val%4) ^ 3];
			val  /= 4;
		}
		System.out.println(new String(res));
		System.out.println(new String(res2));
	}
	
	double[] containmentScore(Read r)
	{
		double[] res = new double[3];
		int common = 0, n = ms.length, m = r.ms.length;
		int i = 0, j = 0;
		while(i < n && j < m)
		{
			long a = ms[i], b = r.ms[j];
			if(a < b) i+=2;
			else if(a > b) j+=2;
			else
			{
				common+=Math.min(ms[i+1], r.ms[j+1]);
				i+=2;
				j+=2;
			}
		}
		
		int totCount = 0;
		for(j = 0; j<m; j += 2)
		{
			totCount += r.ms[j+1];
		}
		
		res[0] = 1.0 * common / totCount;
		
		int idx = 1;
		
		for(long[] x : new long[][]{r.starts, r.ends})
		{
		
			common = 0;
			n = ms.length;
			m = x.length;
			i = j = 0;
			while(i < n && j < m)
			{
				long a = ms[i], b = x[j];
				if(a < b) i+=2;
				else if(a > b) j+=2;
				else
				{
					common+=Math.min(ms[i+1], x[j+1]);
					i+=2;
					j+=2;
				}
			}
			
			totCount = 0;
			for(j = 0; j<m; j += 2)
			{
				totCount += x[j+1];
			}
			
			res[idx] = 1.0 * common / totCount;
			idx++;
		}
		
		return res;
	}
	
	// When sorting, sort by descending read length
	public int compareTo(Read o) {
		return o.len - len;
	}
}
}
