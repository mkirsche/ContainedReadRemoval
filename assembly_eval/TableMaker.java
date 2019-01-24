import java.util.*;
import java.io.*;
public class TableMaker {
public static void main(String[] args) throws Exception
{
	if(args.length != 3)
	{
		System.out.println("Usage: java TableMaker <assemblyDir> <statsDir> <readsetsDir>");
		return;
	}
	String assemblyDir = args[0];
	String statsDir = args[1];
	String readsetsDir = args[2];
	String command = "ls " + readsetsDir;
	Process child = Runtime.getRuntime().exec(command);
    InputStream seqStream = child.getInputStream();
    Scanner seqInput = new Scanner(seqStream);
    HashMap<String, Long> totReadLengthMap = new HashMap<String, Long>();
    HashMap<String, Integer> numReadsMap = new HashMap<String, Integer>();
    while(seqInput.hasNext())
    {
    	// Get lines and total bases
    	String basename = seqInput.next();
    	String fn = readsetsDir + "/" + basename;
    	boolean fastq = false;
    	Scanner sc = new Scanner(new FileInputStream(new File(fn)));
    	long totLength = 0;
    	int numReads = 0;
    	if(sc.hasNext())
    	{
    		String s = sc.nextLine();
        	if(s.startsWith("@")) fastq = true;
        	numReads++;
        	totLength += sc.nextLine().length();
        	if(fastq) for(int i = 0; i<2; i++) sc.nextLine();
        	while(sc.hasNext())
        	{
        		sc.nextLine();
        		numReads++;
        		totLength += sc.nextLine().length();
            	if(fastq) for(int i = 0; i<2; i++) sc.nextLine();
        	}
    	}
    	numReadsMap.put(basename, numReads);
    	totReadLengthMap.put(basename, totLength);
    }
    
	command = "ls " + assemblyDir;
    child = Runtime.getRuntime().exec(command);
    seqStream = child.getInputStream();
    seqInput = new Scanner(seqStream);
    ArrayList<Result> all = new ArrayList<Result>();
    while(seqInput.hasNext())
    {
    	String assemblyFn = seqInput.next();
    	String curDir = statsDir + "/" + assemblyFn + "_stats";
    	String statsFn = curDir + "/assemblystats.txt";
    	Scanner input = new Scanner(new FileInputStream(new File(statsFn)));
    	Result res = new Result();
    	res.name = assemblyFn;
    	String[] tokens = input.nextLine().split(" ");
    	res.numContigs = Integer.parseInt(tokens[tokens.length-1]);
    	if(res.numContigs == 0)
    	{
    		res.length = 0;
    		res.n50 = 0;
    		res.quast = 0;
    		res.busco = 0;
    		all.add(res);
    		continue;
    	}
    	tokens = input.nextLine().split(" ");
    	res.length = Integer.parseInt(tokens[tokens.length-1]);
    	tokens = input.nextLine().split(" ");
    	res.n50 = Integer.parseInt(tokens[tokens.length-1]);
    	String quastFn = curDir + "/quast/report.txt";
    	res.quast = getQuastScore(quastFn);
    	String buscoFn = curDir + "/run_busco/short_summary_busco.txt";
    	res.busco = getBuscoScore(buscoFn);
    	all.add(res);
    }
    
    updateNames(all, totReadLengthMap, numReadsMap);
    System.out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
    		"Name", "n50", "Length", "# Contigs", "Quast", "Busco", "NumReads", "TotalReadLength");
    for(Result r : all)
    {
    	System.out.println(r);
    }
}
static void updateNames(ArrayList<Result> rs, HashMap<String, Long> trlm, HashMap<String, Integer> nrm)
{
	for(Result r : rs)
	{
		for(String s : trlm.keySet())
		{
			if(r.name.contains(s))
			{
				r.totLength = trlm.get(s);
			}
		}
		for(String s : nrm.keySet())
		{
			if(r.name.contains(s))
			{
				r.numReads = nrm.get(s);
			}
		}
	}
}
static double getBuscoScore(String fn) throws IOException
{
	try
	{
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		while(input.hasNext())
		{
			String line = input.nextLine().trim();
			if(line.startsWith("C:") && line.indexOf('%') != -1)
			{
				return 0.01 * Double.parseDouble(line.substring(2, line.indexOf('%')));
			}
		}
	} 
	catch(Exception e)
	{
		return 0;
	}
	return 0;
}
static double getQuastScore(String fn) throws IOException
{
	try
	{
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.startsWith("Genome fraction"))
			{
				String[] tokens = line.split(" ");
				return 0.01 * Double.parseDouble(tokens[tokens.length-1]);
			}
		}
	}
	catch(Exception e)
	{
		return 0;
	}
	return 0;
}

static class Result
{
	String name;
	int n50, length, numContigs;
	double quast, busco;
	int numReads;
	long totLength;
	Result()
	{
		
	}
	Result(String name, int n50, int length, int numContigs, double quast, double busco)
	{
		this.name = name;
		this.n50 = n50;
		this.length = length;
		this.numContigs = numContigs;
		this.busco = busco;
		this.quast = quast;
	}
	public String toString()
	{
		return name + "\t" + n50 + "\t" + length + "\t" + numContigs + "\t" + quast + "\t" + busco
				+ "\t" + numReads + "\t" + totLength;
	}
}
}
