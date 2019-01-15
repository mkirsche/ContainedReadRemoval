import java.util.*;
import java.io.*;
public class TableMaker {
public static void main(String[] args) throws Exception
{
	if(args.length != 2)
	{
		System.out.println("Usage: java TableMaker <assemblyDir> <statsDir>");
		return;
	}
	String assemblyDir = args[0];
	String statsDir = args[1];
	String command = "ls " + assemblyDir;
    Process child = Runtime.getRuntime().exec(command);
    InputStream seqStream = child.getInputStream();
    Scanner seqInput = new Scanner(seqStream);
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
    System.out.println();
    for(Result r : all)
    {
    	System.out.println(r);
    }
}
static double getBuscoScore(String fn) throws IOException
{
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	while(input.hasNext())
	{
		String line = input.nextLine().trim();
		if(line.startsWith("C:") && line.indexOf('[') != -1)
		{
			return 0.01 * Double.parseDouble(line.substring(2, line.indexOf('[')));
		}
	}
	return 0;
}
static double getQuastScore(String fn) throws IOException
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
	return 0;
}
static class Result
{
	String name;
	int n50, length, numContigs;
	double quast, busco;
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
		return name + "\t" + n50 + "\t" + length + "\t" + numContigs + "\t" + quast + "\t" + busco;
	}
}
}
