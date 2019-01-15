/*
 * FilterLengths filters reads, keeping longer reads, until the desired coverage is reached
 */
import java.util.*;
import java.io.*;
public class FilterLengths {
public static void main(String[] args) throws IOException
{
	String fn, ofn;
	int genomeLength, coverage;
	if(args.length == 4)
	{
		fn = args[0];
		genomeLength = Integer.parseInt(args[1]);
		coverage = Integer.parseInt(args[2]);
		ofn = args[3];
	}
	else
	{
		System.out.println("Usage: java FilterLengths <filename> <genomelength> <coverage> <outfilename>");
		return;
	}
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	ArrayList<Integer> lengths = new ArrayList<Integer>();
	while(input.hasNext())
	{
		String s = input.nextLine();
		if(s.startsWith("@"))
		{
			// FASTQ format
			lengths.add(input.nextLine().length());
			input.nextLine();
			input.nextLine();
		}
		else
		{
			lengths.add(input.nextLine().length());
		}
	}
	Collections.sort(lengths);
	Collections.reverse(lengths);
	long tot = (long)genomeLength * coverage;
	long soFar =  0;
	long minLength = lengths.get(0);
	for(long l : lengths)
	{
		soFar += l;
		if(soFar > tot) break;
		minLength = l;
	}
	input = new Scanner(new FileInputStream(new File(fn)));
	PrintWriter out = new PrintWriter(new File(ofn));
	while(input.hasNext())
	{
		String s = input.nextLine();
		String[] cur = new String[s.startsWith("@") ? 4 : 2];
		cur[0] = s;
		for(int i = 1; i<cur.length; i++) cur[i] = input.nextLine();
		int length = cur[1].length();
		if(length >= minLength)
		{
			for(String line : cur) out.println(line);
		}
	}
	out.close();
}
}
