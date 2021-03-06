/*
 * Outputs the total assembly length, number of contigs, and contig N50 of a wtdbg2 assembly
 * Assumes the contig name lines have a second token containing the contig length 
 */
import java.util.*;
import java.io.*;
public class AssemblyStats {
public static void main(String[] args) throws IOException
{
	String fn = args[0];
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	ArrayList<Long> lengths = new ArrayList<Long>();
	long curLen = 0;
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(!line.startsWith(">")) curLen += line.length();
		else
		{
		    if(curLen != 0) lengths.add(curLen);
		    curLen = 0;
		}
	}
	if(curLen > 0) lengths.add(curLen);
	Collections.sort(lengths);
	long total = 0;
	System.out.println("Number of contigs: " + lengths.size());
	for(long x : lengths) total += x;
	System.out.println("Total length: " + total);
	long half = total/2;
	for(int i = lengths.size()-1; i>=0; i--)
	{
		half -= lengths.get(i);
		if(half <= 0)
		{
			System.out.println("N50: " + lengths.get(i));
			break;
		}
	}
	System.out.println();
}
}
