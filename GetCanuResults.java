import java.util.*;
import java.io.*;
public class GetCanuResults {
public static void main(String[] args) throws IOException
{
	String dirname = "/home/mkirsche/ccs/canu_chr22.fastq.uncontained_hash.8_20_0.75_3.fastq/";
	if(args.length > 0)
	{
		dirname = args[0];
	}
	if(!dirname.endsWith("/")) dirname += "/";
	String fn = dirname + "chr22.report";
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	int ng50 = -1;
	int nContigs = -1;
	int length = -1;
	int nReads = -1;
	while(input.hasNext())
	{
		String line = input.nextLine().trim();
		if(line.contains("NG (bp)"))
		{
			while(input.hasNext())
			{
				String s = input.nextLine().trim();
				String[] tokens = s.split("\\s+");
				if(tokens.length > 2 && tokens[0].equals("--") && tokens[1].equals("50"))
				{
					ng50 = Integer.parseInt(tokens[2]);
					break;
				}
			}
		}
		else if(line.contains("contigs:"))
		{
			String[] tokens = line.split("\\s+");
			if(tokens.length > 2) nContigs = Integer.parseInt(tokens[2]);
			if(tokens.length > 6) length = Integer.parseInt(tokens[6]);
		}
		else if(line.contains("reads."))
		{
			String[] tokens = line.split("\\s+");
			if(tokens.length == 4 && tokens[1].equals("Found"))
			{
				nReads = Integer.parseInt(tokens[2]);
			}
		}
	}
	boolean printParams = printParams(dirname);
	if(!printParams) System.out.println("Trial: " + dirname);
	System.out.println("Reads: " + nReads);
	System.out.println("ng50: " + ng50);
	System.out.println("Contigs: " + nContigs);
	System.out.println("Length: " + length);
}
static boolean printParams(String dirname)
{
	if(!dirname.contains("hash")) return false;
	try {
		String[] tokens = dirname.split("_");
		int kmerFreq = (1<<Integer.parseInt(tokens[2].substring(1 + tokens[2].indexOf('.'))));
		System.out.println("Kmer_sampling_frequency: 1/" + kmerFreq);
		System.out.println("Kmer_length: " + Integer.parseInt(tokens[3]));
		String threshold = String.format("%.2f", Double.parseDouble(tokens[4]));
		System.out.println("Containment_indentity_threshold: " + threshold);
		System.out.println("Seed_kmers: " + Integer.parseInt(tokens[5].substring(0, tokens[5].indexOf('.'))));
		return true;
	} catch (Exception e) {
		return false;
	}
}
}
