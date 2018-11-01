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
	}
	System.out.println("Trial: " + dirname);
	System.out.println("ng50: " + ng50);
	System.out.println("Contigs: " + nContigs);
	System.out.println("Length: " + length);
}
}
