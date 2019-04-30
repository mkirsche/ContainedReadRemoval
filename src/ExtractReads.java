import java.util.*;
import java.io.*;
public class ExtractReads {
public static void main(String[] args) throws IOException
{
	String fn = "ERR2173373.fastq.8.14.5.500.10";
	String readsFn = "ERR2173373.fastq";
	if(args.length > 0)
	{
		if(args.length == 1)
		{
			System.out.println("readlistfilename readfilename");
		}
		else
		{
			fn = args[0];
			readsFn = args[1];
		}
	}
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	boolean fastq = readsFn.contains("fastq") || readsFn.contains("fq");
	HashSet<String> rs = new HashSet<String>();
	while(input.hasNext())
	{
		rs.add(input.nextLine());
	}
	int linesPer = fastq ? 4 : 2;
	String[] buf = new String[linesPer];
	BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(readsFn)));
	input = new Scanner(new FileInputStream(new File(readsFn)));
	int count = 0;
	PrintWriter out = new PrintWriter(new File(fn + ".fastq"));
	int ii = 0;
	long time = System.currentTimeMillis();
	while(true)
	{
		try
		{
			for(int i = 0; i<linesPer; i++) buf[i] = in.readLine();
			if(rs.contains(buf[0].substring(1)))
			{
				count++;
				for(int i = 0; i<linesPer; i++) out.println(buf[i]);
			}
		} 
		catch(Exception e)
		{
			break;
		}
	}
	out.close();
	System.out.println(count);
}
}