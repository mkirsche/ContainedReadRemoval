import java.util.*;
import java.io.*;
public class ExtractReads {
public static void main(String[] args) throws IOException
{
	String fn = "/home/mkirsche/ccs/chr22.fastq.uncontained_hash.8_10_0.95";
	String readsFn = "/home/mkirsche/ccs/chr22.fastq";
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
	HashSet<String> rs = new HashSet<String>();
	while(input.hasNext()) rs.add(input.nextLine());
	String[] buf = new String[4];
	input = new Scanner(new FileInputStream(new File(readsFn)));
	int count = 0;
	PrintWriter out = new PrintWriter(new File(fn + ".fastq"));
	while(input.hasNext())
	{
		for(int i = 0; i<4; i++) buf[i] = input.nextLine();
		if(rs.contains(buf[0].substring(1)))
		{
			count++;
			for(int i = 0; i<4; i++) out.println(buf[i]);
		}
	}
	out.close();
	System.out.println(count);
}
}
