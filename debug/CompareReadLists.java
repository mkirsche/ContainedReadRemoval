import java.util.*;
import java.io.*;
public class CompareReadLists {
public static void main(String[] args) throws IOException
{
	String fn1 = "/home/mkirsche/ccs/chr22_refuncontained.txt";
	String fn2 = "/home/mkirsche/ccs/chr22.fastq.uncontained_hash.8_20_0.9";
	if(args.length > 0)
	{
		if(args.length == 1)
		{
			System.out.println("fn1 fn2");
			return;
		}
		fn1 = args[0];
		fn2 = args[1];
	}
	Scanner input1 = new Scanner(new FileInputStream(new File(fn1)));
	Scanner input2 = new Scanner(new FileInputStream(new File(fn2)));
	HashSet<String> r1 = new HashSet<String>();
	int n1 = 0, n2 = 0, inter = 0;
	while(input1.hasNext()) r1.add(input1.nextLine());
	n1 = r1.size();
	while(input2.hasNext())
	{
		String s = input2.nextLine();
		n2++;
		if(r1.contains(s)) inter++;
	}
	System.out.println(fn1 + ": " + n1);
	System.out.println(fn2 + ": " + n2);
	System.out.println("Shared: " + inter);
}
}
