import java.util.*;
import java.io.*;
public class GetBestAssemblies {
public static void main(String[] args) throws IOException
{
	if(args.length != 2)
	{
		System.out.println("Usage java GetBestAssemblies <tablefile> <refname>");
		return;
	}
	String fn = args[0];
	String refName = args[1];
	ArrayList<String> names = getGoodAssemblies(fn, refName);
	for(String s : names) System.out.println(s);
}
static double prop = .5;
static ArrayList<String> getGoodAssemblies(String fn, String refName) throws IOException
{
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	String[] fields = input.nextLine().split("\t"); // parse header
	int nameIdx = -1, n50Idx = -1;
	for(int i = 0; i<fields.length; i++)
	{
		if(fields[i].equalsIgnoreCase("name")) nameIdx = i;
		else if(fields[i].equalsIgnoreCase("n50")) n50Idx = i;
	}
	ArrayList<String> res = new ArrayList<String>();
	if(nameIdx == -1 || n50Idx == -1) return res;
	ArrayList<String> names = new ArrayList<String>();
	ArrayList<Integer> n50s = new ArrayList<Integer>();
	while(input.hasNext())
	{
		String[] cur = input.nextLine().split("\t");
		names.add(cur[nameIdx]);
		n50s.add(Integer.parseInt(cur[n50Idx]));
	}
	int baseline = 0;
	for(int i = 0; i<names.size(); i++)
	{
		if(names.get(i).equalsIgnoreCase(refName))
		{
			baseline = n50s.get(i);
		}
	}
	for(int i = 0; i<names.size(); i++)
	{
		if(n50s.get(i) >= baseline * prop)
		{
			res.add(names.get(i));
		}
	}
	return res;
}
}
