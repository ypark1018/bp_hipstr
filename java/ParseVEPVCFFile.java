import java.io.*;
import java.util.*;
import java.util.zip.*;

public class ParseVEPVCFFile {

    public void printErrorMessage(String error) {
        System.err.println("ERROR_in_ParseVEPVCFFile : "+error);
	System.exit(-1);
    }

    public int findIndex(String line, String field) throws Exception {
	int startIndex = line.indexOf("Format:");
	if (startIndex < 0) {
	    printErrorMessage("cannot find Format: in "+line);
	}
	String allfield = line.substring(startIndex+8, line.length()-2);
	String[] allfieldtoken = allfield.split("\\|");
	for (int i = 0; i < allfieldtoken.length; i++) {
	    if (allfieldtoken[i].compareTo(field) == 0) {
		return i;
	    }
	}
	printErrorMessage("cannot find "+field+" in "+allfield);
	return -1;
    }

    public void splitString(String str, char delim, ArrayList<String> token) throws Exception {
	if (str.length() == 0) {
	    return;
	}
	if (str.charAt(0) == delim) {
	    printErrorMessage(str+" starts with "+delim);
	}
	boolean done = false;
	int curpos = 0;
	while (!done) {
	    StringBuffer sb = new StringBuffer();
	    while (curpos < str.length() && str.charAt(curpos) != delim) {
		sb.append(str.charAt(curpos));
		curpos++;
	    }
	    if (sb.length() == 0) {
		token.add("");
	    }
	    else {
		token.add(sb.toString());
	    }
	    curpos++;
	    if (curpos == str.length()) {
		done = true;
	    }
	}
    }

    public class Annot {
	String allele;
	String symbol;
	String conseq;
	String impact;
	int count;
    }

    public class AnnotComparator implements Comparator<Object>  {
        public int compare(Object o1, Object o2) {
            Annot a1 = (Annot)o1;
            Annot a2 = (Annot)o2;
	    int impact1 = -1;
	    int impact2 = -1;
	    if (a1.impact.compareTo("HIGH") == 0) {
		impact1 = 1;
	    }
	    else if (a1.impact.compareTo("MODERATE") == 0) {
		impact1 = 2;
	    }
	    else if (a1.impact.compareTo("LOW") == 0) {
		impact1 = 3;
	    }
	    else if (a1.impact.compareTo("MODIFIER") == 0) {
		impact1 = 4;
	    }
	    else {
		printErrorMessage("undefined impact "+a1);
	    }
	    if (a2.impact.compareTo("HIGH") == 0) {
		impact2 = 1;
	    }
	    else if (a2.impact.compareTo("MODERATE") == 0) {
		impact2 = 2;
	    }
	    else if (a2.impact.compareTo("LOW") == 0) {
		impact2 = 3;
	    }
	    else if (a2.impact.compareTo("MODIFIER") == 0) {
		impact2 = 4;
	    }
	    else {
		printErrorMessage("undefined impact "+a2);
	    }
	    if (impact1 == impact2) {
		if (a1.symbol.compareTo("NotMapped") == 0 && a2.symbol.compareTo("NotMapped") != 0) {
		    return 1;
		}
		if (a1.symbol.compareTo("NotMapped") != 0 && a2.symbol.compareTo("NotMapped") == 0) {
		    return -1;
		}
		if (a1.count < a2.count) {
		    return 1;
		}
		else {
		    return -1;
		}
	    }
	    else {
		if (impact1 < impact2) {
		    return -1;
		}
		else {
		    return 1;
		}
	    }
	    //return 0;
        }
    }

    // TO DO:
    // Need to output all alleles with HIGH, MODERATE, LOW impact for each variant (for multi-allelic variants)
    public void createFiles(String inputvcffile, String outputfile) throws Exception {
	BufferedReader in = null;
	if (inputvcffile.endsWith(".gz")) {
	    in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputvcffile))));
	}
	else {
	    in = new BufferedReader(new FileReader(inputvcffile));
	}
	String line = in.readLine();
	while (!line.startsWith("##INFO=<ID=CSQ")) {
	    line = in.readLine();
	}
	int alleleIndex = findIndex(line, "Allele");
	int conseqIndex = findIndex(line, "Consequence");
	int impactIndex = findIndex(line, "IMPACT");
	int symbolIndex = findIndex(line, "SYMBOL");

	System.out.println("index of Allele = "+alleleIndex);
	System.out.println("index of Consequence = "+conseqIndex);
	System.out.println("index of IMPACT = "+impactIndex);
	System.out.println("index of SYMBOL = "+symbolIndex);

	while (!line.startsWith("#CHROM")) {
	    line = in.readLine();
	}
	line = in.readLine();
	PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outputfile)));
	while (line != null) {
	    String[] token = line.split("\\s");
	    String chr = token[0];
	    String pos = token[1];
	    String info = token[7];
	    String[] infotoken = info.split(";");
	    String vepannot = infotoken[infotoken.length-1];
	    if (!vepannot.startsWith("CSQ")) {
		printErrorMessage("VEP annotation does not start with CSQ in "+line);
	    }
	    String[] temptoken = vepannot.split("=");
	    String[] allveptoken = temptoken[1].split(",");
	    Hashtable<String,Integer> annotInfo = new Hashtable<String,Integer>();
	    for (int i = 0; i < allveptoken.length; i++) {
		//System.out.println(allveptoken[i]);
		ArrayList<String> veptoken = new ArrayList<String>();
		splitString(allveptoken[i], '|', veptoken);
		//String[] veptoken = allveptoken[i].split("\\|");
		String allele = veptoken.get(alleleIndex);
		String conseq = veptoken.get(conseqIndex);
		String impact = veptoken.get(impactIndex);
		String symbol = veptoken.get(symbolIndex);
		if (allele.length() == 0) {
		    printErrorMessage("allele has zero length "+chr+" "+pos);
		}
		if (conseq.length() == 0) {
		    printErrorMessage("consequence has zero length "+chr+" "+pos);
		}
		if (conseq.compareTo("?") == 0) {
		    System.out.println("unknown consequence "+chr+" "+pos);
		    continue;
		}
		if (impact.length() == 0) {
		    printErrorMessage("impact has zero length "+chr+" "+pos);
		}
		if (symbol.length() == 0) {
		    symbol = "NotMapped";
		}
		String merged = allele+";"+symbol+";"+conseq+";"+impact;
		if (annotInfo.containsKey(merged)) {
		    int count = annotInfo.get(merged);
		    count++;
		    annotInfo.put(merged, count);
		}
		else {
		    annotInfo.put(merged, 1);
		}
		//System.out.println(chr+" "+pos+" "+conseq+" "+impact+" "+symbol);
	    }
	    Enumeration<String> key = annotInfo.keys();
	    TreeSet<Annot> annotSet = new TreeSet<Annot>(new AnnotComparator());
	    while (key.hasMoreElements()) {
		String merged = key.nextElement();
		int count = annotInfo.get(merged);
		String[] mergedtoken = merged.split(";");
		Annot ann = new Annot();
		ann.allele = mergedtoken[0];
		ann.symbol = mergedtoken[1];
		ann.conseq = mergedtoken[2];
		ann.impact = mergedtoken[3];
		ann.count = count;
		annotSet.add(ann);
	    }

	    Iterator<Annot> ite = annotSet.iterator();
	    Annot ann = ite.next();
	    out.println(chr+" "+pos+" "+ann.allele+" "+ann.impact+" "+ann.symbol+" "+ann.conseq);
	    /*
	    Iterator<Annot> ite = annotSet.iterator();
	    while (ite.hasNext()) {
		Annot ann = ite.next();
		System.out.println(chr+" "+pos+" "+ann.allele+" "+ann.impact+" "+ann.symbol+" "+ann.conseq+" "+ann.count);
	    }
	    System.out.println(); 
	    */
	    /*
	    Annot ann = ite.next();
	    if (ann.symbol.compareTo("NotMapped") == 0 && annotSet.size() > 1) {
		System.out.println(chr+" "+pos+" "+ann.impact+" "+ann.symbol+" "+ann.conseq+" "+ann.count);
		while (ite.hasNext()) {
		    ann = ite.next();
		    System.out.println(chr+" "+pos+" "+ann.impact+" "+ann.symbol+" "+ann.conseq+" "+ann.count);
		}
		System.out.println();
	    }
	    */
	    line = in.readLine();
	}
	in.close();
	out.close();
    }

    public static void main(String[] args) throws Exception {
	if (args.length != 2) {
	    System.err.println("[usage] java ParseVEPVCFFile [input VEP VCF file] [output dir]");
	    System.exit(-1);
	}
	ParseVEPVCFFile cc = new ParseVEPVCFFile();
	cc.createFiles(args[0], args[1]);
    }
}

