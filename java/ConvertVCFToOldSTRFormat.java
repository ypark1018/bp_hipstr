import java.io.*;
import java.util.zip.*;
import java.util.*;

public class ConvertVCFToOldSTRFormat {

    public class Allele {
        String allele;
    }

    public class AlleleComparator implements Comparator<Object> {
        public int compare(Object o1, Object o2) {
            Allele ind1 = (Allele)o1;
            Allele ind2 = (Allele)o2;
            if (ind1.allele.length() > ind2.allele.length()) return 1;
            return -1;
        }
    }

    // TO DO:
    // Write a file that for each STR, write # of alleles, length for each allele sorted by length, and then which REF or ALT corresponds to (ATL1 REF ALT2 ALT3)
    // Write a file that encodes index of allele (1-N) not length
    // Write a file for long format

    public void createFiles(String vcffile, String outputfile) throws Exception {
	int gtIndex = 0;
        BufferedReader in = new BufferedReader(new FileReader(vcffile));
	PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outputfile+".allele.length.wide.txt")));
	PrintWriter out1 = new PrintWriter(new BufferedWriter(new FileWriter(outputfile+".str.length.allele.info.txt")));
	PrintWriter out2 = new PrintWriter(new BufferedWriter(new FileWriter(outputfile+".allele.index.wide.txt")));
	PrintWriter out3 = new PrintWriter(new BufferedWriter(new FileWriter(outputfile+".allele.length.long.txt")));
	PrintWriter out4 = new PrintWriter(new BufferedWriter(new FileWriter(outputfile+".allele.index.long.txt")));

	String line = in.readLine();
	while (!line.startsWith("#CHROM")) {
	    line = in.readLine();
	}
	String[] htoken = line.split("\\s");
	int numIndiv = htoken.length-9;
	ArrayList<String> idList = new ArrayList<String>();
	out.print("#CHR:POS");
	out1.println("CHR:POS\t#_alleles\tLength_of_all_alleles\tREF_ALT_allele");
	out2.print("#CHR:POS");
	for (int i = 9; i < htoken.length; i++) {
	    idList.add(htoken[i]);
	    out.print("\t"+htoken[i]+"_first_chr\t"+htoken[i]+"_second_chr");
	    out2.print("\t"+htoken[i]+"_first_chr\t"+htoken[i]+"_second_chr");
	}
	out.println();
	out2.println();
	out3.println("#ID\tSTR\tGenotype");
	out4.println("#ID\tSTR\tGenotype");
	line = in.readLine();
	while (line != null) {
	    TreeSet<Allele> alleleSet = new TreeSet<Allele>(new AlleleComparator());
	    Hashtable<Integer,Integer> alleleLengthInfo = new Hashtable<Integer,Integer>();
	    Hashtable<Integer,String> lengthAlleleInfo = new Hashtable<Integer,String>();
	    String[] token = line.split("\\s");
	    String chr = token[0];
	    String pos = token[1];
	    String snpid = chr+":"+pos;
	    String rallele = token[3];
	    Allele refallele = new Allele();
	    lengthAlleleInfo.put(rallele.length(), "REF");

	    out.print(snpid);
	    out2.print(snpid);

	    refallele.allele = rallele;
	    alleleSet.add(refallele);
	    // 0 is reference allele
	    alleleLengthInfo.put(0, rallele.length());
	    String[] alttoken = token[4].split(",");
	    for (int i = 0; i < alttoken.length; i++) {
		String aallele = alttoken[i];
		if (aallele.compareTo(rallele) == 0) {
		    System.err.println(snpid+" has the same REF and ALT alleles "+refallele);
		    System.exit(-1);
		}
		Allele altallele = new Allele();
		altallele.allele = aallele;
		alleleSet.add(altallele);
		alleleLengthInfo.put(i+1, aallele.length());
		lengthAlleleInfo.put(aallele.length(), "ALT"+(i+1));
	    }

	    out1.print(snpid+"\t"+alleleSet.size());

	    Iterator<Allele> ite = alleleSet.iterator();
	    Hashtable<Integer,Integer> lengthIndexInfo = new Hashtable<Integer,Integer>();
	    int index = 0;
	    while (ite.hasNext()) {
		Allele tempallele = ite.next();
		int length = tempallele.allele.length();
		lengthIndexInfo.put(length, index);
		index++;
		out1.print("\t"+length);
		//System.out.println(length);
	    }

	    ite = alleleSet.iterator();
	    while (ite.hasNext()) {
		Allele tempallele = ite.next();
		int length = tempallele.allele.length();
		if (!lengthAlleleInfo.containsKey(length)) {
		    System.err.println("cannot find ref or alt allele with length = "+length);
		    System.exit(-1);
		}
		out1.print("\t"+lengthAlleleInfo.get(length));
		//System.out.println(length);
	    }
	    out1.println();

	    for (int i = 9; i < token.length; i++) {
		out3.print(idList.get(i-9)+"\t"+snpid);
		out4.print(idList.get(i-9)+"\t"+snpid);
		// if (token[i].compareTo("./.") == 0)  {
		if (token[i].compareTo(".") == 0)  {
		    out.print("\tNA\tNA");
		    out2.print("\tNA\tNA");
		    out3.println("\tNA\tNA");
		    out4.println("\tNA\tNA");
		    // missing
		}
		else {
		    String[] itoken = token[i].split(":");
		    String genotype = itoken[gtIndex];
		    String[] gtoken = genotype.split("\\|"); // "/" "\\|"
		    //System.out.println(gtoken[0]); ////////
		    int fallele = Integer.parseInt(gtoken[0]);
		    int sallele = Integer.parseInt(gtoken[1]);
		    if (!alleleLengthInfo.containsKey(fallele)) {
			System.err.println("cannot find allele for genotype "+genotype+" at "+snpid);
			System.exit(-1);
		    }
		    if (!alleleLengthInfo.containsKey(sallele)) {
			System.err.println("cannot find allele for genotype "+genotype+" at "+snpid);
			System.exit(-1);
		    }
		    int flength = alleleLengthInfo.get(fallele);
		    int slength = alleleLengthInfo.get(sallele);
		    if (!lengthIndexInfo.containsKey(flength)) {
			System.err.println("cannot find index for genotype "+genotype+" at "+snpid);
			System.exit(-1);
		    }
		    if (!lengthIndexInfo.containsKey(slength)) {
			System.err.println("cannot find index for genotype "+genotype+" at "+snpid);
			System.exit(-1);
		    }
		    int findex = lengthIndexInfo.get(flength)+1;
		    int sindex = lengthIndexInfo.get(slength)+1;
		    out.print("\t"+flength+"\t"+slength);
		    out2.print("\t"+findex+"\t"+sindex);
		    out3.println("\t"+flength+"\t"+slength);
		    out4.println("\t"+findex+"\t"+sindex);
		}
	    }
	    line = in.readLine();
	    out.println();
	    out2.println();
	    //System.exit(-1);
	}
	in.close();
	out.close();
	out1.close();
	out2.close();
	out3.close();
	out4.close();
    }

    public static void main(String[] args) throws Exception {
	if (args.length != 2) {
	    System.err.println("[usage] java ConvertVCFToOldSTRFormat [vcf file] [output file]");
	    System.exit(-1);
	}
	ConvertVCFToOldSTRFormat cc = new ConvertVCFToOldSTRFormat();
	cc.createFiles(args[0], args[1]);
    }
}

