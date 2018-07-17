import java.io.*;
import java.util.zip.*;
import java.util.*;

public class RecodeBothOldNewSTR {

    public void findNonNAPos(String[] token1, String[] token2, ArrayList<Integer> nonNAPosList, String strID) throws Exception {
	if (token1.length != token2.length) {
	    System.err.println("different line length");
	    System.out.println(token1);
	    System.out.println(token2);
	    System.exit(-1);
	}
	for (int i = 1; i < token1.length; i=i+2) {
	    if (token1[i].compareTo("NA")==0) {
		if (token1[i+1].compareTo("NA") != 0) {
		    System.err.println("both chromosomes must be NA "+strID);
		    System.exit(-1);
		}
	    }
	    else if (token2[i].compareTo("NA") == 0) {
		if (token2[i+1].compareTo("NA") != 0) {
		    System.err.println("both chromosomes must be NA "+strID);
		    System.exit(-1);
		}
	    }
	    else {
		nonNAPosList.add(i);
	    }
	}
    }

    public void findNewAlleleCoding(String[] token, ArrayList<Integer> nonNAPosList, Hashtable<Integer,Integer> alleleCodeInfo, TreeSet<Integer> alleleSet) throws Exception {
	Hashtable<Integer,Integer> alleleInfo = new Hashtable<Integer,Integer>();
	for (int i = 0; i < nonNAPosList.size(); i++) {
	    int index = nonNAPosList.get(i);
	    int fallele = Integer.parseInt(token[index]);
	    int sallele = Integer.parseInt(token[index+1]);
	    if (!alleleInfo.containsKey(fallele)) {
		alleleInfo.put(fallele,0);
		alleleSet.add(fallele);
	    }
	    if (!alleleInfo.containsKey(sallele)) {
		alleleInfo.put(sallele,0);
		alleleSet.add(sallele);
	    }
	}
	Iterator<Integer> ite = alleleSet.iterator();
	int codeIndex = 1;
	while (ite.hasNext()) {
	    int allele = ite.next();
	    alleleCodeInfo.put(allele, codeIndex);
	    codeIndex++;
	}
    }

    public void convertToNewCode(String[] token, ArrayList<Integer> nonNAPosList, Hashtable<Integer,Integer> alleleCodeInfo, ArrayList<String> outputList) throws Exception {
	String output = token[0];
	for (int i = 1; i < token.length; i=i+2) {
	    if (!nonNAPosList.contains(i)) {
		output = output+" NA NA";
	    }
	    else {
		int fallele = Integer.parseInt(token[i]);
		int sallele = Integer.parseInt(token[i+1]);
		if (!alleleCodeInfo.containsKey(fallele)) {
		    System.err.println("cannot find allele coding for "+fallele);
		    System.exit(-1);
		}
		if (!alleleCodeInfo.containsKey(sallele)) {
		    System.err.println("cannot find allele coding for "+sallele);
		    System.exit(-1);
		}
		output = output+" "+alleleCodeInfo.get(fallele)+" "+alleleCodeInfo.get(sallele);
	    }
	}
	outputList.add(output);
    }

    public void inferAlleleCode(TreeSet<Integer> alleleSet1, Hashtable<Integer,Integer> alleleInfo1, TreeSet<Integer> alleleSet2, Hashtable<Integer,Integer> alleleInfo2) throws Exception {
	TreeSet<Integer> largeSet = null;
	TreeSet<Integer> smallSet = null;
	Hashtable<Integer,Integer> largeInfo = null;
	Hashtable<Integer,Integer> smallInfo = null;
	if (alleleSet1.size() < alleleSet2.size()) {
	    largeSet = alleleSet2;
	    largeInfo = alleleInfo2;
	    smallSet = alleleSet1;
	    smallInfo = alleleInfo1;
	}
	else {
	    largeSet = alleleSet1;
	    largeInfo = alleleInfo1;
	    smallSet = alleleSet2;
	    smallInfo = alleleInfo2;
	}
	ArrayList<Integer> smallList = new ArrayList<Integer>();
	ArrayList<Integer> largeList = new ArrayList<Integer>();
	Iterator<Integer> ite1 = smallSet.iterator();
	while(ite1.hasNext()) {
	    smallList.add(ite1.next());
	}
	Iterator<Integer> ite2 = largeSet.iterator();
	while(ite2.hasNext()) {
	    largeList.add(ite2.next());
	}
	int numMaxMatch = -1;
	int maxIndexSmall = -1;
	int maxIndexLarge = -1;
	for (int n = 0; n < smallList.size(); n++) {
	    for (int i = 0; i < largeList.size(); i++) {
		int numMatch = 0;
		int firstPosSmall = smallList.get(n);
		int firstPosLarge = largeList.get(i);
		for (int j = n+1; j < smallList.size(); j++) {
		    int pos = smallList.get(j);
		    int diff = pos - firstPosSmall;
		    if (largeInfo.containsKey(firstPosLarge+diff)) {
			numMatch++;
		    }
		}
		if (numMatch > numMaxMatch) {
		    numMaxMatch = numMatch;
		    maxIndexLarge = i;
		    maxIndexSmall = n;
		}
	    }
	}
	//System.out.println("numMaxMatch = "+numMaxMatch+" maxIndexSmall = "+maxIndexSmall+" maxIndexLarge = "+maxIndexLarge);
	int newbadindexSmall = -1;
	for (int i = 0; i < maxIndexSmall; i++) {
	    int pos = smallList.get(i);
	    smallInfo.put(pos,newbadindexSmall);
	    newbadindexSmall--;
	}
	int newbadindexLarge = -1;
	for (int i = 0; i < maxIndexLarge; i++) {
	    int pos = largeList.get(i);
	    largeInfo.put(pos,newbadindexLarge);
	    newbadindexLarge--;
	}

	int newgoodindexSmall = 1;
	int maxIndexSmallBP = smallList.get(maxIndexSmall);
	smallInfo.put(maxIndexSmallBP,newgoodindexSmall);
	newgoodindexSmall++;

	int newgoodindexLarge = 1;
	int maxIndexLargeBP = largeList.get(maxIndexLarge);
	largeInfo.put(maxIndexLargeBP,newgoodindexLarge);
	newgoodindexLarge++;

	for (int i = maxIndexSmall+1; i < smallList.size(); i++) {
	    int pos = smallList.get(i);
	    int diff = pos - maxIndexSmallBP;
	    if (largeInfo.containsKey(maxIndexLargeBP+diff)) {
		smallInfo.put(pos, newgoodindexSmall);
		newgoodindexSmall++;
	    }
	    else {
		smallInfo.put(pos, newbadindexSmall);
		newbadindexSmall--;
	    }
	}

	for (int i = maxIndexLarge+1; i < largeList.size(); i++) {
	    int pos = largeList.get(i);
	    int diff = pos - maxIndexLargeBP;
	    if (smallInfo.containsKey(maxIndexSmallBP+diff)) {
		largeInfo.put(pos, newgoodindexLarge);
		newgoodindexLarge++;
	    }
	    else {
		largeInfo.put(pos, newbadindexLarge);
		newbadindexLarge--;
	    }
	}
	/*
	System.out.println("Small List");
	for (int i = 0; i < smallList.size(); i++) {
	    System.out.println("Allele = "+smallList.get(i)+" New allele = "+smallInfo.get(smallList.get(i)));
	}
	System.out.println("Large List");
	for (int i = 0; i < largeList.size(); i++) {
	    System.out.println("Allele = "+largeList.get(i)+" New allele = "+largeInfo.get(largeList.get(i)));
	}
	*/
    }

    public void recodeTwoLines(String line1, String line2, ArrayList<String> outputList, PrintWriter oldout, PrintWriter newout) throws Exception {
	String[] token1 = line1.split("\\s");
	String[] token2 = line2.split("\\s");
	String strID1 = token1[0];
	String strID2 = token2[0];
	if (strID1.compareTo(strID2) != 0) {
	    System.err.println("different STR ID = "+strID1+" "+strID2);
	    System.exit(-1);
	}
	ArrayList<Integer> nonNAPosList = new ArrayList<Integer>();
	findNonNAPos(token1, token2, nonNAPosList, strID1);

	Hashtable<Integer,Integer> alleleInfo1 = new Hashtable<Integer,Integer>();
	TreeSet<Integer> alleleSet1 = new TreeSet<Integer>();
	Hashtable<Integer,Integer> alleleInfo2 = new Hashtable<Integer,Integer>();
	TreeSet<Integer> alleleSet2 = new TreeSet<Integer>();
	findNewAlleleCoding(token1, nonNAPosList, alleleInfo1, alleleSet1);
	findNewAlleleCoding(token2, nonNAPosList, alleleInfo2, alleleSet2);

	if (alleleInfo1.size() != alleleInfo2.size()) {
	    System.out.println(strID1+": #_of_alleless_in_old_STR = "+alleleInfo1.size()+" is NOT the same as #_of_alleles_in_new_STR = "+alleleInfo2.size());
	    inferAlleleCode(alleleSet1, alleleInfo1, alleleSet2, alleleInfo2);
	}
	else {
	    System.out.println(strID1+": #_of_alleless_in_old_STR = "+alleleInfo1.size()+" is the same as #_of_alleles_in_new_STR = "+alleleInfo2.size());
	}

	convertToNewCode(token1, nonNAPosList, alleleInfo1, outputList);
	convertToNewCode(token2, nonNAPosList, alleleInfo2, outputList);

	oldout.print(strID1);
	Iterator<Integer> ite1 = alleleSet1.iterator();
	while (ite1.hasNext()) {
	    int allele = ite1.next();
	    int coding = alleleInfo1.get(allele);
	    oldout.print("\t"+allele+":"+coding);
	}
	oldout.println();

	newout.print(strID1);
	Iterator<Integer> ite2 = alleleSet2.iterator();
	while (ite2.hasNext()) {
	    int allele = ite2.next();
	    int coding = alleleInfo2.get(allele);
	    newout.print("\t"+allele+":"+coding);
	}
	newout.println();
    }

    public void createFiles(String oldstrfile, String newstrfile, String outputoldfile, String outputnewfile, String newalleleoldfile, String newallelenewfile) throws Exception {
	BufferedReader in1 = new BufferedReader(new FileReader(oldstrfile));
	BufferedReader in2 = new BufferedReader(new FileReader(newstrfile));
	PrintWriter out1 = new PrintWriter(new BufferedWriter(new FileWriter(outputoldfile)));
	PrintWriter out2 = new PrintWriter(new BufferedWriter(new FileWriter(outputnewfile)));
	PrintWriter out3 = new PrintWriter(new BufferedWriter(new FileWriter(newalleleoldfile)));
	PrintWriter out4 = new PrintWriter(new BufferedWriter(new FileWriter(newallelenewfile)));
	String line1 = in1.readLine();
	String line2 = in2.readLine();
	line1 = in1.readLine();
	line2 = in2.readLine();
	while (line1 != null) {
	    ArrayList<String> outputList = new ArrayList<String>();
	    recodeTwoLines(line1, line2, outputList, out3, out4);
	    out1.println(outputList.get(0));
	    out2.println(outputList.get(1));
	    line1 = in1.readLine();
	    line2 = in2.readLine();
	}
	in1.close();
	in2.close();
	out1.close();
	out2.close();
	out3.close();
	out4.close();
    }

    public static void main(String[] args) throws Exception {
	if (args.length != 6) {
	    System.err.println("[usage] java RecodeBothOldNewSTR [old str file] [new str file] [output old str file] [output new str file] [output - allele code info file for old STR] [output - allele code info file for new STR]");
	    System.exit(-1);
	}
	RecodeBothOldNewSTR cc = new RecodeBothOldNewSTR();
	cc.createFiles(args[0], args[1], args[2], args[3], args[4], args[5]);
    }
}

