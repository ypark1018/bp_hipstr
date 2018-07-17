import java.io.*;
import java.util.*;
import java.lang.*;
import java.util.zip.*;

public class ConvertGIGIOutputVCF {

    public void createFiles(String gigifile, String origvcffile, String outputfile) throws Exception {
	BufferedReader ingeno = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gigifile))));

	BufferedReader invcf = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(origvcffile))));
	String linevcf = invcf.readLine();
	while (!linevcf.startsWith("#CHROM")) {
	    linevcf = invcf.readLine();
	}
	PrintWriter out = new PrintWriter(new GZIPOutputStream(new FileOutputStream(outputfile)));
	out.println("##fileformat=VCFv4.2");
	out.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	out.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	String linetemp = ingeno.readLine();
	if (!linetemp.startsWith("id")) {
	    System.err.println("genofile does not start with id ");
	    System.exit(-1);
	}
	String[] tokentemp = linetemp.split("\\s");
	for (int j = 1; j < tokentemp.length; j=j+2) {
	    out.print("\t"+tokentemp[j]);
	}
	out.println();
	ArrayList<String> lineList = new ArrayList<String>();
	lineList.add(ingeno.readLine());
	String line = lineList.get(0);
	linevcf = invcf.readLine();
	while (line != null) {
	    String[] token = line.split("\\s");
	    String[] tokenvcf = linevcf.split("\\s");
	    String[] postoken = token[0].split(":");
	    String chr = postoken[0];
	    int pos = Integer.parseInt(postoken[1]);
	    int posvcf = Integer.parseInt(tokenvcf[1]);
	    if (pos != posvcf) {
		System.err.println("positions do not match "+pos+" "+posvcf+" in vcf ");
		System.exit(-1);
	    }
            String refallele = tokenvcf[3];
            String altallele = tokenvcf[4];
	    out.print(chr+"\t"+pos+"\t"+chr+":"+pos+"\t"+refallele+"\t"+altallele+"\t"+tokenvcf[5]+"\t"+tokenvcf[6]+"\t"+tokenvcf[7]+"\tGT");
	    writeGeno(lineList, out);
	    lineList.set(0,ingeno.readLine());
	    line = lineList.get(0);
	    linevcf = invcf.readLine();
	}
	invcf.close();
	ingeno.close();
	out.close();
    }

    public void writeGeno(ArrayList<String> lineList, PrintWriter out) throws Exception {
	for (int i = 0; i < lineList.size(); i++) {
	    String line = lineList.get(i);
	    String[] token = line.split("\\s");
	    for (int j = 1; j < token.length; j=j+2) {
		String fallele = token[j];
		String sallele = token[j+1];
		String nfallele = "";
		boolean ismissing = false;
		if (fallele.compareTo("0") == 0) {
		    nfallele = ".";
		    ismissing = true;
		}
		else {
		    nfallele = (new Integer(Integer.parseInt(fallele)-1)).toString();
		}
		String nsallele = "";
                if (sallele.compareTo("0") == 0) {
		    nsallele = ".";
		    ismissing = true;
		}
                else {
		    nsallele = (new Integer(Integer.parseInt(sallele)-1)).toString();
                }
		if (ismissing) {
		    nfallele = ".";
		    nsallele = ".";
		}
		out.print("\t"+nfallele+"/"+nsallele);
	    }
	}
	out.println();
    }

    public static void main(String[] args) throws Exception {
	if (args.length != 3) {
	    System.err.println("[usage] java ConvertGIGIOutputVCF [gigi file] [orig vcf file] [output file]");
	    System.exit(-1);
	}
	ConvertGIGIOutputVCF cc = new ConvertGIGIOutputVCF();
	cc.createFiles(args[0], args[1], args[2]);
    }
}