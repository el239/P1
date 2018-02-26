import java.io.*;
import java.util.Scanner;
import java.util.*;


public class gibbsSampler{

public static int lineCount = 0;
public static int sequenceNumber;
public static int motifLength;
public static int sequenceLength;
public static String input;

public static void main(String args[]){
	
	   /* to test runtime
	   long startTime = System.nanoTime();
	   String testRunTime = (mostfrequent("variablestringlength",2));
	   long endTime   = System.nanoTime();
	   long totalTime = endTime - startTime;
	   System.out.println(totalTime+ " nanoseconds");
	   */
//	   randMotif("afdjkdsdf",4);
	   sampler("10by100.FASTA", 3, 10, 100);

//	   consensus(text);
} // end main

   public static void sampler(String DNA, int k, int t, int N) {
	   sequenceNumber = t;
	   motifLength = k;
	   sequenceLength = N;
   
   File text = new File(DNA); // file must be present in same directory or given with specified file path
	   String[] sequences = new String[t];
	   String[] motifs = new String[t];
	   try {
	      Scanner s = new Scanner(text);
		  while(s.hasNextLine()) { // cleans FASTA input and generates sequence array
			 lineCount ++;
			 input = s.nextLine();
			 if(lineCount % 2 == 0){
			 sequences[lineCount/2 - 1] = input;
			 } // end if
		  } // end while
		  
	   } // end try
	   catch (FileNotFoundException e) {
		e.printStackTrace();
	   } // end catch	   
	   
	   int i;
	   for (i = 0; i < t; i++) { // creates random motifs of size k
		   motifs[i] = randMotif(sequences[i]);
		   System.out.println(motifs[i]);
	   } // end for
	   
	   String[] motifsAfterRemoval = new String[sequenceNumber - 1]; // new string with memory allocated for once motif less
	   motifsAfterRemoval = arraySlice(motifs); // removes random motif from array
	   
	   int[][] countMatrix = new int[4][motifsAfterRemoval.length]; // creates two-dimensional array to hold count
	   countMatrix = motifCount(motifsAfterRemoval); // makes motif count from spliced array
	   countMatrix = pseudoCount(countMatrix); // overwrites count with pseudocount
   } // end sampler
   
   public static int[][] pseudoCount(int[][] theCountMatrix){
	   int i;
	   int j;
	   int x = 0;
	   for (i = 0; i < 4; i++){ // iterates though elements in OUTSIDE loop
		   for (j = 0; j < motifLength; j++) { // iterates through Strings in INSIDE loop
			   theCountMatrix[i][j] = theCountMatrix[i][j] + 1; // adds 1 to every count
		   } // end for
	   } // end for
	   System.out.println("pseudocount: ");
	   for (int[] motifLength : theCountMatrix){ // nice matrix print
		    switch (x){ // format print
		    case 0:
		    	System.out.print("a: ");
		    	break;
		    case 1:
		    	System.out.print("c: ");
		    	break;
		    case 2:
		    	System.out.print("g: ");
		    	break;
		    case 3:
		    	System.out.print("t: ");
		    	break;
		    } // end switch
		    System.out.println(Arrays.toString(motifLength));
		    x++;
	   }// end for
	   return theCountMatrix;
   } // end pseudoCount
   
   public static String[] arraySlice(String[] allMotifs){ // cuts out a random element motif
	   int i;
	   int j = 0;
	   int rand = (int)(Math.random() * (allMotifs.length - 1));
	   String[] oneLess = new String[sequenceNumber - 1]; // creates string with memory allocated for 1 less motif
	   System.out.println("Motif of sequence corresponding to element #" + rand + " removed.");
	   System.out.println("Motifs after splice: ");
	   for (i = 0; i < rand; i++){ // fills new array with portion before splice 
		   oneLess[i] = allMotifs[i];
		   System.out.println(oneLess[i]);
	   } // end for
	   for (j = i; j < oneLess.length; j++) { // fills new array with portion after splice
		   oneLess[j] = allMotifs[j + 1];
		   System.out.println(oneLess[j]);
	   } // end for
	   return oneLess;
   } // end stringSlice
   
   public static int[][] motifCount(String[] theMotifs){
	   int i;
	   int j;
	   int x = 0;
	   int aCount = 0;
	   int cCount = 0;
	   int gCount = 0;
	   int tCount = 0;
	   int [][] theMotifCount = new int[4][motifLength];
	   String letter = new String("");
	   for (i = 0; i < motifLength; i++){ // iterates though elements in OUTSIDE loop
		   for (j = 0; j < theMotifs.length; j++) { // iterates through Strings in INSIDE loop
			   char base = theMotifs[j].charAt(i);
			   letter = ("");
			   letter += base;
			   switch (letter){ // tallies nucleotides at each motif position 
			   case "a":
				   aCount++;
				   break;
			   case "c":
				   cCount++;
				   break;
			   case "g":
				   gCount++;
				   break;
			   case "t":
				   tCount++;
				   break;
			   default:
				   System.out.println("Invalid base letter in motif"); // shouldn't occur unless there is an unrecognized letter
				   break;
			   } // end switch
		   } // end for
		   
		   // fills two-dimensional array with counts
		   theMotifCount[0][i] = aCount;
		   theMotifCount[1][i] = cCount; 
		   theMotifCount[2][i] = gCount;
		   theMotifCount[3][i] = tCount;	
		   // resets counts for next element position
		   aCount = 0;
		   cCount = 0;	
		   gCount = 0;
		   tCount = 0;		   
	   } // end for
       System.out.println("count: ");
	   for (int[] motifLength : theMotifCount){ // nice matrix print
		    switch (x){ // format print
		    case 0:
		    	System.out.print("a: ");
		    	break;
		    case 1:
		    	System.out.print("c: ");
		    	break;
		    case 2:
		    	System.out.print("g: ");
		    	break;
		    case 3:
		    	System.out.print("t: ");
		    	break;
		    } // end switch
		    System.out.println(Arrays.toString(motifLength));
		    x++;
	   } // end for
	   
	   return theMotifCount;
   } // end motifCount

   public static String randMotif(String sequence) {
	   int start = (int)(Math.random() * (sequence.length() - motifLength + 1)); // find random element between 0 and near String end to begin
	   return sequence.substring(start,start + motifLength); // creates random motif of k elements
   } // end randMotif

   
   
  
   
/*   
   public static void consensus(File fastaInput){

   } // end consensus
*/
} // end class
