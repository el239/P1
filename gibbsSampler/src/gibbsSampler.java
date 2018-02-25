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
		  while(s.hasNextLine()) {
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
	   for (i = 0; i < t; i++) {
		   motifs[i] = randMotif(sequences[i]);
		   System.out.println(motifs[i]);
	   } // end for
	   
	   arraySlice(motifs);
	   motifCount(motifs);
	   
   } // end sampler
   
   public static String[] arraySlice(String[] allMotifs){
	   int rand = (int)(Math.random() * (allMotifs.length - 1));
	   String[] oneLess = new String[sequenceNumber - 1];
	   for ()
//	   System.arraycopy(allMotifs, 0, oneLessllMotifs.length - 1 - rand);
	   System.out.println("motif of sequence: " + rand + " removed.");
	   return allMotifs;
   } // end stringSlice
   
   public static int[][] motifCount(String[] theMotifs){
	   int i;
	   int j;
	   int aCount = 0;
	   int cCount = 0;
	   int gCount = 0;
	   int tCount = 0;
	   int [][] theMotifCount = new int[4][motifLength];
	   String letter = new String("");
	   for (i = 0; i < motifLength; i++) {
		   for (j = 0; j < sequenceNumber; j++) {
			   char base = theMotifs[j].charAt(i);
			   letter = ("");
			   letter += base;
			   switch (letter){
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
				   System.out.println("Invalid base letter in motif");
				   break;
			   } // end switch
		   } // end for

		   theMotifCount[0][i] = aCount;
		   theMotifCount[1][i] = cCount; 
		   theMotifCount[2][i] = gCount;
		   theMotifCount[3][i] = tCount;	
		   aCount = 0;
		   cCount = 0;	
		   gCount = 0;
		   tCount = 0;		   
	   } // end for

	   for (int[] motifLength : theMotifCount){
		    System.out.println(Arrays.toString(motifLength));
	   } // end for
	   
	   return theMotifCount;
   } // end motifCount

   public static String randMotif(String sequence) {
	   int start = (int)(Math.random() * (sequence.length() - motifLength + 1)); // find random element between 0 and near String end to begin
	   
	   return sequence.substring(start,start + motifLength); // creates random motif of k elements
   } // end randMotif

   
   
  
   
/*   
   public static void consensus(File fastaInput){

	   int i = 0;
	   try {
	   Scanner s = new Scanner(fastaInput);
		  while(s.hasNextLine()) {
			 lineCount ++;
			 input = s.nextLine();
			 if (lineCount % 2 == 0) {
				 System.out.println(lineCount);
				 cleanInput[i] = input;
				 i++;
			 } // end if
			 System.out.println(cleanInput[0]);
		  } // end while	
	   } // end try
	   catch (FileNotFoundException e) {
		e.printStackTrace();
	   } // end catch	   	   
	   
   } // end consensus
*/
} // end class













/*	
public static String mostfrequent(String text, int k){ 

count = 0;
String word, block = null, mostWord = new String("");
String finalWord = new String("");

if (k > text.length() || k <= 0) // precondition violated
	throw new IllegalArgumentException("Sequence length must be an integer between 1 and text length");

for (int p = 0; p < text.length()-k+1; p++){ // increments through k number of "slices" of first k-size sequence 
	word = text.substring(p,k+p); // takes sequence as a template
	int innerCount = 0;
	
	for (int i = 0; i < text.length()-k+1; i++){
       block = text.substring(i,i+k); // chops text in chunks of size k
       if (block.equals(word)){
    	   
    	   innerCount ++; // increments when "word" finds matches, min. value 1 for when word and block coincide
    	   i += k-1; // VERY IMPORTANT- skips to end of block, i.e., prevents "aaaa" being counted as 3 instances when k=2
           
    	   if (innerCount > count){  
               count = innerCount; // writes over count when a higher frequency sequence is encountered
               mostWord = block; // writes over the the most frequent "word"
               if (finalWord.length() == 0 && count > 1){ // sets initial final word
            	   finalWord = mostWord;
               } // end if
           } // end if
    		
    	    // checks for repeats in "finalWord" before appending
    		if (innerCount == count && innerCount > 1 && !block.equals(finalWord.substring(finalWord.length()-k, finalWord.length()))){
    		   int j;
    		   int appendCount = 0; // marks repeats
    		   for (j = 0; j < finalWord.length()/k; j++) {
    			   if (block.equals(finalWord.substring(j*k,j*k+k))) { // checks for repeats in final word
    				   appendCount++;
    		       } // end if
    		   } // end for
    		   
    			   if (appendCount == 0){
    	 	       finalWord += block; // appends tie cases, if not already in "finalWord"
    			   } // end if
    		       appendCount = 0;
   		    } // end if   
       } // end if
	} // end for 
	
innerCount = 0; // resets match counter
    
} // end for

if (count == 1) // case for no-repeat sequences
	return ("The sequence has no repeat elements.\n");

if (finalWord.length() == k) {
	System.out.print("The highest frequency sequence is " + "\"" + finalWord + "\"" + ", occurring " + count + " times." + "\n");
    return "";
    } // end if
else {
	System.out.print("The highest frequency sequences are ");
	for (int i = 0; i < finalWord.length()/k - 1; i++){ // print formatting
	   System.out.print("\"" + finalWord.substring(i*k,i*k+k) + "\"" + ", ");
	} // end for
	
	if (finalWord.length()/k > 2){ // "each"
	   System.out.print("and " + "\"" + finalWord.substring(finalWord.length()-k,finalWord.length()) + "\"" + ", each occurring " + count + " times." + "\n");
	return(""); 
	} // end if 
	
	else // "both"
	   System.out.print("and " + "\"" + finalWord.substring(finalWord.length()-k,finalWord.length()) + "\"" + ", both occurring " + count + " times." + "\n");
	return("");
} // end else
} // end mostfrequent
*/



