import java.io.*;
import java.util.Scanner;
import java.util.*;

public class gibbsSampler{

public static int lineCount = 0;
public static int sequenceNumber;
public static int motifLength;
public static int sequenceLength;
public static int rand;
public static String input;
public static String deletedSequence;

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
	   float[][] profileMatrix = new float[4][motifsAfterRemoval.length];
	   profileMatrix = profileCalc(countMatrix);
	   float[] probabilities = new float[sequenceLength - motifLength + 1]; // an element for each possible k-mer
	   deletedSequence = sequences[rand]; // instantiates sequence deleted from beginning of run
	   probabilities = probCalc(profileMatrix); // checks all possible k-mers in deleted strand to make probabilities
	   probabilities = dieMaker(probabilities); // adds a multiplier and overwrites the array so that the sum of elements equals 1
	   int newMotifStartElement = dieRoller(probabilities);
	   motifs[rand] = deletedSequence.substring(newMotifStartElement, newMotifStartElement + motifLength); // sets motif of delete sequence to die result
	   System.out.println("the new motif for the deleted sequence is: " + motifs[rand]);
	   
	   char[] consensusMotif = new char[motifLength];
	   consensusMotif = consensusCalc(motifCount(motifs)); // fills the consensus char array with the most frequent base

	   System.out.println("\nHamming distance score is: " + hammingCalc(motifs, consensusMotif));
   } // end sampler
   
   public static int hammingCalc(String[] theMotifs, char[] theConsensusMotif){
	   int i;
	   int j;
	   int hammingScore = 0;
	   for (i = 0; i < sequenceNumber; i++){ // iterates though elements in OUTSIDE loop
		   for (j = 0; j < motifLength; j++){ // iterates through Strings in INSIDE loop
			   if (theMotifs[i].charAt(j) != theConsensusMotif[j]){ // tallies for every mismatch
				   hammingScore++;
			   } // end if
		   } // end for
   	   } // end for
	   return hammingScore;
   } // end hammingCalc
   
   public static char[] consensusCalc(int[][] theMotifs) {
	   int i;
	   char[] theConsensusMotif = new char[motifLength]; 
	   System.out.println("Consensus sequence: ");
	   for (i = 0; i < motifLength; i++){ // iterates though elements in OUTSIDE loop
		   
		   int max = theMotifs[0][i]; // sets default "winner" base for consensus formation 
		   theConsensusMotif[i] = 'a'; // bias towards a<-c<-g<-t, but to no significant effect overal
		   if (theMotifs[1][i] > max) {
			   max = theMotifs[1][i]; // overwrites 'a' if there is a base with a higher tally
			   theConsensusMotif[i] = 'c'; 
			   } // end if
		   if (theMotifs[2][i] > max) {
			   max = theMotifs[2][i]; // overwrites 'a' if there is a base with a higher tally
			   theConsensusMotif[i] = 'g'; 
			   } // end if
		   if (theMotifs[3][i] > max) {
			   max = theMotifs[3][i]; // overwrites 'a' if there is a base with a higher tally
			   theConsensusMotif[i] = 't'; 
			   } // end if
		   System.out.print(theConsensusMotif[i]);
		   max = 0;
		   } // end for
	   return theConsensusMotif;
   } // end consensus Calc
   
   public static int dieRoller(float[] theProbabilities) {
	   float runningSum = 0;
	   int i;
	   float random = (float)Math.random();
	   System.out.println("\nRANDOM NUMBER IS: " + random);
	   for (i = 0; i < theProbabilities.length; i++){
		   runningSum += theProbabilities[i]; // sequentially adds up values to create an array of increasing floats
		   if (random <= runningSum){
			   System.out.println("Running sum: " + runningSum);
			   System.out.println("element: " + i);
			   return i;
		   } // end if
	   } // end for
	   System.out.println("Error selecting motif from deleted sequence (dieRoller function");
	   return -1;
   } // end theProbabilities
   
   public static float[] dieMaker(float[] theProbabilities) {
	   int i;
	   float sum = 0;
	   System.out.println("THE PROBABILITIES");
	   for(i = 0; i < theProbabilities.length; i++){
		   sum += theProbabilities[i]; // gets sum of all probabilities to make denominator with
	   } // end for
	   for(i = 0; i < theProbabilities.length; i++){
		   theProbabilities[i] = theProbabilities[i] / sum;
		   System.out.print(theProbabilities[i] + "  ");
	   } // end for
	   return theProbabilities;
   } // end dieMaker
   
   public static float[] probCalc(float[][] theProfileMatrix){
	   float[] theProbabilities = new float[sequenceLength - motifLength + 1];
	   int i;
	   int j;
	   int p;
	   float kmerProbability = 1;
	   float baseChance = 0;
	   String kmer = new String("");
	   for (p = 0; p < sequenceLength-motifLength + 1; p++){ // increments through k-size "slices" of deleted sequence 
			kmer = deletedSequence.substring(p, motifLength + p);
			for (i = 0; i < 3; i++) {
				if (kmer.charAt(i)=='a') {
					baseChance = theProfileMatrix[0][i];
				}
				else if (kmer.charAt(i)=='c') {
					baseChance = theProfileMatrix[1][i];
				}
				else if (kmer.charAt(i)=='g'){
					baseChance = theProfileMatrix[2][i];
				}
				else if (kmer.charAt(i)=='t') {
					baseChance = theProfileMatrix[3][i];
				}
				else {
					System.out.println("Error encountered fetching probabilities"); // debug safety
				} // end else
				kmerProbability *= baseChance;
				baseChance = 0;
			} // end for
			theProbabilities[p] = kmerProbability;
			System.out.print(theProbabilities[p] + "  ");
			kmerProbability = 1;
	   } // end for 
	   return theProbabilities;
   } // end probCalc
   
   public static float[][] profileCalc(int[][] theCountMatrix){
	   int i;
	   int j;
	   int x = 0;
	   float[][] theProfileMatrix = new float[4][motifLength];
	   for (i = 0; i < 4; i++){ // iterates though elements in OUTSIDE loop
		   for (j = 0; j < motifLength; j++){ // iterates through Strings in INSIDE loop
			   theProfileMatrix[i][j] = (float)theCountMatrix[i][j]/ (sequenceNumber - 1 + 4); // adjusts denominator
		   } // end for
   	   } // end for
	   System.out.println("profile: ");
	   for (float[] motifLength : theProfileMatrix){ // nice matrix print
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
	   return theProfileMatrix;
   } // end profileCalc

   public static int[][] pseudoCount(int[][] theCountMatrix){
	   int i;
	   int j;
	   int x = 0;
	   for (i = 0; i < 4; i++){ // iterates though elements in OUTSIDE loop
		   for (j = 0; j < motifLength; j++){ // iterates through Strings in INSIDE loop
			   theCountMatrix[i][j] = (theCountMatrix[i][j] + 1); // adds 1 to every count, divides by 8
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
	   rand = (int)(Math.random() * allMotifs.length - 1);
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

} // end class
