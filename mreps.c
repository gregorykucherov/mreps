/* Computing exact repetitions */

#include <stdlib.h>
#include <stdio.h>
#include "defs.h"

extern void print_rep(int rinitpos, int rendpos, int rlength, int rperiod, int rnumerr);
extern void print_score(int rinitpos,int rendpos,int rlength, int rperiod, float rscore);
extern void finalrepets (int prvFacBgn, int curFacBgn);
extern void SortFirstTypeRepeats(void);
extern void EmergSortFirstTypeRepeats(void);  
extern void FindSecondTypeExactRepeats(void);
void  Showmreps(void);
void print_mrepeat(int rinitpos, int rendpos, int rperiod);
extern void RecodeSmalSeq(int chkLength);
extern void RecodeCapSeq(int chkLength);
extern int *computeLIx(char *pattern, int patLen, int dir, int maxPerBnd);
extern void mainrepets(int prvFacBgn, int curFacBgn, int nxtFacBgn);

extern sfactorization  FactorizeforGDR(void);

extern void guessDepth(int length);

extern int step, LastWindow;
extern int maxPer, oldPerBnd;
extern char *seq_factor_copy;
extern int FACTYPE;
extern listreps *presortLists, *sortedLists, *maxIntPtr;
extern char  nextLetterCheck, prevLetterCheck ;
extern limits lim ;

extern FILE *output_file;                 // output file (for xml output)
extern char * output_file_name;
extern int xmloutput;

int lenWrd;    // declare it here!
int tooBigReps ;

char *seq_original, *seq_original1;   // treated window of original sequence

int verbose = 1;         // defined as extern in defs.h
int verbose_long = 1 ;   // defined as extern in defs.h

void  Savemreps(void);
int check_mrepeat(listreps checkRepeat);
void ComputePerBound(void);
extern int numFact, *factorBgns;
extern listreps *allReps;


void maxreps (char *seq_parameter, int length_parameter)
     /* seaching for exact reps in the current window */
{
  int i;  
  sfactorization sfact ;
  seq_original=seq_parameter; // remember the sequence address in global variable 'seq_original'
  seq_original1=seq_original-1;
  lenWrd=length_parameter;    // remember the sequence length in global variable 'lenWrd'

  VERB(2) printf("input : %s\n", seq_original); 
  VERB(2) printf("length : %d\n", length_parameter);

  /* validating maximal possible period */ 
  if (LastWindow)
     maxPer=(lim.max_period>length_parameter/2) ? length_parameter/2 :lim.max_period;
  else
     maxPer=(lim.max_period>step/2) ? step/2 :lim.max_period;
  VERB(1)  printf("Maximal possible period of repeats is %d\n\n", maxPer);

  /* s-factorization (sfact.c) */

  // create a copy of the sequence to process it in RecodeSeq()
  seq_factor_copy = (char *)malloc(length_parameter+2);
  seq_factor_copy[0]='\0';
/*   strcpy (seq_factor_copy+1, seq_original); */

  if (*seq_original>='a')
     RecodeSmalSeq(length_parameter);
  else
     RecodeCapSeq(length_parameter);

  VERB(2) 
    {
      printf("recoded seq: ");
      for (i=1;i<=length_parameter;i++) printf("%1d",seq_factor_copy[i]);
      printf("\n");
    }

  FACTYPE=WITHOUT_EXTRA_SYMBOL;
  guessDepth(length_parameter);
  sfact = FactorizeforGDR();
  
  if (lim.err_number)
     ComputePerBound();
  else
     free(seq_factor_copy);
 
  VERB(1) printf("s-factorization : %d factors\n", sfact->nbfactors);
  
  /* repets of type 1 (mainrepets.c) */

  if( ( presortLists=(listreps *) calloc ((2+lenWrd), sizeof(listreps)) ) == NULL )
    {  
      fprintf(stderr,"Not enough memory for saving the repetitions\n");
      RECOVER_XML_OUTPUT(output_file,output_file_name,20);
      exit(20);  
    }

  for (i=1; i< sfact -> nbfactors ; i++)
    {

      VERB(1) printf(" s-fact: [%d, %d]\n", 
                     sfact->factor[i], sfact->factor[i+1]) ;
      mainrepets (sfact->factor[i-1], sfact->factor[i], sfact->factor[i+1]); 
	  
    }       
  finalrepets (sfact->factor[i-1], sfact->factor[i]);

  VERB(1) printf("repetitions of type 1 : OK.\n");
  
  /* Now we sort repetitions of type 1 */

/*  if ( (sortedLists=(listreps *)calloc(lenWrd+2, sizeof(listreps))) )
 *   SortFirstTypeRepeats();
 *  else
 */
    EmergSortFirstTypeRepeats();  
 
  VERB(1) printf("sorting m-repets of type 1 : OK.\n");

  VERB(2) printf ("now repetitions of type 2 \n");

  /* repetitions of type 2 */

  FindSecondTypeExactRepeats();

  VERB(1) printf("finding repetitions of type 2 : OK.\n");

  free (sfact);

  /* Now we print the appropriate repetitions */
  if (lim.err_number==0)
     Showmreps();
  else
     Savemreps();
}


void  Showmreps(void)
     /* filtering out exact reps which start in the second half of the window (unless it is the last one),
	or extend beyond the window;
	call print_mrepeat on those reps which go through filtering */
{ 
  int endPstn, initPos, repPer;
  listreps current_rep, *current_pointer, remListRep;

/*   if (sortedLists+step<=maxIntPtr) */
/*   if (step > 0 && sortedLists+step<=maxIntPtr)     // if it is not the last window (?) */
  if (LastWindow==NO)
    maxIntPtr=sortedLists+step-1;                  // output repetition starting in the first half only

  current_pointer=sortedLists;


  if ( (current_rep=(*current_pointer)) )  // if there are reps starting at the first position
                                  // treat them separately (because of prevLetterCheck)
    do
      {  
	repPer=current_rep->rep.period;
	if (seq_original1[repPer]!=prevLetterCheck) // repeat does not extend to the left
	  {  
	    if ( (endPstn=current_rep->rep.endpos)<lenWrd ) // it does not touch the end
	      print_mrepeat(0, endPstn, repPer);
	    else
	      {  
		if (seq_original[lenWrd-repPer]!=nextLetterCheck)
		  print_mrepeat(0, lenWrd, repPer);
		else
		  tooBigReps=YES;
	      }
	  } 
	current_rep=(remListRep=current_rep)->next;
	free(remListRep); 
      }  
    while(current_rep);

  while(++current_pointer<=maxIntPtr)     
    {  
      if ((current_rep=(*current_pointer)))
	{  
	  initPos=current_rep->rep.initpos;
	  do
	    {  
	      repPer=current_rep->rep.period;
	      if ( (endPstn=current_rep->rep.endpos)<lenWrd ) // the repeat does not
		// touch the end of the window
		print_mrepeat(initPos, endPstn, repPer);
	      else
		{  
		  if (seq_original[lenWrd-repPer]!=nextLetterCheck)
		    print_mrepeat(initPos, lenWrd, repPer);
		  else
		    tooBigReps=YES;
		}
	      current_rep=(remListRep=current_rep)->next;
	      free(remListRep); 
	    }  
	  while(current_rep);
	}
    }

  free(sortedLists);
}

void  Savemreps(void)
     /* filtering out exact reps which start in the second half of the window (unless it is the last one),
	or extend beyond the window;
	save those reps which go through filtering */
{ 
  listreps current_rep, *current_pointer, next_rep;
  
  if ( (allReps=(listreps *)calloc(oldPerBnd, sizeof(listreps)))==NULL )
    {  
      printf("Error: Not enough memory for storing exact repeats\n");
      RECOVER_XML_OUTPUT(output_file,output_file_name,20);
      exit(20);
    }

  if (LastWindow==NO)
    maxIntPtr=sortedLists+step-1;                  // output repetition starting in the first half only

  for(current_pointer=maxIntPtr;current_pointer>sortedLists;current_pointer--)   
    {  
      if ((current_rep=(*current_pointer)))
	{ while ( (next_rep=current_rep->next) )
	  {  if ( check_mrepeat(current_rep) )
	     {  current_rep->next=allReps[current_rep->rep.period];
                allReps[current_rep->rep.period]=current_rep;
	     }
             else
	        free(current_rep);
	     current_rep=next_rep;
	  }
          if ( (current_rep->rep.endpos<lenWrd) ||              // the repeat does not
               (seq_original[lenWrd-current_rep->rep.period]    // cross the end of the window
                != nextLetterCheck) )
	  {  if ( check_mrepeat(current_rep) )
	     {  current_rep->next=allReps[current_rep->rep.period];
                allReps[current_rep->rep.period]=current_rep;
	     }
             else
	        free(current_rep);
	  }
	  else
          {  free(current_rep);
             tooBigReps=YES;
          }
	}  
    }

  if ( (current_rep=(*current_pointer)) )  // if there are reps starting at the first position
                                           // treat them separately (because of prevLetterCheck)
  {  while ( (next_rep=current_rep->next) )
     {  if (seq_original1[current_rep->rep.period]!=prevLetterCheck) 
                                           // repeat does not extend to the left
	{  if ( check_mrepeat(current_rep) )
	   {  current_rep->next=allReps[current_rep->rep.period];
              allReps[current_rep->rep.period]=current_rep;
	   }
           else
	      free(current_rep);
	}
        else
	   free(current_rep);
        current_rep=next_rep;
     }
     if (seq_original1[current_rep->rep.period]!=prevLetterCheck) 
                                           // repeat does not extend to the left
     {  if ( (current_rep->rep.endpos<lenWrd) ||              // the repeat does not
             (seq_original[lenWrd-current_rep->rep.period]    // cross the end of the window
              != nextLetterCheck) )
        {  if ( check_mrepeat(current_rep) )
	   {  current_rep->next=allReps[current_rep->rep.period];
              allReps[current_rep->rep.period]=current_rep;
	   }
           else
	      free(current_rep);
        }
        else
        {  free(current_rep);
           tooBigReps=YES;
        }
     }
     else
	free(current_rep);
  }

  free(sortedLists);
}


int check_mrepeat(listreps checkRep)
     /* check if the period of the given exact rep is primitive and pass it further for 
        saving in structure "allReps" */
{ int rinitpos, rendpos;
  int rlength, j, *LP;
  
  rinitpos=checkRep->rep.initpos++;
  rendpos=checkRep->rep.endpos++;
  rlength=rendpos-rinitpos;

  if (lim.min_period>1)
  {  
/*   compute the minimal period */
     LP=computeLIx(seq_original+rinitpos, rlength, 1, lim.min_period);
     for (j = 1 ; j<lim.min_period ; j++)
     {  
        if ( j+LP[j]>=rlength )
        { free(LP);       // throw away the repetition if its minimal period
          return(NO); }   // is smaller than the user-specified minimal period
     }
     free(LP);
  }

  return(YES);
}     


void ComputePerBound()
{   int fctIndx, sumFact;
    if (numFact<3)
    {
      printf("Error: number of errors is too big for this sequence\n");
      printf("       (or the sequence has too regular structure)\n");
      printf("       Run with a smaller number of errors\n");
      RECOVER_XML_OUTPUT(output_file,output_file_name,71);
      exit(71);
    }
    oldPerBnd=factorBgns[fctIndx=2];
    while(++fctIndx<=numFact)
    {  if ( (sumFact=factorBgns[fctIndx]-factorBgns[fctIndx-2])>oldPerBnd )
          oldPerBnd=sumFact;
    }
}


void print_mrepeat(int rinitpos, int rendpos, int rperiod)
     /* check if the period of the given exact rep is primitive and pass it further for printing 
	to function 'print_rep' (for exe version) or 'printWeb' (for web version)
	(both defined in file printOutput.c) */
{  
  int rlength, j, *LP;
  float rexp;	 

/*   compute the minimal period */
  rlength=rendpos-rinitpos;               // the length of the repeat
  if (lim.min_period>1)
    {  
      LP=computeLIx(seq_original+rinitpos, rlength, 1, lim.min_period);
      for (j = 1 ; j<lim.min_period ; j++)
	{  
	  if ( j+LP[j]>=rlength )
	  { free(LP);   // throw away the repetition if its minimal period
            return; }   // is smaller than the user-specified minimal period
	}
      free(LP);
    }

  rexp = (float)rlength/rperiod;          // the exp of the repeat     
    if ( (rlength >= lim.min_size) && (rlength <= lim.max_size)    // LIMIT for size
         && (rexp >= lim.min_exponent) )                           // LIMIT for exponent

  //This part is for the exe version of the program. Comment it out for web version.

/*        print_rep(rinitpos+1, rendpos, rlength,  rperiod, 0); */
       print_score(rinitpos+1, rendpos, rlength,  rperiod, 0.0);


  //This part is for the web version which should be compiled on a SUN machine. Comment it out if used for the exe version.
  //printWeb(start_pstn, rinitpos, rendpos, rlength, rperiod)
   
}

