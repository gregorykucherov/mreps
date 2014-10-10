/* mainrepets.c : computing repetitions of the first type (exact case) */

#include <stdlib.h>
#include <stdio.h>
#include "defs.h"

#define MINDIFLENPER       9
extern int allowsmall;
extern int *computeLIx (char *pattern, int patLen, int dir, int maxPerBnd);
extern int *computeLOx (char *pattern, int patLen, int dir, int maxPerBnd);
/* extern void NewFirstTypeRepeat(int newRepBgn, int newRepEnd); */
extern listreps *presortLists;
/* extern listreps *sortedLists; */
extern int maxPer;
extern char *seq_original;
/* extern float minExp; */
extern limits lim;
extern int step, LastWindow;

extern FILE *output_file;                 // output file (for xml output)
extern char * output_file_name;
extern int xmloutput;


void NewFirstTypeExactRepeat(int prRepBgn, int prRepEnd);
void startrepets(int curFacBgn, int nxtFacBgn);
void finalrepets (int prvFacBgn, int curFacBgn);

int Period;                            // declare it here! (rather than in searchforHeadGDR.c)
int prvFacLen, maxLftLen, curFacLen, maxLftRoot;


void mainrepets (int prvFacBgn, int curFacBgn, int nxtFacBgn)
     /* searching for exact reps of the first type crossing the 
	current border between factors in Crochemore's s-factorisation,
	and located within the current window;
	called by maxreps */
{ 
  int rightLen;
  int *LPu, *LSt,  *LPt,  *LSu;

  curFacLen=nxtFacBgn-curFacBgn;
  prvFacLen=curFacBgn-prvFacBgn;
  maxLftRoot=curFacLen+prvFacLen;
  if ( (maxLftLen=maxLftRoot+prvFacLen)>curFacBgn )
    {
      maxLftLen = curFacBgn ;
      if (maxLftRoot>curFacBgn)
	{  
	  maxLftRoot = curFacBgn+1;
	  if (maxLftRoot<=curFacLen)
	    {  
	      startrepets(curFacBgn, nxtFacBgn);
	      return;  
	    }
	}
     
    }

  VERB(1) printf(" s-factor : [prvFacLen=%d, maxLftLen=%d, maxLftRoot=%d, curFacLen=%d]\n", 
		 prvFacLen, maxLftLen, maxLftRoot, curFacLen ) ;

  /* computing Main's extension functions */
  LPu = computeLIx (seq_original+curFacBgn, curFacLen, 1, curFacLen);
  LSt = computeLIx (seq_original+curFacBgn-1, maxLftLen, -1, maxLftRoot);
 
  LPt = computeLOx (seq_original+curFacBgn, curFacLen, 1, maxLftRoot);
  LSu = computeLOx (seq_original+curFacBgn-1, maxLftLen, -1, curFacLen);

  Period=lim.min_period;
  while( Period<curFacLen && Period<=maxPer )    // LIMIT
    {
      /* reps with a full period in previous factors */
      VERB(3)
	printf("  Period=%d *2* :  LPt[%d]=%d, LSt[%d]=%d \n", Period, 
	       Period, LPt[Period], Period, LSt[Period]);
	    
      if (LPt[Period] + LSt[Period] >= Period)
	{  
	  if ( (rightLen=LPt[Period]) >= Period)
	    {  
	      if ( (rightLen<curFacLen) ||
		   (seq_original[nxtFacBgn] != seq_original[nxtFacBgn-Period]) )
		NewFirstTypeExactRepeat( curFacBgn-(LSt[Period]+Period), 
					 curFacBgn+rightLen );
	      Period++;
	      continue;  
	    }
	  else
	    {  
	      if ( ( rightLen>0 ) ||
		   ( LSt[Period]+Period<prvFacLen ) )
		NewFirstTypeExactRepeat( curFacBgn-(LSt[Period]+Period), 
					 curFacBgn+rightLen );
	    }
	}

      /* m-repets with a full period in current factor */
      VERB(3) 
	printf("  Period=%d *1* :  LSu[%d]=%d, LPu[%d]=%d \n", Period ,
	       Period, LSu[Period], Period, LPu[Period]) ;

      if (LSu[Period] + LPu[Period] >= Period)
	{  
	  if (  (LPu[Period]+Period<curFacLen) ||
		(seq_original[nxtFacBgn] != seq_original[nxtFacBgn-Period]) )
	    NewFirstTypeExactRepeat( curFacBgn-LSu[Period],
				     curFacBgn+LPu[Period]+Period );
	}
      Period++;
    }

  while( Period<maxLftRoot && Period<=maxPer )    // LIMIT
    {
      /* reps with a full period in previous factors only */
      VERB(3)
	printf("  Period=%d *2* :  LPt[%d]=%d, LSt[%d]=%d \n", Period, 
	       Period, LPt[Period], Period, LSt[Period]);
	    
      if (LPt[Period] + LSt[Period] >= Period)
	{  
	  if ( (rightLen=LPt[Period])<curFacLen )
	    {  
	      if ( ( rightLen>0 ) ||
		   ( LSt[Period]+Period<prvFacLen ) )
		NewFirstTypeExactRepeat( curFacBgn-(LSt[Period]+Period), 
					 curFacBgn+rightLen );
	    }
	  else
	    {  
	      if (seq_original[nxtFacBgn] != seq_original[nxtFacBgn-Period])
		NewFirstTypeExactRepeat( curFacBgn-(LSt[Period]+Period), 
					 curFacBgn+rightLen );
	    }
	}
      Period++;
    }
  free (LPu);
  free (LSu);
  free (LPt);
  free (LSt);  
}



void NewFirstTypeExactRepeat(int newRepBgn, int newRepEnd)
     /* called from mainrepets, startrepets and finalrepets;
	processing an exact rep of the first type. For a found rep,
	creates the corresponding structure 'listreps' and 
	inserts it into the corresponing list from 'presortLists'
	(array of lists of reps of the first type found in the 
	current window and sorted by end positions) */
{ 
  listreps newListRep, nextRep;
  repetition *newRep;
  int rsize ; 

  rsize = newRepEnd-newRepBgn;         // the size of the rep 

  if ( ((rsize-Period>=MINDIFLENPER) || allowsmall) &&                          
       (step<0 || newRepBgn < VIRTUAL_STEP)                        // GK
     )             
    {  
      if ((nextRep=presortLists[newRepEnd]))
	{  
	  if (nextRep->rep.initpos == newRepBgn)
	    {  
	      VERB(3) printf (": *NO* \n") ;
	      return;
	    }
	}
      if ( (newListRep=(listreps)malloc(sizeof(struct s_listreps)))==NULL )
	{  
	  fprintf(stderr,"The output is too big.\n");
	  fprintf(stderr,"Narrow period bounds.\n");
	  RECOVER_XML_OUTPUT(output_file,output_file_name,20);
	  exit(20);
	}
      newRep=&(newListRep->rep);
      newRep->initpos=newRepBgn;
      newRep->endpos=newRepEnd;
      newRep->period=Period;
      newListRep->next=nextRep;
      presortLists[newRepEnd]=newListRep;
      VERB(3) printf (": *OK* \n") ;
    }
  else
    VERB(3) printf ("\n");
}
 
 
void startrepets(int curFacBgn, int nxtFacBgn)
     /* called from mainrepets;
	similar to 'mainrepets' for the case when the reps can
	contain the beginning of the current window */
{  
  int rightLen;
  int *LPu, *LSt,  *LPt,  *LSu;

  VERB(1) printf(" s-factor : [prvFacLen=%d, maxLftLen=%d, maxLftRoot=%d, curFacLen=%d]\n", 
		 prvFacLen, maxLftLen, maxLftRoot, curFacLen ) ;

  /* computing longest extension functions */
  LPu = computeLIx (seq_original+curFacBgn, curFacLen, 1, curFacLen);
  LSt = computeLIx (seq_original+curFacBgn-1, maxLftLen, -1, maxLftRoot);
 
  LPt = computeLOx (seq_original+curFacBgn, curFacLen, 1, maxLftRoot);
  LSu = computeLOx (seq_original+curFacBgn-1, maxLftLen, -1, curFacLen);

  Period=lim.min_period;
  while( Period<maxLftRoot && Period<=maxPer )    // LIMIT
    {
      /* reps with a full period in previous factors */
      VERB(3) printf("  Period=%d *2* :  LPt[%d]=%d, LSt[%d]=%d \n", Period, 
		     Period, LPt[Period], Period, LSt[Period]);
	    
      if (LPt[Period] + LSt[Period] >= Period)
	{  
	  if ( (rightLen=LPt[Period]) >= Period)
	    {  
	      if ( (rightLen<curFacLen) ||
		   (seq_original[nxtFacBgn] != seq_original[nxtFacBgn-Period]) )
		NewFirstTypeExactRepeat( curFacBgn-(LSt[Period]+Period), 
					 curFacBgn+rightLen );
	      Period++;
	      continue;  
	    }
	  else
	    {  
	      if ( ( rightLen>0 ) ||
		   ( LSt[Period]+Period<prvFacLen ) )
		NewFirstTypeExactRepeat( curFacBgn-(LSt[Period]+Period), 
					 curFacBgn+rightLen );
	    }
	}

      /* reps with a full period in current factor */
      VERB(3) 
	printf("  Period=%d *1* :  LSu[%d]=%d, LPu[%d]=%d \n", Period ,
	       Period, LSu[Period], Period, LPu[Period]) ;

      if (LSu[Period] + LPu[Period] >= Period)
	{  
	  if (  (LPu[Period]+Period<curFacLen) ||
		(seq_original[nxtFacBgn] != seq_original[nxtFacBgn-Period]) )
	    NewFirstTypeExactRepeat( curFacBgn-LSu[Period],
				     curFacBgn+LPu[Period]+Period );
	}
      Period++;
    }

  while( Period<curFacLen && Period<=maxPer )    // LIMIT
    {
      /* reps with a full period in current factor only */
      VERB(3)
	printf("  Period=%d *1* :  LSu[%d]=%d, LPu[%d]=%d \n", Period ,
	       Period, LSu[Period], Period, LPu[Period]) ;

      if (LSu[Period] + LPu[Period] >= Period)
	{  
	  if (  (LPu[Period]+Period<curFacLen) ||
		(seq_original[nxtFacBgn] != seq_original[nxtFacBgn-Period]) )
	    NewFirstTypeExactRepeat( curFacBgn-LSu[Period],
				     curFacBgn+LPu[Period]+Period );
	}
      Period++;
    }
  free (LPu);
  free (LSu);
  free (LPt);
  free (LSt);  
}

     
void finalrepets (int prvFacBgn, int curFacBgn)
     /* called from maxreps; 
	searching for exact reps of the first type 
	containing the end of the current window, 
	but not the last border between s-factors (?)
	(in other words, all exact reps of the first 
	type in the end of word, which 
	have not been found by mainrepets) */
{ 
  int *LSt ;
  prvFacLen=curFacBgn-prvFacBgn;
  maxLftRoot=(1+prvFacLen)/2;

  VERB(1) printf(" final repetitions : [prvFacLen=%d, maxLftRoot=%d, curFacLen=%d]\n", 
		 prvFacLen, maxLftRoot, curFacLen ) ;

  /* computing longest extension functions */
  LSt = computeLIx (seq_original+curFacBgn-1, prvFacLen, -1, maxLftRoot);

  Period=lim.min_period;
  while( Period<maxLftRoot && Period<=maxPer )    // LIMIT
    {
      /* reps in previous factor completely */
      VERB(3)
	printf("  Period=%d *2* :  LSt[%d]=%d \n", Period, Period, LSt[Period]);
	    
      if (LSt[Period] >= Period)
	{  
	  if ( LSt[Period]+Period<prvFacLen ) 
	    NewFirstTypeExactRepeat( curFacBgn-(LSt[Period]+Period), 
				     curFacBgn );
	}
      Period++;
    }
  free (LSt);  
}
