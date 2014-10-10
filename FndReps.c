/* computing approximate repetitions */

/* #include "FndReps.h" */

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include <sys/types.h>
#include <sys/stat.h>
#define TAKEALLNUMERR
#define MIN_PER_COEF       2
#define MINDIFLENPER       9
#define CHECKLENBASE      18
#define NUMERRCOEFLEN      4
#define CHECKSCOREBASE  0.24
#define SMALLSCOREBASE  0.12
#define NUMERRCOEFSCOR  0.01
#define MAXCHECKSCORE    0.4
#define SHOWSCORE
#define NUMERRCOEF       1.2
#define SHOWJUNK 
extern limits lim ;
extern int allowsmall;
int endPstn, initPos, repPeriod, repLength, repNumErr;
extern char *seq_working_copy, *seq_original, *seq_original1, *seq_factor_copy;
extern int lenWrd, actPerBnd, maxNumErr, maxNumErr2;
/* int *numRep=0;  */         // was never initialized ?? [GK]
extern int maxPer, minPer, dblMinPer, oldPerBnd  ;
/* extern float minExp; */
extern listreps *sortedLists, *perLists, *initLists, *maxIntPtr; 
/* extern int PRIMEDEPTH; */
extern int FACTYPE, numFact;
extern int start_pstn;       // window start position
/* extern int final_pstn; */
/* int numRun , numRun2 ; */
/* extern int noprint; */
extern int step, LastWindow, maxInitPstn;
listreps *initLists1, *allReps;
extern FILE *output_file;                 // output file (for xml output)
extern char * output_file_name;
extern int xmloutput;
extern int *copyEnds, *factorEnds;
extern void maxreps (char *seq_parameter, int length_parameter);
extern void guessDepth(int);
extern void print_rep(int,int,int,int,int);
extern void print_score(int,int,int,int,float);
extern sfactorization FactorizeforGDR(void);
extern void ComputeGlobalFactors(void);
extern void mainstage(void);
extern void SortFirstTypeRepeats(void);
extern void EmergSortFirstTypeRepeats(void);
extern void FindSecondTypeRepeats(void);
extern void SortAllRepeats(void);
extern void EmergSortAllRepeats(void);
extern void FindAllRuns(void);
extern void FindAllStepRuns(void);
void Per2InitLists(void);
void RecodeSmalSeq(int chkLength);
void RecodeCapSeq(int chkLength);
void showByInit(void);
int FindNumErr(void);
void CheckMrepsPeriod(void);
void CheckMrepScore(void);
int CheckRun1(listreps checkedRun);
void AddNewReps(void);
float FindScore(void);
int FixNumErrProgram(void);
int FiltReps(listreps);
float repScore, checkScore;
//extern int numBadReps, numGoodReps;
float maxCheckLen, minCheckLen, midCheckLen, chkLenCoef, smalLenCoef;
#ifdef SHOWJUNK
extern int showjunk ;
extern FILE *junkfile;
#endif


void mismatchProgram(char *seq_parameter, int length_parameter) 
     /* computing approximate reps in the current window */
{

  maxreps (seq_parameter, length_parameter);

  seq_working_copy = seq_parameter-1; // for compatibility with main stage of algorithm
  /* validating maximal possible period */ 
  if (LastWindow)
    maxPer=length_parameter/2;
  else
    maxPer=step/2;

  VERB(1)  printf("Maximal possible period of repeats is %d\n\n", maxPer);

  guessDepth(length_parameter);
  FACTYPE=WITH_EXTRA_SYMBOL;
  FactorizeforGDR();
  free(seq_factor_copy);

#ifdef TAKEALLNUMERR
  for(maxNumErr=1; maxNumErr<=lim.err_number; maxNumErr++)
    if ( FixNumErrProgram()==NO )
       break;
#else
  maxNumErr=1;
  while(maxNumErr < lim.err_number)
    {  
      if ( FixNumErrProgram()==NO )
	{  
	  maxNumErr=0;
	  break;
	}
      if (maxNumErr<(int)(NUMERRCOEF*maxNumErr))
        maxNumErr=(int)(NUMERRCOEF*maxNumErr);
      else
        maxNumErr++;
    }
  if (maxNumErr)
    {  
      maxNumErr=lim.err_number;
      FixNumErrProgram();  
    }
#endif
  free(copyEnds);
  free(factorEnds);

  Per2InitLists();
  showByInit();  
}

int FixNumErrProgram(void)
{ 
#ifndef  TAKEALLNUMERR
  printf("maxNumErr=%d\n", maxNumErr);
#endif 
  
  minPer=((int)(MIN_PER_COEF*maxNumErr)>maxNumErr) ? 
    (int)(MIN_PER_COEF*maxNumErr) : maxNumErr+1;
  minPer=(minPer>lim.min_period) ? minPer : lim.min_period;
  dblMinPer=2*minPer;
  maxNumErr2=maxNumErr+2;
  if ( (lenWrd<dblMinPer) || (numFact<maxNumErr2) )
    {
/*       printf("Error: number of errors is too big for this sequence\n"); */
/*       printf("       (or the sequence has too regular structure)\n"); */
/*       printf("       Run with a smaller number of errors\n"); */
      printf("Warning: resolution is too big for this sequence\n");
      return(NO);
    }
  ComputeGlobalFactors();
  mainstage();

  EmergSortFirstTypeRepeats();  
      
  FindSecondTypeRepeats();
  
  if ( (perLists=(listreps *)calloc(actPerBnd, sizeof(listreps))) )
    SortAllRepeats();
  else
    EmergSortAllRepeats();
  
  /*   if (step>=lenWrd) */
  /*   if (step < 0 || VIRTUAL_STEP >= lenWrd)  // GK */
      /*      FindAllRuns(); */
      /*   else */
      /*      FindAllStepRuns() ; */
      if (LastWindow)
	FindAllRuns();
      else
	FindAllStepRuns() ;

      CheckMrepsPeriod();
      CheckMrepScore();

      AddNewReps();
      return(YES); 
}

/* void howMany ( ){ */
/*  printf("\n***Total number of global runs is %d ***\n\n", numRun2); */
/* } */


void RecodeCapSeq(int checkLen)
/* recode sequence 'seq_factor_copy' starting from position 1 (!)
   set the global length variable */
{
  char curSmbl;
  int pstn=1;
  int number_of_N=0;
  
  while((curSmbl=seq_original1[pstn])!='\0')
    { 
      switch (curSmbl) 
	{
	case 'a':
          seq_original1[pstn]='A';
	case 'A':
	  seq_factor_copy[pstn++]=1;
	  break;
	case 'c':
          seq_original1[pstn]='C';
	case 'C':
	  seq_factor_copy[pstn++]=2;
	  break;
        case 'g':
          seq_original1[pstn]='G';
	case 'G':
	  seq_factor_copy[pstn++]=3;
	  break;
	case 't':
          seq_original1[pstn]='T';
	case 'T':
	  seq_factor_copy[pstn++]=4;
	  break;
	  /* 	commented out by GK, 17/10/01 */
	  /*       case '\n': */
	  /* 	  seq_factor_copy[pstn]=0; */
	  /* 	  break; */
	case 'n':
	case 'N':
	  {
	    number_of_N++;
	    if (1.00*number_of_N / checkLen > LIMIT_N_PROPORTION)
	      {
		fprintf(stderr,"Error: Too many N's in the window\n");
		RECOVER_XML_OUTPUT(output_file,output_file_name,17)
		exit (17);
	      }
	    seq_factor_copy[pstn]=(rand() % 4)+1;
            switch (seq_factor_copy[pstn])
	      {  case 1:
                   seq_original1[pstn++]='A';
                   break;
                 case 2:
                   seq_original1[pstn++]='C';
                   break;
                 case 3:
                   seq_original1[pstn++]='G';
                   break;
                 case 4:
                   seq_original1[pstn++]='T';
                   break;
	      }
/* 	    printf("zamenyaem na %d\n",seq_factor_copy[pstn-1]); */
	    break;
	  }
	default :
	  {
	    fprintf(stdout, "Error: input sequence has incorrect symbol %c at position %d\n", 
		    curSmbl, start_pstn+pstn);
	    RECOVER_XML_OUTPUT(output_file,output_file_name,16)
	    exit(16);
	  }
	}
    }
/*   lenWrd=pstn-1; */
  if (pstn-1 != checkLen)
    {  
      fprintf(stdout, "Error: input sequence has special zero symbol at position %d\n",start_pstn+pstn);
      exit(-9);
    }
  seq_factor_copy[pstn]='\0';

  if (number_of_N)
    { printf ("Warning: symbols N in your sequence has been replaced randomly\n"); 
      printf ("by {A,C,G,T}. This might have created artefact repetitions.\n\n"); 
    }
  VERB(3)  printf("Length of input window : %d\n", lenWrd);
}


void RecodeSmalSeq(int checkLen)
/* recode sequence 'seq_factor_copy' starting from position 1 (!)
   set the global length variable */
{
  char curSmbl;
  int pstn=1;
  int number_of_N=0;
  
  while((curSmbl=seq_original1[pstn])!='\0')
    { 
      switch (curSmbl) 
	{
	case 'A':
          seq_original1[pstn]='a';
	case 'a':
	  seq_factor_copy[pstn++]=1;
	  break;
	case 'C':
          seq_original1[pstn]='c';
	case 'c':
	  seq_factor_copy[pstn++]=2;
	  break;
        case 'G':
          seq_original1[pstn]='g';
	case 'g':
	  seq_factor_copy[pstn++]=3;
	  break;
	case 'T':
          seq_original1[pstn]='t';
	case 't':
	  seq_factor_copy[pstn++]=4;
	  break;
	  /* 	commented out by GK, 17/10/01 */
	  /*       case '\n': */
	  /* 	  seq_factor_copy[pstn]=0; */
	  /* 	  break; */
	case 'n':
	case 'N':
	  {
	    number_of_N++;
	    if (1.00*number_of_N / checkLen > LIMIT_N_PROPORTION)
	      {
		fprintf(stderr,"Error: Too many N's in the window\n");
		RECOVER_XML_OUTPUT(output_file,output_file_name,17)
		exit (17);
	      }
	    seq_factor_copy[pstn]=(rand() % 4)+1;
            switch (seq_factor_copy[pstn])
	      {  case 1:
                   seq_original1[pstn++]='a';
                   break;
                 case 2:
                   seq_original1[pstn++]='c';
                   break;
                 case 3:
                   seq_original1[pstn++]='g';
                   break;
                 case 4:
                   seq_original1[pstn++]='t';
                   break;
	      }
/* 	    printf("zamenyaem na %d\n",seq_factor_copy[pstn-1]); */
	    break;
	  }
	default :
	  {
	    fprintf(stdout, "Error: input sequence has incorrect symbol %c at position %d\n", 
		    curSmbl, start_pstn+pstn);
	    RECOVER_XML_OUTPUT(output_file,output_file_name,16)
	    exit(16);
	  }
	}
    }
/*   lenWrd=pstn-1; */
  if (pstn-1 != checkLen)
    {  
      fprintf(stdout, "Error: input sequence has special zero symbol at position %d\n",start_pstn+pstn);
      exit(-9);
    }
  seq_factor_copy[pstn]='\0';

  if (number_of_N)
    { printf ("Warning: symbols n in your sequence has been replaced randomly\n"); 
      printf ("by {a,c,g,t}. This might have created artefact repetitions.\n\n"); 
    }
  VERB(3)  printf("Length of input window : %d\n", lenWrd);
}



void Per2InitLists(void)
{  
  listreps nxtrep, current;
  if (LastWindow)
    maxInitPstn=lenWrd-lim.min_period;
  else
    maxInitPstn=step+oldPerBnd;
  if ( (initLists=(listreps *)calloc(maxInitPstn, sizeof(listreps)))==NULL )
    {  
      printf("Error: Not enough memory for storing repeats\n");
      RECOVER_XML_OUTPUT(output_file,output_file_name,20)
      exit(20);
    }
  initLists1=initLists-1;
  for(repPeriod=oldPerBnd-1; repPeriod>=lim.min_period; repPeriod--)
    {   
      if ((current=allReps[repPeriod]))
	{  
	  do
	    {  
	      nxtrep=current->next;
	      current->next=initLists1[current->rep.initpos];
	      initLists1[current->rep.initpos]=current;
	    }  
	  while((current=nxtrep));
	}
    }
  free(allReps);
}

/* Output results according to initial position */

void showByInit(void) 
     /* processing the runs stored in the array initLists.
	the runs verifying all parameters are output by function
	'print_rep' */
{  
  listreps current, remListRep;

  for(initPos=1;initPos<=maxInitPstn;initPos++)
    {  
      if ((current=initLists1[initPos]))
	{  
	  do
	    {  
	      endPstn=current->rep.endpos;
	      repLength=endPstn-initPos;
	      repPeriod=current->rep.period;
	      if ( (repLength>=lim.min_size) && (repLength<=lim.max_size) && 
		   ((float)repLength/repPeriod>=lim.min_exponent) )
		{ 	     
#ifdef SHOWSCORE
		  repScore=FindScore();
		  print_score(initPos, endPstn-1, repLength, repPeriod, repScore);
#else
		  repNumErr=FindNumErr();
		  print_rep(initPos, endPstn-1, repLength, repPeriod, repNumErr);
#endif
		}

	      current=(remListRep=current)->next;
	      free(remListRep); 
	    }  
	  while(current);
	}
    }
  free(initLists);
}

#ifndef SHOWSCORE
int FindNumErr(               )
{   
  char *frstSmb, *scndSmb, *thrdSmb, *endPtr;
  int nmbErr=0, maxNmbErr;
  scndSmb=seq_original1+initPos;
  frstSmb=scndSmb+repPeriod;
  if (repLength<=2*repPeriod)
    {  
      endPtr=seq_original1+endPstn;
      while(frstSmb<endPtr)
	{  
	  if ( *(frstSmb++)!=*(scndSmb++) )
	    nmbErr++;  
	}
      return(nmbErr);
    }
  else
    {  
      endPtr=frstSmb+repPeriod;
      while(frstSmb<endPtr)
	{  
	  if ( *(frstSmb++)!=*(scndSmb++) )
	    nmbErr++;  
	}
      maxNmbErr=nmbErr;
      endPtr=seq_original1+endPstn;
      thrdSmb=seq_original1+initPos;
      do
	{  
	  if ( *(frstSmb++)!=*scndSmb )
	    {  
	      if ( *(scndSmb++)==*(thrdSmb++) )
		{  
		  if ( ++nmbErr>maxNmbErr )
		    maxNmbErr=nmbErr;  
		}
	    }
          else
	    {  
	      if ( *(scndSmb++)!=*(thrdSmb++) )
                nmbErr--;
	    }
	}  
      while(frstSmb<endPtr);
      return(maxNmbErr);
    }
}
#endif

void CheckMrepsPeriod(void)
{  listreps current, prevrep, remListRep;
   int realPer;
#ifdef SHOWJUNK
   int pstn;
#endif
   initLists1=initLists-1;
   actPerBnd=(actPerBnd>lim.max_period) ? lim.max_period+1 : actPerBnd ;
   if ( (perLists=(listreps *)calloc(actPerBnd, sizeof(listreps)))==NULL )
   {  printf("Not enough memory for checking repeats\n");
      exit(44);
   }

   for(initPos=1;initPos<=maxInitPstn;initPos++)
   {  
      if ((current=initLists1[initPos]))
      {  do
	{   realPer=CheckRun1(current);
            if ( (realPer>=lim.min_period) && (realPer<=lim.max_period) )
            {  if ((prevrep=perLists[realPer]))
               {  if (endPstn>prevrep->rep.endpos)
		  {  if ( (initPos<=prevrep->rep.endpos-2*realPer)
                          || (initPos==prevrep->rep.initpos) )
                     {  prevrep->rep.endpos=endPstn;
                        current=(remListRep=current)->next;
	                free(remListRep);
                     }
                     else
                     {  current->rep.initpos=initPos;
                        current->rep.period=realPer;
                        current=(remListRep=current)->next;
                        remListRep->next=prevrep;
                        perLists[realPer]=remListRep;
                     }
                  }
                  else
                  {  current=(remListRep=current)->next;
	             free(remListRep);
                  }
               }
               else
               {  current->rep.initpos=initPos;
                  current->rep.period=realPer;
                  current=(remListRep=current)->next;
                  remListRep->next=NULL;
                  perLists[realPer]=remListRep;
               }
            }
            else
            {  
#ifdef SHOWJUNK
	      if(showjunk)
              {  for(pstn=initPos;pstn<endPstn;pstn++)
                    fputc(seq_original1[pstn], junkfile);
                 fputc('\n', junkfile);
                 fprintf(junkfile, "is of period %d\n\n", realPer);
	      }
#endif               
               current=(remListRep=current)->next;
	       free(remListRep);
            }
         }  while(current);
      }
   }
   free(initLists);
}


void CheckMrepScore(void)
{  listreps current, prevrep, nextrep;

   minCheckLen=MINDIFLENPER+1;
   midCheckLen=(float)(CHECKLENBASE+NUMERRCOEFLEN*maxNumErr);
   checkScore=CHECKSCOREBASE+NUMERRCOEFSCOR*maxNumErr;
   chkLenCoef=checkScore/(midCheckLen-minCheckLen);
   maxCheckLen=minCheckLen+MAXCHECKSCORE/chkLenCoef;
   if (allowsmall)
      smalLenCoef=(float)SMALLSCOREBASE/midCheckLen;

   for(repPeriod=actPerBnd-1; repPeriod>=lim.min_period; repPeriod--)
   {  if ((prevrep=perLists[repPeriod]))
      {  if ((current=prevrep->next))
         {  prevrep->next=NULL;
            while((nextrep=current->next))
            {  if ( FiltReps(prevrep) )
                  current->next=prevrep;
               else
               {  current->next=prevrep->next;
                  free(prevrep);
               }
               prevrep=current;
               current=nextrep;
            }
            if ( FiltReps(prevrep) )
               current->next=prevrep;
            else
            {  current->next=prevrep->next;
               free(prevrep);
            }
            if ( FiltReps(current) )
               perLists[repPeriod]=current;
            else
            {  perLists[repPeriod]=current->next;
               free(current);
            }
         }
         else
         {  if ( FiltReps(prevrep) )
               perLists[repPeriod]=prevrep;
            else
            {  perLists[repPeriod]=NULL;
               free(prevrep);
            }
         }
      }
   }
}


void AddNewReps(void)
{  listreps *newReps, nxtrep, fstrep, scdrep;
   int newPerBnd, dblPer;
   if (oldPerBnd>=actPerBnd)
   {  newReps=perLists;
      newPerBnd=actPerBnd;
   }
   else
   {  newReps=allReps;
      allReps=perLists;
      newPerBnd=oldPerBnd;
      oldPerBnd=actPerBnd;
   }
   for(repPeriod=newPerBnd-1; repPeriod>=lim.min_period; repPeriod--)
   {  if ((fstrep=newReps[repPeriod]))
      {  if ((scdrep=allReps[repPeriod]))
         {  dblPer=2*repPeriod;
            if ( fstrep->rep.initpos < scdrep->rep.initpos )
               allReps[repPeriod]=fstrep;
            else
	    {  fstrep=scdrep;
               scdrep=newReps[repPeriod];
            }
            while((nxtrep=fstrep->next))
            {  if ( nxtrep->rep.initpos <= scdrep->rep.initpos )
                  fstrep=nxtrep;
               else
               {  if ( scdrep->rep.endpos>fstrep->rep.endpos )
                  {  if ( (scdrep->rep.initpos<=fstrep->rep.endpos-dblPer)
                          || (scdrep->rep.initpos==fstrep->rep.initpos) )
                     {  fstrep->rep.endpos=scdrep->rep.endpos;
                        fstrep->next=scdrep->next;
                        free(scdrep);
                        scdrep=nxtrep;
                     }
                     else
                     {  fstrep->next=scdrep;
                        fstrep=scdrep;
                        scdrep=nxtrep;
                     }
                  }
                  else
                  {  nxtrep=scdrep->next;
                     free(scdrep);
                     if ( (scdrep=nxtrep)==NULL )
                        break;
                  }
               }
            }
            if (scdrep)
            {  while ( scdrep->rep.endpos<=fstrep->rep.endpos )
               {  nxtrep=scdrep->next;
                  free(scdrep);
                  if ( (scdrep=nxtrep)==NULL )
                     break;
               }
               if (scdrep)
               {  if ( (scdrep->rep.initpos<=fstrep->rep.endpos-dblPer)
                       || (scdrep->rep.initpos==fstrep->rep.initpos) )  
                  {  fstrep->rep.endpos=scdrep->rep.endpos;
                     fstrep->next=scdrep->next;
                     free(scdrep);
                  }
                  else
                     fstrep->next=scdrep;
               }
            }
         }
         else
            allReps[repPeriod]=fstrep;                 
      }
   }
   free(newReps);
}


int FiltReps(listreps actrep)
{   float limScore;
#ifdef SHOWJUNK
    int pstn;
#endif
    initPos=actrep->rep.initpos;
    endPstn=actrep->rep.endpos;
    repLength=endPstn-initPos;
    if (allowsmall)
    {  repScore=FindScore();
       if (repLength>=midCheckLen)
       {  if (repLength<=maxCheckLen)
             limScore=chkLenCoef*(repLength-minCheckLen);
          else
             limScore=MAXCHECKSCORE;
       }
       else
       {  if ( (limScore=checkScore-(midCheckLen-repLength)*smalLenCoef)
                >MAXCHECKSCORE ) limScore=MAXCHECKSCORE;
       }
    }
    else
    {  if (repLength-repPeriod>=MINDIFLENPER)
       {  repScore=FindScore();
          if (repLength<=maxCheckLen)
	     limScore=chkLenCoef*(repLength-minCheckLen);
          else
             limScore=MAXCHECKSCORE;
       }
       else
       { 
#ifdef SHOWJUNK
          if (showjunk)
          {  for(pstn=initPos;pstn<endPstn;pstn++)
                fputc(seq_original1[pstn], junkfile);
             fputc('\n', junkfile);
             fprintf(junkfile, "is of period %d and length %d\n\n", 
                               repPeriod, repLength);
          }
#endif
          return(NO);
       } 
    }          
    if (repScore<=limScore)
       return(YES);
    else
    {
#ifdef SHOWJUNK
       if (showjunk)
       {  for(pstn=initPos;pstn<endPstn;pstn++)
             fputc(seq_original1[pstn], junkfile);
          fputc('\n', junkfile);
          fprintf(junkfile, "is of period %d and relative score %f\n\n", 
                            repPeriod, limScore-repScore);
       }
#endif  
       return(NO);
    }
     
}


float FindScore(void)
{  
  int *chkErr, *chkPtr, realNumErr=0;
  char *frstSmb, *scndSmb, *endPtr;
  chkErr=(int *)malloc(repLength*sizeof(int));
  scndSmb=seq_original1+initPos;
  frstSmb=scndSmb+repPeriod;
  if (repLength<=2*repPeriod)
    {  
      endPtr=seq_original1+endPstn;
      do
	{  
	  if ( *(frstSmb++)!=*(scndSmb++) )
            realNumErr++;  
	}  
      while(frstSmb<endPtr);
    }
  else
    {  
      chkPtr=chkErr;
      endPtr=frstSmb+repPeriod;
      do
	{  
	  if ( *(frstSmb++)!=*(scndSmb++) )
	    {  
	      realNumErr++;
	      *(chkPtr++)=0;
	    }
	  else
            chkPtr++;
	}  
      while(frstSmb<endPtr);
      chkPtr=chkErr;
      endPtr=seq_original1+endPstn;
      do
	{  
	  if ( *frstSmb!=*scndSmb )
	    {  
	      if ( (*frstSmb!=*(scndSmb-repPeriod)) || *chkPtr )
		{  
		  realNumErr++;
		  *(chkPtr+repPeriod)=0;
		}
	      else
		*(chkPtr+repPeriod)=1;
	    }
	  chkPtr++;
	  scndSmb++;
	}  
      while(++frstSmb<endPtr);
    }
  free(chkErr);
  return( (float)realNumErr/(repLength-repPeriod) );
}


int CheckRun1(listreps chkRun)
{  int chkPer, realPer, numErr, *chkErr, *chkPtr, realNumErr=0;
   char *frstSmb, *scndSmb, *endPtr;
   endPstn=chkRun->rep.endpos;
   repLength=endPstn-initPos;
   realPer=chkRun->rep.period;
   chkErr=(int *)malloc(repLength*sizeof(int));
   scndSmb=seq_original1+initPos;
   frstSmb=scndSmb+realPer;
   if (repLength<=2*realPer)
   {  endPtr=seq_original1+endPstn;
      do
      {  if ( *(frstSmb++)!=*(scndSmb++) )
            realNumErr++;  
      }  while(frstSmb<endPtr);
   }
   else
   {  chkPtr=chkErr;
      endPtr=frstSmb+realPer;
      do
      {  if ( *(frstSmb++)!=*(scndSmb++) )
         {  realNumErr++;
            *(chkPtr++)=0;
         }
         else
            chkPtr++;
      }  while(frstSmb<endPtr);
      chkPtr=chkErr;
      endPtr=seq_original1+endPstn;
      do
      {  if ( *frstSmb!=*scndSmb )
         {  if ( (*frstSmb!=*(scndSmb-realPer)) || *chkPtr )
            {  realNumErr++;
               *(chkPtr+realPer)=0;
            }
            else
               *(chkPtr+realPer)=1;
         }
         chkPtr++;
         scndSmb++;
      }  while(++frstSmb<endPtr);
   }
   chkPer=(realPer>repLength/2) ? repLength/2 : realPer-1;
   do
   {  numErr=0;
      scndSmb=seq_original1+initPos;
      frstSmb=scndSmb+chkPer;
      chkPtr=chkErr;
      endPtr=frstSmb+chkPer;
      do
      {  if ( *(frstSmb++)!=*(scndSmb++) )
         {  numErr++;
            *(chkPtr++)=0;
         }
         else
            chkPtr++;
      }  while(frstSmb<endPtr);
      chkPtr=chkErr;
      endPtr=seq_original1+endPstn;
      while(frstSmb<endPtr)
      {  if ( *frstSmb!=*scndSmb )
         {  if ( (*frstSmb!=*(scndSmb-chkPer)) || *chkPtr )
            {  numErr++;
               *(chkPtr+chkPer)=0;
            }
            else
               *(chkPtr+chkPer)=1;
         }
         chkPtr++;
         frstSmb++;
         scndSmb++;
      }
      if (realNumErr)
      {  if ( (numErr*(repLength-realPer))/(repLength-chkPer)<=realNumErr )
         {  realPer=chkPer;
            realNumErr=numErr;
         }
      }
      else
      {  if (numErr==0)
            realPer=chkPer;
      }
   }  while(--chkPer>=1);
   free(chkErr);
   return(realPer);
}


