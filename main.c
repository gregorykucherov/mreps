/* main.c : Main interface */
/***************************/

#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <time.h>

#include "defs.h"
#define INPUT_ARG    0
#define INPUT_STDIN  1   // input from stdin; not implemented
#define INPUT_FILE   2

#define MAX_STEP        5000000
#define MAX_INPUT_SHOW  300
#define MAX_LINE_LENGTH 80

#define SHOWJUNK 
extern int tooBigReps;         // flag if there exist repetitions which don't fit to the window
int allowsmall = NO;    // for indicating small repetitions which may be aleatory

extern void maxreps(char *seq, int length);
extern void mismatchProgram(char *seq,int length);
extern void printNbOfReps(void);     // printing the number of repetitions found


/* Parameters */
limits lim ;
int from = 1 ;
int to = -1 ;
int toWasSpecified = NO;
int maxsizeWasSpecified = NO;
int maxperiodWasSpecified = NO;
int start_pstn;                    // window start position
int step = -1;                     // window size
int LastWindow = NO;               // flag indicating whether the window treated is the last one
int maxPer, minPer, dblMinPer, oldPerBnd, maxNumErr ;
extern int actPerBnd;              // this variables are to insure
                                   // interface with Roman's part
int noprint = 0;                   // flag whether reps themselves should be printed
FILE *output_file;                 // output file (for xml output)
char * output_file_name;
int xmloutput = NO;                // no xml output by default
time_t now;

/* Variables */
int inputLength;       // sequence length or file size
int nrRep ;            // number of repetitions output
int situation;         // situation code
int maxPer, minPer, dblMinPer, oldPerBnd  ;
char nextLetterCheck='\0', prevLetterCheck='\0';  // letters immediately following and preceding
                                                  // the current window
//int numBadReps, numGoodReps;
#ifdef SHOWJUNK
int showjunk = NO;
FILE *junkfile;
#endif

/* Little helps */

void wrongcall(char *prog)
{
  fprintf(stderr,"Usage: %s [ options ] { sequencefile | -s sequence }\n", prog);
  fprintf(stderr,"Try %s -help\n", prog);
  exit(1);
}

void toofew(char *prog)
{
  fprintf(stderr,"Too few arguments\n") ;
  wrongcall (prog) ;
}

void init_limits ()
{
  lim. err_number = 0 ;
  lim. min_period = lim. min_exponent = lim. min_size = 1 ;
  lim. max_period = lim. max_exponent = lim. max_size = -1 ;
}
    
void logo()
{
  printf("\n *****************************************************************************\n");
  printf(  " *                              mreps 2.6                                    *\n");
  printf(  " *                                                                           *\n");
  printf(  " *                Finding tandem repeats in DNA sequences                    *\n");
  printf(  " *                                                                           *\n");
  printf(  " *                      http://mreps.univ-mlv.fr/                            *\n");
  printf(  " *****************************************************************************\n\n");
}

void help()     
{ 
  
  printf( "Usage:\n");
  printf( "       mreps [ <options> ] { <sequencefile> | -s <sequence> }\n");
  printf( "finds tandemly repeated fragments in a DNA sequence \n\n");
  printf( "The options are :\n");
  //  printf( " -b : batch mode (verbose current state on stderr)\n");
/*   printf( " -v : verbose (and so -vv, -vvv)\n\n"); */
  printf( " -s <string>  : specifies the sequence in command line\n");
  printf( " -fasta       : allows DNA sequences in FASTA format \n\n");
  printf( " -res n       : \"resolution\" (error level)\n");
  printf( " -from n      : starting position n\n");
  printf( " -to n        : end position n\n");
  printf( " -minsize n   : repeats whose size is at least n\n");
  printf( " -maxsize n   : repeats whose size is at most n\n");
  printf( " -minperiod n : repeats whose period is at least n\n");
  printf( " -maxperiod n : repeats whose period is at most n\n");
  printf( " -exp x       : repeats whose exponent is at least x\n");
/*   printf( " -nacgt       : allow the occurrence of non-acgt \n"); */
/*   printf( "                characters in a DNA sequence.\n"); */
/*   printf( "                (all the non-acgt characters will be considered\n"); */
/*   printf( "                as the same character)\n"); */
/*   printf( "                (option used only for ACGT sequence or file)\n"); */
  printf( " -allowsmall  : output small repeats that can occur randomly\n\n");
/*   printf( "               (All separators will be considered as the\n"); */
/*   printf( "               same character: blank.)\n"); */
  printf( " -win n       : process by sliding windows of size 2*n overlaping by n\n");
  printf( " -xmloutput <file> : outputs to <file> in xml format\n");
  printf( " -noprint     : if specified, the repetition sequences will not be output \n");
  printf( "\nExample:\n");
  printf( "mreps -res 3 -exp 3.0 -from 10000 -to 12000 ecolim52.fas\n\n");
  exit(1);
}

/* Three macros for argument parsing  !!! Test errors !! */

#define ARG_INT(option,name)          \
{                                     \
  if (!strcmp (*argv, option))        \
    {                                 \
      argc-- ; argv++ ;               \
      if (argc == 1) toofew (prog) ;  \
      name = atoi(*argv);             \
      argc-- ; argv++ ;               \
      if (name <= 0)                  \
        {                             \
      fprintf(stderr,"Error: Argument %s (=%d) must be greater or equal to 1!\n",option, name);\
      exit (8);                     \
        }                             \
      continue;                       \
    }                                 \
}

#define ARG_FLOAT(option,name)        \
{                                     \
  if (!strcmp (*argv, option))        \
    {                                 \
      argc-- ; argv++ ;               \
      if (argc == 1) toofew (prog) ;  \
      name = atof(*argv);             \
      argc-- ; argv++ ;               \
      if (name <= 1.0)                \
        {                             \
      fprintf(stderr,"Error: Argument %s (=%f) must be greater than 1.0!\n",option, name);\
      exit (9);                    \
        }                             \
      continue;                       \
    }                                 \
}

#define ARG_FLAG(option,action)       \
{                                     \
  if (!strcmp (*argv, option))        \
    {                                 \
      action ;                        \
      argc-- ; argv++ ;               \
      continue;                       \
    }                                 \
}

#define min(a,b) ((a)>(b) ? b : a)
 
#define PROCESS_ACGT(word,length)     \
{ dblMinPer=2*lim.min_period;         \
  if (length>=dblMinPer)              \
  {  if (lim.err_number == 0)         \
       maxreps (word,length);         \
     else                             \
       mismatchProgram(word,length);  \
     printNbOfReps();                 \
  }                                   \
  else                                \
     printf("Processed sequence is too short\n\n");\
}

/* ************************************************************************* */

char RetrieveChar(int sit , FILE * inputFile)
/* retrieves the next character in the ACGT case */
/* (assumes sit==ACGT_FILE_FASTA or sit==ACGT_FILE_PLAIN) */
/* skips end-of-line */
{
  char c;
  while ((c=fgetc(inputFile))=='\n');
  return c;
}

void printNbOfReps()
{ 
  if (xmloutput==NO)
    if (nrRep)
      {
	printf(" ---------------------------------------------------------------------------------------------\n");
	if (nrRep == 1)
	  printf("RESULTS: There is 1 repeat in the processed sequence\n\n");
	else
	  printf("RESULTS: There are %d repeats in the processed sequence\n\n", nrRep);
      }
    else
      printf("RESULTS: There are no repeats in the processed sequence\n\n");
  else
    fprintf(output_file,"<nbofreps>%d</nbofreps>\n",nrRep);

    if (tooBigReps)
      {
	if (xmloutput==NO)
	  printf("Warning: repeats spanning beyond the window have been detected\n\n");
	else
	   fprintf(output_file,"<warning-big-reps-detected></warning-big-reps-detected>\n");
      }
}


/* ************************************************************************* */

char ProcessSeqFromFile(FILE * inputFile)
     /* The function processes the next sequence from the input file */
     /* returns the first character after the sequence */
     /* normally '>', if the sequence ends with '>' (case ACGT_FILE_FASTA) */
     /* or EOF, if there is nothing after */
{ 
  int seq_pstn=0;      // sequence position counter
  char * buff;          // MAKE IT GLOBAL
  char workchar;
  int counter, buffLen;

  prevLetterCheck = '\0';  // in case it is not the first sequence (fasta)
  nrRep = 0; 

/* Skipping sequence to FROM position */
  while(seq_pstn < from)
    {
      workchar = RetrieveChar(situation,inputFile);
      if (workchar==EOF || (situation==ACGT_FILE_FASTA && workchar=='>')) 
	{
	  printf("Warning: The sequence is shorter than FROM parameter (%d)\n\n", from+1); 
	  return workchar;
	}
      seq_pstn++;
    }

/* Allocate buffer */
  if (step > 0)
    buffLen=2*step;
  else 
    buffLen=min(inputLength,to-from);
  buff = malloc((3+buffLen)*sizeof(char));
  buff[0]=buff[1]='\0';
  buff+=2;
  start_pstn=seq_pstn;  // remember window start position in sequence

/* Read sequence after from */
  counter=0;
  workchar= RetrieveChar(situation,inputFile);
  while ((workchar != EOF) && (situation!=ACGT_FILE_FASTA || workchar!='>') && seq_pstn < to)
    {
      buff[counter++]=workchar;
      workchar=RetrieveChar(situation,inputFile);
      seq_pstn++;
      if (counter == buffLen && seq_pstn != to && workchar != EOF)  // buffer is full but TO is not achived 
	                                                            // and EOF is not reached
	                                                            // (can happen only with win option)
	{
	  buff[counter]='\0';
	  nextLetterCheck=workchar;
	  nrRep = 0; 
	  tooBigReps = NO;
	  LastWindow=NO;
	  /* Process current block */
	  if (xmloutput == NO)
	    printf("\n* Processing window [%d : %d] *\n\n", start_pstn+1, start_pstn+counter);
	  else 
	    {
	      fprintf(output_file,"<window>\n");
	      fprintf(output_file,"<windowstart>%d</windowstart>\n",start_pstn+1);
	      fprintf(output_file,"<windowend>%d</windowend>\n",start_pstn+counter);
	    }
	  PROCESS_ACGT(buff, counter);
	  if (xmloutput==YES)
	    fprintf(output_file,"</window>\n");

	  /* Shift the second half of the block onto the first */
	  prevLetterCheck=buff[step-1];
	  counter=0;
	  while (counter < step)
	    {
	      buff[counter]=buff[counter+step];
	      counter++;
	    }
	  start_pstn+=step;
	}
    }

  /* Processing tail block */
  if (xmloutput == NO)
    printf("\n* Processing window [%d : %d] *\n\n", start_pstn+1, start_pstn+counter);
  else 
    {
      fprintf(output_file,"<window>\n");
      fprintf(output_file,"<windowstart>%d</windowstart>\n",start_pstn+1);
      fprintf(output_file,"<windowend>%d</windowend>\n",start_pstn+counter);
    }
  buff[counter]='\0';
  nextLetterCheck='\0';
  nrRep = 0; 
  tooBigReps = NO;
  LastWindow=YES;
  PROCESS_ACGT(buff, counter);
  if (xmloutput==YES)
    fprintf(output_file,"</window>\n");
  free(buff-2);

  /* Check the situation */
  if (seq_pstn == to)  // TO has been reached
      /* skip the rest of the sequence */
      while ((workchar=RetrieveChar(situation,inputFile))!='>' && workchar!=EOF);
  else   // EOF or '>' has been reached
    if (toWasSpecified)
      printf("Warning: TO(=%d) is beyond the end of sequence\n",to);
  return workchar;
}
	
      
/* ************************************************************************* */

void ProcessSeqFromCommandLine(char * inputSeq)
/* The function processes the sequence read from the command line */
{ 
  int final_pstn;

  prevLetterCheck = '\0';
  start_pstn=from;
  
  if (step > 0) 
    while( (inputLength-start_pstn > 2*step) && (toWasSpecified==NO || start_pstn+2*step < to))
      {
	final_pstn=start_pstn+2*step;
	nextLetterCheck=inputSeq[final_pstn];
        inputSeq[start_pstn-1]=inputSeq[start_pstn-2]='\0';
	inputSeq[final_pstn]='\0';
	nrRep = 0; 
	tooBigReps = NO;
	LastWindow=NO;
	if (xmloutput == NO)
	  printf("\n* Processing window [%d : %d] *\n\n", start_pstn+1, final_pstn);
	else 
	  {
	    fprintf(output_file,"<window>\n");
	    fprintf(output_file,"<windowstart>%d</windowstart>\n",start_pstn+1);
	    fprintf(output_file,"<windowend>%d</windowend>\n",final_pstn);
	  }
	PROCESS_ACGT(inputSeq+start_pstn, 2*step);
	if (xmloutput==YES)
	  fprintf(output_file,"</window>\n");
	inputSeq[final_pstn]=nextLetterCheck;
	start_pstn+=step;
	prevLetterCheck = inputSeq[start_pstn-1];
      };

  /* Processing tail block */
  final_pstn = (toWasSpecified)? to : inputLength;
  nextLetterCheck = '\0';
  inputSeq[start_pstn-1]=inputSeq[start_pstn-2]='\0';
  inputSeq[final_pstn]='\0';
  nrRep = 0; 
  tooBigReps = NO;
  LastWindow=YES;
  if (xmloutput == NO)
    printf("\n* Processing window [%d : %d] *\n\n", start_pstn+1, final_pstn);
  else 
    {
      fprintf(output_file,"<window>\n");
      fprintf(output_file,"<windowstart>%d</windowstart>\n",start_pstn+1);
      fprintf(output_file,"<windowend>%d</windowend>\n",final_pstn);
    }
  PROCESS_ACGT(inputSeq+start_pstn, final_pstn-start_pstn);
      fprintf(output_file,"</window>\n");
}
  


/* ************************************************************************* */

char ReadSeqName(FILE * inputFile)
/* The function reads the sequence name (ACGT_FILE_FASTA case only).
   Only the first MAX_LINE_LENGTH characters are significant.
   The function prints sequence name and 
   returns '\n' if an end-of-line is reached, and 
   EOF if EOF is unexpectedly reached
*/
{
  char workchar;
  char * seqName;
  int counter=0;

  seqName = malloc ((1 + MAX_LINE_LENGTH)*sizeof(char));

  while ( (workchar = fgetc(inputFile))!='\n' && workchar!=EOF && counter<MAX_LINE_LENGTH )
    seqName[counter++] = workchar;
  seqName[counter] = '\0';
  if (xmloutput == NO)
    printf("Processing sequence '%s'\n", seqName+strspn(seqName," "));
  else
    fprintf(output_file,"<sequence-name>%s</sequence-name>\n", seqName+strspn(seqName," "));

  /* Check the situation */
  if (counter==MAX_LINE_LENGTH)
    while (workchar!='\n' && workchar!=EOF) workchar=fgetc(inputFile);
  if (workchar==EOF) 
    {
      printf("Warning: unexpected end-of-file after sequence name\n");
      free(seqName);
      return workchar;
    }
  /* workchar=='\n' */
  free(seqName);
  return workchar;
}

/* ************************************************************************* */

/* The main program */

int main (int argc,char * argv[])
{ 
  char *prog = argv[0] ;
  int input_type = INPUT_FILE ;   // type of input (seq or file)
  //  int batch = 0 ;
/*   int opt_space = 0;    // separators are different */
  int opt_fasta = NO;   // a plain sequence
  struct stat st ;      // file descriptor
  char workchar;        // working char
  FILE *input_file;            // input file

  verbose = 0;          // do not verbose

  /* let us start */

  init_limits ();       // initialize borderline values 

  logo();               // print the logo of the program 
  
  /* 1. Parsing the arguments of the command line */

  if (argc == 1)       // if no argument at all
    toofew (prog) ;
 
  argv++ ;             // look at the first argument

  if (!strcmp (*argv, "-h") || !strcmp (*argv, "-help")) /* help */ 
    help() ; 

/*   while (argc >= 3 || (argc == 2 && input_type == INPUT_STDIN) ) */
  while (argc >= 3 && input_type != INPUT_ARG)
    {
      /* 1.1 Which input ? */

/*       ARG_FLAG("-",    input_type = INPUT_STDIN) ; */
      ARG_FLAG("-s",   input_type = INPUT_ARG) ;   // must be followed by a sequence

      /* 1.2 Do I verbose ? */

      ARG_FLAG("-v",   verbose = 1) ;
      ARG_FLAG("-vv",  verbose = 2) ;
      ARG_FLAG("-vvv", verbose = 3) ;

      /* 1.3 Is it a batch ? 

	 ARG_FLAG("-b",   batch = 1) ; */

      /* 1.4 Allow non acgt characters ? */

/*       ARG_FLAG("-nacgt", allow_non_acgt = YES) ; */
/*       this option is suspended */
/*       it is currently only used in LET_TO_NUM which is in turn used only  */
/* 	 in the ASCII case, which makes this option meaningless */
      
      /* 1.5 Allow FASTA format ? */

      ARG_FLAG("-fasta", opt_fasta=YES) ;

     /* 1.6 Separators are the same character? */

/*       ARG_FLAG("-sspace", opt_space=YES) ; */

      /* 1.7 Limits on the sequence ? */

      ARG_INT("-from",   from);
      ARG_INT("-to",     to);

      ARG_INT("-w",   step);
      ARG_INT("-win",   step);
      ARG_INT("-window",   step);

      /* 1.8 Do I allow errors ? */

      ARG_INT("-err",    lim.err_number );
      ARG_INT("-res",    lim.err_number );
      ARG_INT("-resolution",    lim.err_number );

/*       ARG_INT("-k" , k); */

      /* 1.9 Limits on the repetitions ? */

/*       ARG_INT("-minP", lim.min_period ); */
      ARG_INT("-minp", lim.min_period );
/*       ARG_INT("-minPeriod", lim.min_period ); */
      ARG_INT("-minperiod", lim.min_period );

/*       ARG_INT("-maxP", lim.max_period); */
      ARG_INT("-maxp", lim.max_period);
/*       ARG_INT("-maxPeriod", lim.max_period ); */
      ARG_INT("-maxperiod", lim.max_period );

      ARG_INT("-minsize",   lim.min_size );
/*       ARG_INT("-minSize",   lim.min_size ); */
      ARG_INT("-maxsize",   lim.max_size );
/*       ARG_INT("-maxSize",   lim.max_size ); */

      ARG_FLOAT("-exp",  lim.min_exponent ); 

      /* 1.10 noprint option */

      ARG_FLAG("-noprint" , noprint=1);
      
      /* 1.11 allowsmall option */

      ARG_FLAG("-allowsmall" , allowsmall=YES);

     /* 1.12 Output in xml format? */

      if (!strcmp (*argv, "-xmloutput"))
	{
	  xmloutput=YES;
	  argc-- ; argv++ ;
	  output_file_name=*argv;
	  if ((output_file = fopen(*argv,"w")) == NULL)
	    {
	      fprintf(stderr,"Error: cannot open xml output file\n");
	      exit (2);
	    }
	  argc-- ; argv++ ;
	  continue;
	}

     /* 1.121 Conserver the unsorted mreps? */
#ifdef SHOWJUNK
      if (!strcmp (*argv, "-filter"))
	{ showjunk=YES;
	  argc-- ; argv++ ;
	  if ((junkfile = fopen(*argv,"w")) == NULL)
	    {
	      fprintf(stderr,"Error: cannot open file %s\n", *argv);
	      exit (22);
	    }
	  argc-- ; argv++ ;
	  continue;
	}
#endif

      /* 1.13 What's that? */
      fprintf(stderr,"Error: unrecognized option : %s.\n", *argv) ;
      wrongcall (prog) ;
    }

  if (argc == 1) 
    toofew (prog) ;

  /* 2. Validate parameters  */

  /* 2.0 print xml header */

  if (xmloutput==YES)
    {
      fprintf(output_file,"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
      fprintf(output_file,"\t<mreps>\n");
      now = time(NULL);
      fprintf(output_file,"\t\t<time>%s\t\t</time>\n",ctime(&now));
    }

  /* 2.1 Compute the length of the input */

  if (input_type == INPUT_FILE)   // input from file
    { 
      /* Open the file */

      if ((input_file = fopen(*argv,"rt")) == NULL) 
	{
	  fprintf(stderr,"Error: Cannot open file '%s'!\n",*argv);
	  PROCESS_ERROR(output_file,output_file_name,3)
	  exit (3);
	}

      /* Set variables */

      stat (*argv, &st);
      inputLength = (int) (st.st_size);   // take the length of the file
    }  
  else  // input from command line
    inputLength = strlen(*argv);

  /* 2.2 Validate MIN/MAX PERIOD and MIN/MAX SIZE*/

  if (lim.max_period == -1)       // MAXPERIOD was not specified
    lim.max_period=inputLength/2;
  else 
    {
      maxperiodWasSpecified = YES;
      if (lim.max_period>inputLength/2) lim.max_period=inputLength/2;
      if (lim.max_period < lim.min_period) 
	{
	  fprintf(stderr,"Error: Maximal period must be greater or equal to the minimal period! \n");
	  PROCESS_ERROR(output_file,output_file_name,4)
	  exit (4);
	}
    }

  if (lim.max_size == -1)         // MAXSIZE was not specified
    lim.max_size=inputLength;
  else
    {
      maxsizeWasSpecified = YES;
      if (lim.max_size < lim.min_size) 
	{
	  fprintf(stderr,"Error: Maximal size must be greater than or equal to the minimal size! \n");
	  PROCESS_ERROR(output_file,output_file_name,5)
	  exit (5);
	}
    }

  /* 2.3 Validate FROM and TO */
  
  if (from > inputLength)   
    { 

      fprintf(stderr,"Error: FROM(=%d) is too big!\n",from);
      PROCESS_ERROR(output_file,output_file_name,6)
      exit (6);
    }

  if (to > 0)             // TO was specified
    {
      toWasSpecified=YES;
      if (to < from) 
	{
	  fprintf(stderr,"Error: TO(=%d) must be greater or equal than FROM(=%d)\n",to,from); 
	  PROCESS_ERROR(output_file,output_file_name,7)
	  exit (7);
	}
      if (to > inputLength && input_type == INPUT_ARG)
	{
	  printf("Warning: TO(=%d) exceeds sequence length\n",to);
	  to = inputLength;
	}
    }
  else                  // TO was not specified
    to = inputLength;
  
  /* Decrement FROM to translate 'natural enumeration' to C indexing */

  from--;

  /* 2.4 Validate STEP */

  if (step > MAX_STEP) step = MAX_STEP;
  
  /* 2.5 Validate number of possible errors */

  VERB(1)  printf("Resolution parameter is %d\n", lim.err_number);

  /* 2.6 Validate ERR */
  if (lim.err_number>0 && lim.err_number >= inputLength/2-1)
    {
      fprintf(stderr,"Error: number of errors (=%d) is too big with respect to the sequence length\n",lim.err_number);
      PROCESS_ERROR(output_file,output_file_name,71)
      exit (71);
    }

  /* 2.7 Validate MINPERIOD and MAXPERIOD */
   if ( lim.min_period<1 )
      lim.min_period=1;

  VERB(1)  printf("Minimal possible period of repeats is %d\n", lim.min_period);  
  
/*   if (lim.max_size == -1) lim.max_size  = 2*step ; */

  /* 3. Identify the situation */

/*   if (opt_space) */
/*     { */
/*       fprintf(stderr,"Error: Option -sspace is meaningless in acgt version\n"); */
/*       PROCESS_ERROR(output_file,output_file_name,10) */
/*       exit (10); */
/*     } */
  if (opt_fasta)
    {
      if (input_type == INPUT_ARG)   // the input is an ACGT sequence in Fasta format
	{ 
	  fprintf(stderr,"Error: Option -fasta is not possible.\nWrite sequence in a file.\n");
	  PROCESS_ERROR(output_file,output_file_name,11)
	  exit (11);
	}

      else                   // the input is an ACGT file in Fasta format
	  situation = ACGT_FILE_FASTA;
    }

  if ((step==-1 && inputLength > MAX_INPUT_SHOW)|| (2*step > MAX_INPUT_SHOW))
    verbose_long = -1;

  /* 4. Print parameters into xml output file */

  if (xmloutput==YES)
    {
      fprintf(output_file,"\t\t<parameters>\n");
      /* 	  fprintf(output_file,"\t<type-of-input>"\n); */
      switch (input_type)
	{
	case INPUT_ARG: 
	  fprintf(output_file,"\t\t\t<type-of-input>command line</type-of-input>\n");
	  fprintf(output_file,"\t\t\t<input-sequence>%s</input-sequence>\n",*argv);
	  break;
	case INPUT_FILE:
	  if (situation == ACGT_FILE_FASTA)
	    fprintf(output_file,"\t\t\t<type-of-input>file in fasta format</type-of-input>\n");
	  else 
	    fprintf(output_file,"\t\t\t<type-of-input>file in plain format</type-of-input>\n");
	  break;
	default:
	  fprintf(output_file,"\t\t\t<type-of-input>unknown input</type-of-input>\n");
	}
      /* 	  fprintf(output_file,"</type-of-input>\n"\n); */
      fprintf(output_file,"\t\t\t<err>%d</err>\n",lim.err_number);
      fprintf(output_file,"\t\t\t<from>%d</from>\n",from+1);
      fprintf(output_file,"\t\t\t<to>%d</to>\n",(toWasSpecified)? to : -1);
      fprintf(output_file,"\t\t\t<win>%d</win>\n",step);
      fprintf(output_file,"\t\t\t<minsize>%d</minsize>\n",lim.min_size);
      fprintf(output_file,"\t\t\t<maxsize>%d</maxsize>\n",(maxperiodWasSpecified)? lim.max_size : -1);
      fprintf(output_file,"\t\t\t<minperiod>%d</minperiod>\n",lim.min_period);
      fprintf(output_file,"\t\t\t<maxperiod>%d</maxperiod>\n",(maxperiodWasSpecified)? lim.max_period : -1);
      fprintf(output_file,"\t\t\t<minexponent>%.2f</minexponent>\n",lim.min_exponent);
      fprintf(output_file,"\t\t</parameters>\n");
      fprintf(output_file,"\t\t<results>\n");
    }
	

  /* 5. Process the input */

  if (input_type == INPUT_ARG)
     ProcessSeqFromCommandLine(*argv);
  else
    if (situation == ACGT_FILE_FASTA)
      {
	if ((workchar = getc(input_file)) != '>')
	  {
	    fprintf(stderr,"Error: fasta file must start with '>'\n");
	    TERMINATE_XML_OUTPUT(output_file,output_file_name,15);
	    exit (15);
	  }
	while (workchar == '>')
	  {
	    workchar=ReadSeqName(input_file);
	    if (workchar != EOF)
	      {
		nrRep = 0; 
		tooBigReps=NO;
		if (xmloutput==YES) fprintf(output_file,"<repetitions>\n");
		workchar=ProcessSeqFromFile(input_file);
		if (xmloutput==YES) fprintf(output_file,"</repetitions>\n");
	      }
	  }
      }
    else
      workchar=ProcessSeqFromFile(input_file);	

#ifdef SHOWJUNK
  if (showjunk)
     fclose(junkfile);
#endif
  TERMINATE_XML_OUTPUT(output_file,output_file_name,0);
  return 0;
}

