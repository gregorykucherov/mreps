/* defs.h : Macros */
/*******************/

/* All global definitions are here */

#include <stdio.h>
#include <stdlib.h>

#ifndef __REPS_DEFS
#define __REPS_DEFS

#define ACGT_FILE_FASTA   1
/* #define ACGT_SEQ_PLAIN    2 */
/* #define ACGT_FILE_PLAIN   3 */
#define ASCII_FILE_SSPACE 5 
/* #define ASCII_SEQ         6 */
/* #define ASCII_FILE        7 */

/* Constants */
/* #define MAXDISPLAY 2000 */

/* Factorization options */

#define WITH_OVERLAP 0
#define WITHOUT_OVERLAP 1

/* Used for parameter FACTYPE in roman/factorizeforGDR.c */
/* Lempel-Ziv factorization for k-mismatch */
#define WITH_EXTRA_SYMBOL 0

/* s-factorization for exact case */
#define WITHOUT_EXTRA_SYMBOL 1


/* Globals */

extern int verbose ;         /* This externs are set in main.c */
extern int verbose_long ;

/* extern int allow_non_acgt ;  // This one is in dawg.c  */

#define VERB(i)   if ( verbose >= i ) 
#define VERBL(i)  if ( verbose_long >= i ) 
#define YES  1
#define NO   0

#define VIRTUAL_STEP ((LastWindow) ? 2*step : step)

/* Types */

struct s_limits
{
  int err_number ;
  int min_period ;
  int max_period ;        
  float min_exponent ;
  int max_exponent ;       /* Not implemented */
  int min_size ;
  int max_size ;           
} ;

typedef struct s_limits limits ;

struct s_sfactorization
{
  int nbfactors ;
  int *factor ;
  int *previous ;
} ;

typedef struct s_sfactorization *sfactorization ;

struct s_repetition
{
  char type ;
  int period ;
  int initpos ;
  int endpos ;
} ;

typedef struct s_repetition repetition ;

struct s_listreps
{
  repetition rep ;
  struct s_listreps *next ;
  
} ;

typedef struct s_listreps *listreps ;

/* for accounting repetitions crossing the bounds
 * between steps
 */
struct s_transreps
{
  int period ;
  struct s_transreps *next ;
  
} ;

typedef struct s_transreps *transreps ;


/* Now the functions */
/* ----------------- */

/* let2num */

#define LET_TO_NUM(c) ((unsigned char) c)

#endif /* __REPS_DEFS */

#define TERMINATE_XML_OUTPUT(output_file,output_file_name,error_code)   \
{                                                            \
  if (xmloutput == YES)                                      \
    {                                                        \
      /*fprintf(output_file,"</window>\n");                    */\
      fprintf(output_file,"</results>\n");                   \
      fprintf(output_file,"<errorcode>%d</errorcode>\n",error_code);  \
      fprintf(output_file,"</mreps>\n");                      \
      printf("results in xml format are written in file %s\n",output_file_name);  \
    }                                                        \
}

#define RECOVER_XML_OUTPUT(output_file,output_file_name,error_code)   \
{                                                            \
  if (xmloutput == YES)                                      \
    {                                                        \
      fprintf(output_file,"</window>\n");                   \
      fprintf(output_file,"</repetitions>\n");                   \
      fprintf(output_file,"</results>\n");                   \
      fprintf(output_file,"<errorcode>%d</errorcode>\n",error_code);  \
      fprintf(output_file,"</mreps>\n");                      \
      printf("results in xml format are written in file %s\n",output_file_name);  \
    }                                                        \
}

#define PROCESS_ERROR(output_file,output_file_name,error_code)   \
{                                                            \
  if (xmloutput == YES)                                      \
    {                                                        \
      fprintf(output_file,"<errorcode>%d</errorcode>\n",error_code);  \
      fprintf(output_file,"</mreps>\n");                      \
      printf("error %d is registered in xml file %s\n",error_code,output_file_name);  \
    }                                                        \
}

#define LIMIT_N_PROPORTION 0.05
