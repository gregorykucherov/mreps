#include "defs.h"
#define EXPERIMENT
extern int start_pstn, nrRep, noprint;
extern char * seq_original;
extern int xmloutput;
extern FILE *output_file;


void print_rep(int rinitpos, int rendpos, int rlength, int rperiod, int rnumerr)
     /* printing an exact or approximate rep */
{    
  if ((nrRep++)==0 && xmloutput==NO)
    {  
      printf("   from   ->       to  : \t size\t <per.>\t [exp.]\t errors: \trepetition\n") ;
      printf(" ---------------------------------------------------------------------------------------------\n");
    }

  if (xmloutput==YES)
    fprintf (output_file,
	     "\t<repeat>\n\t\t<start>%d</start>\n\t\t<end>%d</end>\n\t\t<length>%d</length>\n\t\t<period>%d</period>\n\t\t<exponent>%.2f</exponent>\n\t\t<errors>%d</errors>\n",
	     start_pstn+rinitpos, start_pstn+rendpos, rlength, rperiod,((float) rlength) / rperiod, rnumerr) ; 
  else
    printf("%8d  ->  %8d :   \t %d \t <%d> \t [%.2f] \t %d \t",
	   start_pstn+rinitpos, start_pstn+rendpos, rlength, rperiod, ((float) rlength) / rperiod, rnumerr
	   ) ;

  if (!noprint)
    {
      int per_start,pos_in_per;
/*       int re = (rendpos<rinitpos + MAXDISPLAY)? rendpos : rinitpos+(rperiod-1); */
      if (xmloutput==YES)
	fprintf(output_file,"\t\t<sequence>\n"); 

      //output the repetitions with a space after each period
/*       for (position = rinitpos; position <= re; position++) */
/*       for (position = rinitpos; position <= rendpos; position++) */
/* 	{ */
/* 	  if (xmloutput==YES) */
/* 	    fputc(seq_original[position-1],output_file) */
/* 	  else */
/* 	    printf ("%c",seq_original[position-1]); */
/* 	  if (((position+1-rinitpos)%rperiod)==0) */
/* 	    printf(" "); */
/* 	}; */

      for (per_start = rinitpos; per_start <= rendpos; per_start+=rperiod)
	{
	  if (xmloutput==YES)
	    fprintf(output_file,"\t\t\t<unit>"); 
	  for (pos_in_per=0; pos_in_per<rperiod && per_start+pos_in_per<=rendpos; pos_in_per++)
	    {
	      if (xmloutput==YES)
		fputc(seq_original[per_start+pos_in_per-1],output_file);
	      else 
		printf ("%c",seq_original[per_start+pos_in_per-1]);
	    }
	  if (xmloutput==YES)
	    fprintf(output_file,"</unit>\n");
	  else 
	    printf(" ");
	}
      if (xmloutput==YES)
	fprintf(output_file,"\t\t</sequence>\n");

    }
  if (xmloutput==YES)
    fprintf(output_file,"\t</repeat>\n"); 
  else
    printf ("\n");
}

#ifdef EXPERIMENT
void print_score(int rinitpos, int rendpos, int rlength, int rperiod, float rscore)
     /* printing an exact or approximate rep */
{    
  if ((nrRep++)==0 && xmloutput==NO)
    {  
      printf("   from   ->       to  : \t size\t <per.>\t [exp.]\t\t err-rate \tsequence\n") ;
      printf(" ---------------------------------------------------------------------------------------------\n");
    }

  if (xmloutput==YES)
    fprintf (output_file,
	     "\t<repeat>\n\t\t<start>%d</start>\n\t\t<end>%d</end>\n\t\t<length>%d</length>\n\t\t<period>%d</period>\n\t\t<exponent>%.2f</exponent>\n\t\t<score>%.3f</score>\n",
	     start_pstn+rinitpos, start_pstn+rendpos, rlength, rperiod,((float) rlength) / rperiod, rscore) ; 
  else
    printf("%8d  ->  %8d :   \t %d \t <%d> \t [%.2f] \t %.3f \t\t",
	   start_pstn+rinitpos, start_pstn+rendpos, rlength, rperiod, ((float) rlength) / rperiod, rscore
	   ) ;

  if (!noprint)
    {
      int per_start,pos_in_per;
/*       int re = (rendpos<rinitpos + MAXDISPLAY)? rendpos : rinitpos+(rperiod-1); */
      if (xmloutput==YES)
	fprintf(output_file,"\t\t<sequence>\n"); 

      for (per_start = rinitpos; per_start <= rendpos; per_start+=rperiod)
	{
	  if (xmloutput==YES)
	    fprintf(output_file,"\t\t\t<unit>"); 
	  for (pos_in_per=0; pos_in_per<rperiod && per_start+pos_in_per<=rendpos; pos_in_per++)
	    {
	      if (xmloutput==YES)
		fputc(seq_original[per_start+pos_in_per-1],output_file);
	      else 
		printf ("%c",seq_original[per_start+pos_in_per-1]);
	    }
	  if (xmloutput==YES)
	    fprintf(output_file,"</unit>\n");
	  else 
	    printf(" ");
	}
      if (xmloutput==YES)
	fprintf(output_file,"\t\t</sequence>\n");

    }
  if (xmloutput==YES)
    fprintf(output_file,"\t</repeat>\n"); 
  else
    printf ("\n");
}
#endif

/* void printLocal(int start_pstn , int rinitpos, int rendpos, int rlength, int rperiod) { */

/*    printf ("%8d  ->  %8d :   \t %d \t <%d> \t [%.2f] \t", */
/* 	   start_pstn+rinitpos, start_pstn+rendpos-1, rlength, rperiod, */
/* 	   ((float) rlength) / rperiod */
/* 	  ) ;  */
/* } */

/* void printLocalMismatch(listreps current) { */

/*   printf( "%8d  ->  %8d :   \t %d \t <%d> \t [%.2f] \t" , (current->rep).initpos+1 , (current->rep).endpos-1 , (current->rep).endpos-(current->rep).initpos-1, (current->rep).period, ((float)((current->rep).endpos-(current->rep).initpos-1)/(current->rep).period));   */

/* } */

/* void printWeb(int start_pstn , int rinitpos, int rendpos, int rlength, int rperiod) { */
/*   fprintf (stdout, "%d %d %d %d %.2f ", */
/* 	   start_pstn+rinitpos, start_pstn+rendpos-1, rlength, rperiod, */
/* 	   ((float) rlength) / rperiod); */
/* } */

/* void printWebMismatch(listreps current) { */

/*   printf( "%d %d %d %d %.2f " , (current->rep).initpos+1 , (current->rep).endpos-1 , (current->rep).endpos-(current->rep).initpos-1, (current->rep).period, ((float)((current->rep).endpos-(current->rep).initpos-1)/(current->rep).period)); */
/* } */
