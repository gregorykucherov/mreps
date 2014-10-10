#include <math.h>
#include "../defs.h"

/* #include "factorizeforGDR.h" */

#define SECARRSIZ 10000

extern int lenWrd, maxNumErr;

void InitPrTree(void);
void MovePrFacNode(void);
void MakeSufTree(void);
void EnterSecTree(void);
void TreatInitSecNode(void);
void MakeInitSecNode(void);
void GoInitFacNode(void);
void MoveSecFacNode(void);
void TreatOldSecNode(void);
void MakeNewSecNode(void);
void GoOldSecNode(void);
void ComputeFinFactor(void);
void ComputeGlobalFactors(void);   // called only in FndReps.c
void MakePowof4(void);
void guessDepth(int length);   // called in mreps.c FndReps.c
int ComputeSizeFact(void);


int pow4[16], sumpow4[16], prPow4, prSumpow4; 
char *seq_factor_copy;
int maxNumErr2, actPerBnd;
int *prTree, *lvList, prLeaf;
int (**secTree)[5], numSecNod, numSecArr;
int blockSize, finBlockSize, *maxSumFact, sizeFact;
int *copyEnds, *copyBgns, *factorEnds;
int *factorBgns;    // factorBgns[i] is the last (head) symbol 
                    // (in natural enumeration) of i-th factor
int numFact, numGlobFact, numFullBlock;
int curPrNod, curFacNod, curPstn, curSmb;
int facSecDep, facSecLeaf, *facSecNod;  
int newList[5], newSecLeaf[5], nxtSecLeaf, curSecLeaf;
int prevCopyPstn, minSecDep;
int PRIMEDEPTH, FACTYPE, prDepth1;


void guessDepth(int length) {

  double  logSize ;
  int depth;
/*   length /= 1000; */
  logSize = log10((double)(length /1000.0));

  //depth is approximately linear with logSize. depth = 2x + 2;
 
  depth = (int)(logSize*2+2);
  if (depth > 10) depth = 10;
  else if (depth < 4) depth = 4;
 
  PRIMEDEPTH = depth;
  prDepth1=PRIMEDEPTH+1;
}


sfactorization  FactorizeforGDR(void) 
{
  sfactorization fStruct;
  int indx;
  fStruct= (sfactorization) malloc (sizeof (struct s_sfactorization)) ;
  MakePowof4();
  
  sizeFact=ComputeSizeFact();
  if (FACTYPE==WITH_EXTRA_SYMBOL)
    {  
      if ((fStruct->previous=copyEnds=
	      (int *)malloc((sizeFact+1)*sizeof(int))) == NULL )
	{  
	  printf("Not enough memory for factorization\n");
	  exit(29); 
	}
      if ((fStruct->factor=factorEnds=
	    (int *)malloc( (sizeFact+1)*sizeof(int))) == NULL )
	{  
	  printf("Not enough memory for factorization\n");
	  exit(30); 
	}
      *(factorEnds++)=(*(copyEnds++))=0;
      prevCopyPstn=(*(factorEnds++))=(*(copyEnds++))=1;
    }
  else
    {  
      if ((fStruct->previous=copyBgns=
	      (int *)malloc( (sizeFact+1)*sizeof(int))) == NULL )
	{  
	  printf("Not enough memory for factorization\n");
	  exit(29); 
	}
    if ((fStruct->factor=factorBgns=
	    (int *)malloc( (sizeFact+1)*sizeof(int))) == NULL )
      {  
	printf("Not enough memory for factorization\n");
	exit(30); 
      }
    *(factorBgns++)=(*(copyBgns++))=0;
    *factorBgns=1;
    }
  numFact=2;
  if ( (prTree=(int *)calloc(sumpow4[PRIMEDEPTH+2], sizeof(int)))==NULL )
    {  
      printf("Not enough memory for factorization\n");
      exit(9); 
    }
  secTree=(int (**)[5])calloc(lenWrd/SECARRSIZ+1, sizeof(int (*)[5]));
  if ((secTree[numSecArr=0]=(int (*)[5])calloc(SECARRSIZ, sizeof(int [5])))==NULL )
    {  
      printf("Can't allocate memory for first secondary array\n");
      exit(11);
    }
  if ( (lvList=(int *)calloc(lenWrd+1, sizeof(int)))==NULL )
    {  
      printf("Can't allocate memory for factorization\n");
      exit(10); 
    }
  facSecDep=0;
  curFacNod=0;
  numSecNod=0;
  InitPrTree();
  while(curPstn<=lenWrd)
    {  
      MakeSufTree();
      curPstn++;
    }
  ComputeFinFactor();
  free(lvList);
  for(indx=0;indx<=numSecArr;indx++)
    free(secTree[indx]);
  free(secTree);
  free(prTree);
  fStruct->nbfactors=numFact;
  if (FACTYPE==WITH_EXTRA_SYMBOL) 
    {  
    copyEnds=fStruct->previous;
    factorEnds=fStruct->factor;
    }
  else
    {  
      copyBgns=fStruct->previous;
      factorBgns=fStruct->factor;
    }
  return(fStruct);
}


void InitPrTree(void) {

  int sufNode, sufLen, maxPstn;
  prTree[0]=-1;
  prTree[curPrNod=(int)seq_factor_copy[1]]=1;
  maxPstn=(PRIMEDEPTH<lenWrd) ? PRIMEDEPTH : lenWrd;
  
  for(curPstn=2;curPstn<=maxPstn;curPstn++) {
    
    curSmb=(int)seq_factor_copy[curPstn];
    MovePrFacNode();
    curPrNod=4*curPrNod+curSmb;
    prTree[curPrNod]=curPstn;
    sufLen=curPstn-1;
    sufNode=(curPrNod-sumpow4[sufLen])%pow4[sufLen]+sumpow4[sufLen];
    while(prTree[sufNode]==0) {
      prTree[sufNode]=curPstn;
      sufLen--;
      sufNode=(sufNode-sumpow4[sufLen])%pow4[sufLen]+sumpow4[sufLen];
    }
  }
}


void MovePrFacNode(void) {
  
  int curCopyPstn, curFacBgn;
  curFacNod=4*curFacNod+curSmb;
  if ((curCopyPstn=prTree[curFacNod]))
    prevCopyPstn=curCopyPstn;
  else {
    if (FACTYPE==WITH_EXTRA_SYMBOL) {
      *(copyEnds++)=prevCopyPstn+1;
      *(factorEnds++)=prevCopyPstn=curPstn;
      numFact++;
      curFacNod=0;
    }
    else {
      if ( (curFacBgn=curPstn-1)>(*factorBgns) ) {
	*(copyBgns++)=prevCopyPstn+(*factorBgns)-curFacBgn;
	*(++factorBgns)=curFacBgn;
	numFact++;
	if ( (prevCopyPstn=prTree[curFacNod=curSmb])==0 )
	  {  *(copyBgns++)=curFacBgn;
	  *(++factorBgns)=curPstn;
	  numFact++;
	  curFacNod=0;
	  }
      }
      else {
	*(copyBgns++)=curFacBgn;
	*(++factorBgns)=curPstn;
	numFact++;
	curFacNod=0;
      }
    }
  }
}


void MakeSufTree(void) {
  int sufNode, sufLen;
  curSmb=(int)seq_factor_copy[curPstn];
  prLeaf=4*curPrNod+curSmb;
  if (curFacNod<prSumpow4) {
    if (curFacNod<0)
      MoveSecFacNode();
    else
      MovePrFacNode();
    lvList[curPstn]=prTree[prLeaf];
    prTree[prLeaf]=curPstn;   
  }
  else
    EnterSecTree();  
  
  sufNode=curPrNod=(prLeaf-prSumpow4)%prPow4+prSumpow4;
  sufLen=PRIMEDEPTH;
     while(prTree[sufNode]==0) {
       prTree[sufNode]=curPstn;
       sufLen--;
       sufNode=
	 (sufNode-sumpow4[sufLen])%pow4[sufLen]+sumpow4[sufLen];
     }
}


void EnterSecTree(void) {
  int indx, actSmb;
  if ( (curFacNod=prTree[prLeaf])>0 ) {
    for(indx=1;indx<=4;indx++)
      newList[indx]=0;
    do {
      curSecLeaf=curFacNod;
      actSmb=seq_factor_copy[curSecLeaf+1];
      if (newList[actSmb])
	newSecLeaf[actSmb]=lvList[newSecLeaf[actSmb]]=
	  curSecLeaf;
      else
	newList[actSmb]=newSecLeaf[actSmb]=curSecLeaf;
    }  while( (curFacNod=lvList[curFacNod])>0 );
    if (curFacNod<0)
      TreatInitSecNode();
    else
      MakeInitSecNode();
  }
  else {
    if (curFacNod<0)
      GoInitFacNode();
    else { 
      if (FACTYPE==WITH_EXTRA_SYMBOL) {
	*(copyEnds++)=prevCopyPstn+1;
	*(factorEnds++)=prevCopyPstn=
	  prTree[prLeaf]=curPstn;
	numFact++;
	curFacNod=0;
      }
      else {
	*(copyBgns++)=prevCopyPstn-PRIMEDEPTH;
	*(++factorBgns)=curPstn-1;
	prTree[prLeaf]=curPstn;
	numFact++;
	if ( (prevCopyPstn=prTree[curFacNod=curSmb])==0 ) {
	  *(copyBgns++)=(*factorBgns);
	  *(++factorBgns)=curPstn;
	  numFact++;
	  curFacNod=0;
	}
      }
    }
  }
}


void TreatInitSecNode(void) { 
  
  int indx;
  facSecNod=(int *)
    (secTree[(-curFacNod)/SECARRSIZ]+(-curFacNod)%SECARRSIZ);
  for(indx=1;indx<=4;indx++) { 
    if (newList[indx]) {
      lvList[newSecLeaf[indx]]=facSecNod[indx];
      facSecNod[indx]=newList[indx];
    }
  }
  prTree[prLeaf]=curFacNod;
  facSecLeaf=curPstn;
  facSecDep=1;
}


void MakeInitSecNode(void)
{  int nodIndx, indx;
   if ( (nodIndx=(++numSecNod)%SECARRSIZ)==0 )
   {  if ( ( secTree[++numSecArr]=(int (*)[5])calloc(SECARRSIZ, sizeof(int [5])) )
             ==NULL )
      {  printf("Can't allocate memory for %d-th secondary array\n", numSecArr+1);
         exit(11);
      }
   }
   facSecNod=(int *)(secTree[numSecArr]+nodIndx);
   *facSecNod=curSecLeaf;
   for(indx=1;indx<=4;indx++)
   {  if (newList[indx])
      {  facSecNod[indx]=newList[indx];
         lvList[newSecLeaf[indx]]=0;
      }
   }
   curFacNod=prTree[prLeaf]=-numSecNod;
   facSecLeaf=curPstn;
   facSecDep=1;
}


void GoInitFacNode(void)
{  facSecNod=(int *)
   (secTree[(-curFacNod)/SECARRSIZ]+(-curFacNod)%SECARRSIZ);
   facSecLeaf=curPstn;
   facSecDep=1;
}


void MoveSecFacNode(void)
{
  int indx, actSmb;
  if ( (nxtSecLeaf=facSecNod[curSmb])>0 )
    {  facSecDep++;
    for(indx=1;indx<=4;indx++)
      newList[indx]=0;
    do
      { 
	curSecLeaf=nxtSecLeaf;
	actSmb=seq_factor_copy[curSecLeaf+facSecDep];
	if (newList[actSmb])
	  newSecLeaf[actSmb]=lvList[newSecLeaf[actSmb]]=
	    curSecLeaf;
	else
	  newList[actSmb]=newSecLeaf[actSmb]=curSecLeaf;
      }  while( (nxtSecLeaf=lvList[nxtSecLeaf])>0 );
    if (nxtSecLeaf<0)
      TreatOldSecNode();
    else
      MakeNewSecNode();
    }
  else
    {
      if (nxtSecLeaf<0)
	GoOldSecNode();
      else
	{ 
	  facSecNod[curSmb]=facSecLeaf;
	  if (FACTYPE==WITH_EXTRA_SYMBOL)
         {  *(copyEnds++)=(*facSecNod)+facSecDep;
	 *(factorEnds++)=prevCopyPstn=curPstn;
	 numFact++;
	 curFacNod=0;
         }
	  else
         {
	   *(copyBgns++)=(*facSecNod)-prDepth1;
	   *(++factorBgns)=curPstn-1;
	   numFact++;
	   if ( (prevCopyPstn=prTree[curFacNod=curSmb])==0 )
	     {
	       *(copyBgns++)=(*factorBgns);
               *(++factorBgns)=curPstn;
               numFact++;
               curFacNod=0;
	     }
         }
	}
    }
}


void TreatOldSecNode(void)
{  int indx, *nxtSecNod;
   nxtSecNod=(int *)
   (secTree[(-nxtSecLeaf)/SECARRSIZ]+(-nxtSecLeaf)%SECARRSIZ);
   for(indx=1;indx<=4;indx++)
   {  if (newList[indx])
      {  lvList[newSecLeaf[indx]]=nxtSecNod[indx];
         nxtSecNod[indx]=newList[indx];
      }
   }
   facSecNod[curSmb]=nxtSecLeaf;
   facSecNod=nxtSecNod; 
}


void MakeNewSecNode(void)
{  int nodIndx, indx, *nxtSecNod;
   if ( (nodIndx=(++numSecNod)%SECARRSIZ)==0 )
   {  if ( ( secTree[++numSecArr]=(int (*)[5])calloc(SECARRSIZ, sizeof(int [5])) )
             ==NULL )
      {  printf("Can't allocate memory for %d-th secondary array\n", numSecArr+1);
         exit(12);
      }
   }
   nxtSecNod=(int *)(secTree[numSecArr]+nodIndx);
   *nxtSecNod=curSecLeaf;
   for(indx=1;indx<=4;indx++)
   {  if (newList[indx])
      {  nxtSecNod[indx]=newList[indx];
         lvList[newSecLeaf[indx]]=0;
      }
   }
   facSecNod[curSmb]=-numSecNod;
   facSecNod=nxtSecNod; 
}


void GoOldSecNode(void)
{  int *nxtSecNod;
   nxtSecNod=(int *)
   (secTree[(-nxtSecLeaf)/SECARRSIZ]+(-nxtSecLeaf)%SECARRSIZ);
   facSecNod=nxtSecNod; 
   facSecDep++;
}


void ComputeFinFactor(void)
{  if (FACTYPE==WITH_EXTRA_SYMBOL)  
   {  if (curFacNod<0)
         *copyEnds=(*facSecNod)+facSecDep;
      else
         *copyEnds=prevCopyPstn+1;
      *factorEnds=curPstn;
   }
   else
   {  if (curFacNod)
      {  if (curFacNod<0)
            *copyBgns=(*facSecNod)-prDepth1;
         else
            *copyBgns=prevCopyPstn+(*factorBgns)-lenWrd;
         *(++factorBgns)=lenWrd;
      }
      else
         numFact--;
   }
}


void ComputeGlobalFactors(void)
{    
  int blockSizeBnd, fctIndx, subFctIndx, blockIndx, sumFact;
  actPerBnd=0;
  blockSize=2;

  if ( (blockSizeBnd=maxNumErr+1)>2 )
    {
      while( 2*blockSize<=blockSizeBnd )
	blockSize*=2;
    }
VERB(1)  printf("Global factor has %d Lempel-Ziv factors\n", blockSize);
  numFullBlock=(numFact-1)/blockSize;
  numGlobFact=numFullBlock+1;
VERB(1)  printf("Input sequence has %d global factors\n\n", numGlobFact);
  
  finBlockSize=numFact-numFullBlock*blockSize;
  
  maxSumFact=(int *)calloc( numGlobFact+1, sizeof(int) );
  fctIndx=(maxNumErr2>blockSize) ? maxNumErr2 : blockSize;
  subFctIndx=fctIndx-maxNumErr2;
  
  while(fctIndx<=numFact)
    {
      sumFact=factorEnds[fctIndx]-factorEnds[subFctIndx];
      if (sumFact>actPerBnd)
	actPerBnd=sumFact;
      blockIndx=fctIndx/blockSize;
      while(blockIndx*blockSize>subFctIndx)
	{
	  if (sumFact>maxSumFact[blockIndx])
	    maxSumFact[blockIndx]=sumFact;
	  blockIndx--;
	}
      fctIndx++;
      subFctIndx++;
    }
  maxSumFact[numGlobFact]=
    factorEnds[numFact]-factorEnds[numFact-maxNumErr2];
}


void MakePowof4(void)
{ 
  int indx;
  sumpow4[0]=0;
  pow4[0]=1;
  for(indx=1;indx<16;indx++)
    { 
      pow4[indx]=4*pow4[indx-1];
      sumpow4[indx]=sumpow4[indx-1]+pow4[indx-1];
    }
  prPow4=pow4[PRIMEDEPTH]; 
  prSumpow4=sumpow4[PRIMEDEPTH];
}


int ComputeSizeFact(void)
{    
  int depth, sumLen=0;
  if (FACTYPE==WITH_EXTRA_SYMBOL)
    {  
      for(depth=1;depth<16;depth++)
	{  
	  if ( (sumLen+depth*pow4[depth])>=lenWrd )
	    return(sumpow4[depth]+(lenWrd-sumLen)/depth);
	  else
	    sumLen+=depth*pow4[depth];
	}
    }
  else
    {  
      for(depth=0;depth<15;depth++)
	{ 
	  if ( (sumLen+depth*pow4[depth+1])>=lenWrd )
	    return(sumpow4[depth+1]+(lenWrd-sumLen)/depth);
	  else
	    sumLen+=depth*pow4[depth+1];
	}
    }
  /* this situation cannot happen unless the sequence length is about 20 10^9 */
  fprintf(stderr,"the input sequence is too long\n");
  exit(100);
} 
