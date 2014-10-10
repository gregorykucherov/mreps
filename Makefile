# Makefile


#CC	=	gcc  -Wall -g -lm  
#CC	=	gcc  -g -lm  
CC	=	gcc  -O3
#LDFLAGS = 	 -L/usr/openwin/lib""	
#LDFLAGS =       -L/usr/X11R6/lib  
#LIBS	=       -lm -lX11

##########

all: mreps

##########


mreps:	defs.h mreps_acgt.o mainrepets.o kmp.o main_acgt.o FndReps.o factorizeforGDR.o mainsearchforGDR.o searchforHeadGDR.o finalstageforGDR.o printOutput_acgt.o
	$(CC) -o mreps mreps_acgt.o mainrepets.o kmp.o main_acgt.o  FndReps.o factorizeforGDR.o mainsearchforGDR.o searchforHeadGDR.o finalstageforGDR.o printOutput_acgt.o -lm


########## 

main_acgt.o:	defs.h main.c
	$(CC) -c -o $@ main.c

mreps_acgt.o: defs.h mreps.c
	$(CC) -c -o $@ mreps.c

kmp.o: defs.h  kmp.c
	$(CC) -c -o $@ kmp.c

printOutput_acgt.o: defs.h printOutput.c
	$(CC) -c -o $@ printOutput.c

FndReps.o: defs.h FndReps.c 
	$(CC) -c -o $@ FndReps.c

factorizeforGDR.o: defs.h  factorizeforGDR.c
	$(CC) -c -o $@ factorizeforGDR.c

mainsearchforGDR.o: defs.h mainsearchforGDR.c
	$(CC) -c -o $@ mainsearchforGDR.c

searchforHeadGDR.o: defs.h searchforHeadGDR.c
	$(CC) -c -o $@ searchforHeadGDR.c

finalstageforGDR.o: defs.h finalstageforGDR.c
	$(CC) -c -o $@ finalstageforGDR.c


mainrepets.o: defs.h mainrepets.c
	$(CC) -c -o $@ mainrepets.c

##########
