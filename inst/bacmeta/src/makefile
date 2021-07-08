CC=       g++
STD=      -std=c++0x
CFLAGS=   -Wall $(STD) -O3   -Wno-unused-result   
OBJS=     Bac.o Gene.o Para.o Pop.o Simu.o Superpop.o

all:    simuli

simuli: $(OBJS)
	$(CC) $(OBJS) -o ../simu  

Bac.o:  Bac.cpp Bac.h Para.h 
	$(CC) $(CFLAGS) -c Bac.cpp

Gene.o: Gene.cpp Gene.h Para.h Superpop.h
	$(CC) $(CFLAGS) -c Gene.cpp

Para.o: Para.cpp Para.h
	$(CC) $(CFLAGS) -c Para.cpp

Pop.o:  Pop.cpp Pop.h Bac.h Para.h
	$(CC) $(CFLAGS) -c Pop.cpp

Simu.o: Simu.cpp
	$(CC) $(CFLAGS) -c Simu.cpp

Superpop.o:     Superpop.cpp Superpop.h Pop.h Para.h 
	$(CC) $(CFLAGS) -c Superpop.cpp

clean:
	rm -f *.o

clobber:        clean
	rm simu
