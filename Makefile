CC=mpicxx
eigenbin = /home/ld7/bin	
CFLAGS=-c -Wall -I${eigenbin}
LFLAGS=-limf -lm
all: floquet
floquet: main.o bdg.o dist.o prx.o
	$(CC) *.o -o floquet $(LFLAGS)

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

bdg.o: bdg.cpp
	$(CC) $(CFLAGS) bdg.cpp

prx.o: prx.cpp
	$(CC) $(CFLAGS) prx.cpp

dist.o: dist.cpp
	$(CC) $(CFLAGS) dist.cpp



touch: 
	touch *.cpp *.h
clean:
	rm *.o floquet *~ *#

