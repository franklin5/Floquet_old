CC=mpicxx
eigenbin = /home/ld7/bin	
CFLAGS=-c -Wall -I${eigenbin}
LDFLAGS=-limf -lm
SRC=main.cpp bdg.cpp prx.cpp dist.cpp 
OBJ=$(SRC:.cpp=.o)
EXE=floquet

all: $(SRC) $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)
.cpp.o: 
	$(CC) $(CFLAGS) $< -o $@

touch: 
	touch *.cpp *.h
clean:
	rm *.o floquet *~ *#

