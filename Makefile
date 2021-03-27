EXEC=Indenteur
CC=g++ -std=c++11

S=GRAPH

BITS=$(getconf LONG_BIT)
ifeq ($(BITS),64)
	FBITS=-m64
else
	FBITS=
endif

FLAGS=$(FBITS) -O3

all : Graph

Graph : main.o problems.o algos.o messages.o graph.o 
	$(CC) $(FLAGS) -o bin/Graph main.o problems.o algos.o messages.o graph.o 

main.o : $(S)/main.cpp $(S)/graph.h $(S)/messages.h $(S)/algos.h $(S)/problems.h
	$(CC) -c $(FLAGS) $(S)/main.cpp

problems.o : $(S)/problems.cpp $(S)/problems.h 
	$(CC) -c $(FLAGS) $(S)/problems.cpp

algos.o : $(S)/algos.cpp $(S)/algos.h 
	$(CC) -c $(FLAGS) $(S)/algos.cpp

messages.o : $(S)/messages.cpp $(S)/messages.h 
	$(CC) -c $(FLAGS) $(S)/messages.cpp

graph.o : $(S)/graph.cpp $(S)/graph.h 
	$(CC) -c $(FLAGS) $(S)/graph.cpp

#per cancellare i file oggetto fai make clean
clean :
	@rm *.o

rmexec :
	@rm $(EXEC)