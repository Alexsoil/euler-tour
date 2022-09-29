BOOSTDIR = '/usr/include/BOOSTDIR'
LEDALIB = 'LEDA/LEDA'
LEDAINCL = 'LEDA/LEDA/incl'
LIBFLAGS = -lleda -lm

f = src/main4.cpp

compile: $(f)
	g++ $(f) -o exec -O3 -std=c++0x -fno-strict-aliasing -lboost_system -fopenmp -I $(LEDAINCL) -L $(LEDALIB) $(LIBFLAGS)

clean:
	rm -f exec

run:
	./exec

leda:
	g++ src/leda_tree.cpp -o exec2 -O3 -std=c++0x -fno-strict-aliasing -lboost_system -I $(LEDAINCL) -L $(LEDALIB) $(LIBFLAGS)	

runleda:
	./exec2