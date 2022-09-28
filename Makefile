BOOSTDIR = '/usr/include/BOOSTDIR'
LEDALIB = 'LEDA/LEDA'
LEDAINCL = 'LEDA/LEDA/incl'
LIBFLAGS = -lleda -lm

f = src/main3.cpp

compile: $(f)
	g++ $(f) -o exec -O3 -std=c++0x -fno-strict-aliasing -lboost_system -fopenmp -I $(LEDAINCL) -L $(LEDALIB) $(LIBFLAGS)

clean:
	rm -f exec

run:
	./exec