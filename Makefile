BOOSTDIR = '/usr/include'
LEDALIB = 'LEDA/LEDA'
LEDAINCL = 'LEDA/LEDA/incl'
LIBFLAGS = -lleda -lm

f = src/main.cpp

compile: $(f)
	g++ $(f) -o exec -O3 -std=c++0x -fno-strict-aliasing -I $(BOOSTDIR) -I $(LEDAINCL) -L $(LEDALIB) $(LIBFLAGS)

clean:
	rm -f exec

run:
	./exec