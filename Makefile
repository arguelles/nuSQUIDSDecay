# FLAGS

CFLAGS= -O3 -fPIC -std=c++11 -g
CFLAGS+= -I./inc `pkg-config --cflags nusquids`
LDFLAGS+= `pkg-config --libs nusquids`

all: mains/exCross.o mains/run_decay

mains/exCross.o : inc/exCross.h mains/exCross.cpp
	@$(CXX) $(CFLAGS) -c mains/exCross.cpp -o $@

mains/run_decay : mains/run_decay.cpp mains/exCross.o
	@echo Compiling run_decay
	@$(CXX) $(CFLAGS) mains/run_decay.cpp mains/exCross.o $(LDFLAGS) -o $@

.PHONY: clean
clean:
	rm -rf ./mains/run_decay ./mains/exCross.o



