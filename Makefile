# FLAGS
CFLAGS= -O3 -fPIC -std=c++11
CFLAGS+= -I./inc `pkg-config --cflags squids nusquids hdf5`
LDFLAGS+= `pkg-config --libs squids nusquids hdf5` -lpthread

all: mains/exCross.o  mains/partial_rate_example mains/couplings_example

mains/exCross.o : inc/exCross.h mains/exCross.cpp
	@ $(CXX) $(CFLAGS) -c mains/exCross.cpp -o $@

mains/partial_rate_example : mains/partial_rate_example.cpp mains/exCross.o inc/nusquids_decay.h
	@echo Compiling partial_rate_example
	@ $(CXX) $(CFLAGS) mains/partial_rate_example.cpp mains/exCross.o  -o $@ $(LDFLAGS)

mains/couplings_example : mains/couplings_example.cpp mains/exCross.o inc/nusquids_decay.h
	@echo Compiling couplings_example
	@ $(CXX) $(CFLAGS) mains/couplings_example.cpp mains/exCross.o  -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf ./mains/partial_rate_example ./mains/couplings_example ./mains/exCross.o

