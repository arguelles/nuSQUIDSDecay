# FLAGS
CFLAGS= -O3 -fPIC -std=c++11
CFLAGS+= -I./include `pkg-config --cflags squids nusquids hdf5`
LDFLAGS+= `pkg-config --libs squids nusquids hdf5` -lpthread

all: examples/exCross.o  examples/partial_rate_example examples/couplings_example examples/test

examples/exCross.o : include/exCross.h examples/exCross.cpp
	@ $(CXX) $(CFLAGS) -c examples/exCross.cpp -o $@

examples/partial_rate_example : examples/partial_rate_example.cpp examples/exCross.o include/nusquids_decay.h
	@echo Compiling partial_rate_example
	@ $(CXX) $(CFLAGS) examples/partial_rate_example.cpp examples/exCross.o  -o $@ $(LDFLAGS)

examples/couplings_example : examples/couplings_example.cpp examples/exCross.o include/nusquids_decay.h
	@echo Compiling couplings_example
	@ $(CXX) $(CFLAGS) examples/couplings_example.cpp examples/exCross.o  -o $@ $(LDFLAGS)

examples/test: examples/test.cpp examples/exCross.o include/nusquids_decay.h
	@echo Compiling test
	@ $(CXX) $(CFLAGS) examples/test.cpp examples/exCross.o  -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf ./examples/partial_rate_example ./examples/couplings_example ./examples/test  ./examples/exCross.o
