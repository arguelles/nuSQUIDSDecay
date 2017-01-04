# FLAGS

CFLAGS= -O3 -fPIC -std=c++11
#CFLAGS+= -g -fsanitize=address
CFLAGS+= -I./inc `pkg-config --cflags squids nusquids hdf5`
LDFLAGS+= `pkg-config --libs squids nusquids hdf5`

#all: mains/exCross.o mains/run_decay mains/run_decay_constant_density mains/run_decay_constant_density_with_decay
all: mains/exCross.o mains/run_decay 

mains/exCross.o : inc/exCross.h mains/exCross.cpp
	@ $(CXX) $(CFLAGS) -c mains/exCross.cpp -o $@

mains/run_decay : mains/run_decay.cpp mains/exCross.o inc/nusquids_decay.h
	@echo Compiling run_decay
	@ $(CXX) $(CFLAGS) mains/run_decay.cpp mains/exCross.o  -o $@ $(LDFLAGS)

#mains/run_decay_constant_density : mains/run_decay_constant_density.cpp mains/exCross.o inc/nusquids_decay.h
#	@echo Compiling run_decay_constant_density
#	@ $(CXX) $(CFLAGS) mains/run_decay_constant_density.cpp mains/exCross.o  -o $@ $(LDFLAGS)

#mains/run_decay_constant_density_with_decay : mains/run_decay_constant_density_with_decay.cpp mains/exCross.o inc/nusquids_decay.h
#	@echo Compiling run_decay_constant_density_with_decay
#	@ $(CXX) $(CFLAGS) mains/run_decay_constant_density_with_decay.cpp mains/exCross.o  -o $@ $(LDFLAGS)
#
.PHONY: clean
clean:
	rm -rf ./mains/run_decay ./mains/exCross.o
#	rm -rf ./mains/run_decay ./mains/run_decay_constant_density_with_decay ./mains/run_decay_constant_density ./mains/exCross.o



