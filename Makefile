# FLAGS

CFLAGS= -O3 -fPIC -std=c++11 -g
CFLAGS+= -I./inc `pkg-config --cflags nusquids`
LDFLAGS+= `pkg-config --libs nusquids`

ifeq ($(UNAME_S),Linux)
	LDFLAGS+=-lsupc++
endif

all: mains/exCross.o mains/run_musquids mains/run_musquids_simple mains/run_musquids_less_simple

mains/exCross.o : inc/exCross.h mains/exCross.cpp
	@$(CXX) $(CFLAGS) -c mains/exCross.cpp -o $@

mains/run_musquids : mains/run_musquids.cpp mains/exCross.o
	@echo Compiling run_musquids
	@$(CXX) $(CFLAGS) mains/run_musquids.cpp mains/exCross.o $(LDFLAGS) -o $@

mains/run_musquids_simple : mains/run_musquids_simple.cpp mains/exCross.o
	@echo Compiling run_musquids_simple
	@$(CXX) $(CFLAGS) mains/run_musquids_simple.cpp mains/exCross.o $(LDFLAGS) -o $@

mains/run_musquids_less_simple : mains/run_musquids_less_simple.cpp mains/exCross.o
	@echo Compiling run_musquids_less_simple
	@$(CXX) $(CFLAGS) mains/run_musquids_less_simple.cpp mains/exCross.o $(LDFLAGS) -o $@

.PHONY: clean
clean:
	rm -rf ./mains/run_musquids ./mains/run_musquids_simple ./mains/run_musquids_less_simple



