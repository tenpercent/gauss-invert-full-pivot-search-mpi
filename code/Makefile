CXX = mpicxx
ARGS = -O0 -lm
ARGS_FAST = -O3 -lm -Wall
FILES = main.cpp io.cpp tools.cpp gauss-invert.cpp multiplication.cpp simpleinvert.cpp residual.cpp

all: release debug
release:
	$(CXX) $(ARGS_FAST) $(FILES) -o test.out
debug:
	$(CXX) -g $(ARGS) $(FILES) -o testd.out
clean:
	rm *.out 
archive:
	tar -cf backup.tar .
