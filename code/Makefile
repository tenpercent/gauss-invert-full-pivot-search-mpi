CXX = mpicxx
ARGS = -O0 -lmpi -lm
ARGS_FAST = -O3 -lmpi -lm -Wall
FILES = main.cpp io.cpp tools.cpp gauss-invert.cpp multiplication.cpp simpleinvert.cpp

all: release debug
release:
	$(CXX) $(ARGS_FAST) $(FILES) -o test.out
debug:
	$(CXX) -g $(ARGS) $(FILES) -o testd.out
clean:
	rm *.out 
archive:
	tar -cf backup.tar .
