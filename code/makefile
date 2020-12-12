Compile=g++ -g -std=c++17
dtmc_lc: main.o dtmc_lc.o metropolis.o measurement.o experiment.o
	${Compile} -o dtmc_lc main.o dtmc_lc.o metropolis.o  measurement.o experiment.o
main.o: main.cpp dtmc_lc.h
	${Compile} -o main.o -c main.cpp
dtmc_lc.o: dtmc_lc.cpp dtmc_lc.h
	${Compile} -o dtmc_lc.o -c dtmc_lc.cpp
metropolis.o: metropolis.cpp dtmc_lc.h
	${Compile} -o metropolis.o -c metropolis.cpp
measurement.o: measurement.cpp dtmc_lc.h
	${Compile} -o measurement.o -c measurement.cpp
experiment.o: experiment.cpp dtmc_lc.h
	${Compile} -o experiment.o -c experiment.cpp
clean:
	rm *.o
