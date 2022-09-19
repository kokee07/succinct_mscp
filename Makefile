CPP=g++ -std=c++11
CPPFLAGS=-O3 -Wall -DVERBOSE
INCLUDES=-I./include/
VEXHA=./include/BasicCDS.cpp VIDCA_EXHAUSTIVE.cpp
VHEU=./include/BasicCDS.cpp VIDCA_HEURISTIC.cpp
all: vidca_1 vidca_2

vidca_1: VIDCA_EXHAUSTIVE.cpp
	@echo " Building vidca Heuristic with solution provided"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VEXHA) -o VIDCA_EXHAUSTIVE

vidca_2: VIDCA_HEURISTIC.cpp
	@echo " Building vidca Heuristic with no solution provided"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VHEU) -o VIDCA_HEURISTIC


clean:
	-rm *~ *.o *.bak
