CPP=g++ -std=c++11
CPPFLAGS=-O3 -Wall -DVERBOSE
INCLUDES=-I./include/
VEXHA=./include/BasicCDS.cpp VIDCA_EXHAUSTIVE.cpp
VHEU=./include/BasicCDS.cpp VIDCA_HEURISTICS.cpp
VNAT=./include/BasicCDS.cpp VIDCA_NATURAL.cpp
VPAR=./include/BasicCDS.cpp mscp_parallelism.cpp

all: clean vidca_1 vidca_2 natural parallel
vidca_1: VIDCA_EXHAUSTIVE.cpp
	@echo " Building vidca Heuristic with solution provided"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VEXHA) -o VIDCA_EXHAUSTIVE

vidca_2: VIDCA_HEURISTICS.cpp
	@echo " Building vidca Heuristic with no solution provided: "
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VHEU) -o VIDCA_HEURISTICS

natural: VIDCA_NATURAL.cpp
	@echo " Building vidca Heuristic for natural sets"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VNAT) -o VIDCA_NATURAL

parallel:
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VPAR) -fopenmp -o runMSCP
clean:
	-rm *~ *.o *.bak
