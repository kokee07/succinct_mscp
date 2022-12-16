CPP=g++ -std=c++11
CPPFLAGS=-O3 -Wall -DVERBOSE
INCLUDES=-I./include/
VEXHA=./include/BasicCDS.cpp VIDCA_EXHAUSTIVE.cpp
VHEU=./include/BasicCDS.cpp VIDCA_HEURISTICS.cpp
VNAT=./include/BasicCDS.cpp VIDCA_NATURAL.cpp
VSPP=./include/BasicCDS.cpp VIDCA_NATURAL_sppwn.cpp
VSTEP=./include/BasicCDS.cpp VIDCA_BIGSTEP.cpp
Vtest=./include/BasicCDS.cpp test_lectura.cpp
all: test
vidca_1: VIDCA_EXHAUSTIVE.cpp
	@echo " Building vidca Heuristic with solution provided"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VEXHA) -o VIDCA_EXHAUSTIVE

vidca_2: VIDCA_HEURISTICS.cpp
	@echo " Building vidca Heuristic with no solution provided: "
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VHEU) -o VIDCA_HEURISTICS

natural: VIDCA_NATURAL.cpp
	@echo " Building vidca Heuristic for natural sets"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VNAT) -o VIDCA_NATURAL

SPP: VIDCA_NATURAL_sppwn.cpp
	@echo " Building vidca Heuristic for natural sets SPP"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VSPP) -o VIDCA_SPP
bigstep: VIDCA_BIGSTEP.cpp
	@echo " Building vidca big step heuristic for natural sets"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(VSTEP) -o VIDCA_BIGSTEP

test: test_lectura.cpp
	@echo " Building test"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(Vtest) -o test
clean:
	-rm *~ *.o *.bak
