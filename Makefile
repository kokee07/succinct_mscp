CPP=g++ -std=c++11
CPPFLAGS=-O3 -Wall -DVERBOSE
INCLUDES=-I./include/
SOURCES=./include/BasicCDS.cpp VIDCA.cpp
SOURCESS=./include/BasicCDS.cpp VIDCA_WIP.cpp
test=./include/BasicCDS.cpp vidca_test_perf.cpp
MEM=./include/BasicCDS.cpp testMemory.cpp
SOURCESSS=./include/BasicCDS.cpp test.cpp
SOURCESSSS=./include/BasicCDS.cpp VIDCA_FROM_FILE_5.cpp
unary:
	@echo " Building vidca Heuristic" $(SOURCESSSS)
	@$(CPP) $(CPPFLAGS) $(INCLUDES) $(SOURCESSSS) -o vidca3
