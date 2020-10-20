CXX = g++ 
CXXFLAGS = -std=c++17 -Wall -g
IPTH = -lpthread


COMMON= BFS-PGCn.cpp BFS-PGICn.cpp  SSSP-PGCn.cpp SSSP-PGICn.cpp BC-PGCn.cpp BC-PGICn.cpp 

ALL= main-BFS-PGCn main-BFS-PGICn main-SSSP-PGCn main-SSSP-PGICn main-BC-PGCn main-BC-PGICn 

all: $(ALL)

$(ALL): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(IPTH)

.PHONY : clean

clean :
	rm -f *.o $(ALL)
	

