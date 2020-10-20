ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

#compilers
ifdef CILK
PCC = g++ -std=c++17
PCFLAGS = -fcilkplus -lcilkrts -g -DCILK $(INTT) $(INTE)

else ifdef MKLROOT
PCC = icpc
PCFLAGS = -g -DCILKP $(INTT) $(INTE)


else ifdef OPENMP
PCC = g++ -std=c++17
PCFLAGS = -fopenmp -DOPENMP $(INTT) $(INTE)

else
PCC = g++ -std=c++17
PCFLAGS = -g $(INTT) $(INTE)
endif

COMMON= ligra.h graph.h utils.h IO.h parallel.h gettime.h quickSort.h parseCommandLine.h

ALL= BFSLigra SSSPLigra BCLigra

all: $(ALL)

% : %.C $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $< 

.PHONY : clean

clean :
	rm -f *.o $(ALL)
