######################################################################
# Project: TMVA - Toolkit for Multivariate Data Analysis             #
# Code   : Example usage for TMVA analyses                           #
###################################################################### 
BINS = \
	TMVAClassification_ele \
	TMVAClassification_pho \


CXX = g++
CCFLAGS = $(shell root-config --cflags) -ggdb -fopenmp -O3 -g


LD = g++
LDFLAGS = 

INCLUDE = -I./

LIBS = $(shell root-config --libs) -lMLP -lMinuit -lTreePlayer -lTMVA -lTMVAGui -lXMLIO  -lMLP -lm


default : $(BINS)

$(BINS): % : %.C 
	@echo -n "Building $@ ... "
	$(CXX) $(CCFLAGS) $(INCLUDE) $< $(LIBS) -o $@
	@echo "Done"


clean: 
	rm -f *.o	$(BINS)

