ROOTDIR=$(shell pwd)
OUTPUTDIR=$(ROOTDIR)/paraview-output/

.PHONY: all cleanall clean clean_paraview
all: assignment-gcc assignment-icpc

# Target to be used with the GNU Compiler Collection.
assignment-gcc: CXX=g++
assignment-gcc: CXXFLAGS=-fopenmp -O3 -march=native -std=c++0x
assignment-gcc: assignment-code.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# Target to be used with the Intel C++ compiler.
# In order to use this compiler on Hamilton, you should first add the
# corresponding module with
#     $ module add intel/2020.4
assignment-icpc: CXX=icpc
assignment-icpc: CXXFLAGS=-O3 -xHost -std=c++0x
assignment-icpc: assignment-code.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

cleanall: clean clean_paraview

clean:
	rm -rf $(ROOTDIR)/assignment-gcc $(ROOTDIR)/assignment-icpc

clean_paraview:
	find $(OUTPUTDIR) -iname "result-*.vtp" -delete
	rm -rf $(OUTPUTDIR)/result.pvd

test:
	./validate.sh
	python3 validate.py
