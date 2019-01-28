#
# Makefile -- FemtoFitter
#

BUILD_DIR ?= build

.PHONY: all config

all: $(BUILD_DIR)/Makefile
	+make -C $(BUILD_DIR)

config: $(BUILD_DIR)/Makefile
	ccmake $(BUILD_DIR)

clean:
	+make -C $(BUILD_DIR) clean

cleanall:
	rm -r $(BUILD_DIR)/*

help:
	@echo
	@echo " Meta-Make FemtoFitter package"
	@echo
	@echo "Usage:"
	@echo "  make all      # build with standard settings"
	@echo "  make config   # use ccmake to edit build configuration "
	@echo "  make clean    # clean output"
	@echo
	@echo 'Example Environment Variables:'
	@echo '  BUILD_DIR=$$HOME/.build/femtofitter make     # build in specific directory'
	@echo '  CPATH=/opt/root/include:/other/path make    # add extra include directories'
	@echo '  ROOTSYS=/opt/root/include make              # set ROOT location'
	@echo


$(BUILD_DIR)/Makefile: femtofitter/CMakeLists.txt
	cmake -S femtofitter -B $(BUILD_DIR)

docs-html:
	+make -C docs html

bwalk: macros/BayesianWalk.cpp
	#g++ -I. -Ifemtofitter/src $(shell root-config --cflags) -O3 -funroll-loops -ffast-math -march=native -fopenmp -o $@ $<  build/libFemtoFitter.so $(shell root-config --libs)  -lMinuit
	g++ -I. -Ifemtofitter/src $(shell root-config --cflags) -O3 -ffast-math -march=native -fopenmp -o $@ $< build/libFemtoFitter.so -L${ROOTSYS}/lib -lMinuit -lCore -lHist -lRIO -lGpad
