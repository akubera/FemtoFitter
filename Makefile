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
	make -C $(BUILD_DIR) clean

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


$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake ..

