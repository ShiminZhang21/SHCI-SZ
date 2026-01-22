# Just change SHCI_DIR to the location of the SHCI directory
SHCI_DIR := ../..
FLAG := -qopenmp -qmkl -include stdio.h -include stdexcept -include cstdint -include cmath
LIB_DIR := $(SHCI_DIR)/lib
PARENT_DIR := $(SHCI_DIR)/..
SRC := calc_corr.cc
UTIL := $(SHCI_DIR)/src/util.cc
HDET := $(SHCI_DIR)/src/det/half_det.cc

all : calc_corr.cc
	icpx -I$(LIB_DIR) -I$(PARENT_DIR) -I$(SHCI_DIR) $(UTIL) $(HDET) $(SRC) $(FLAG) -o calc_corr
clean :
	$(RM) calc_corr
