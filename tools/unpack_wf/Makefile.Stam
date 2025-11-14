# Just change SHCI_DIR to the location of the SHCI directory
SHCI_DIR := ../..
FLAG := -qopenmp -qmkl -include stdio.h -include stdexcept -include cstdint
LIB_DIR := $(SHCI_DIR)/lib
PARENT_DIR := $(SHCI_DIR)/..
SRC := unpack_wf.cc
UTIL := $(SHCI_DIR)/src/util.cc
HDET := $(SHCI_DIR)/src/det/half_det.cc

all : unpack_wf.cc
	icpx -I$(LIB_DIR) -I$(PARENT_DIR) -I$(SHCI_DIR) $(UTIL) $(HDET) $(SRC) $(FLAG) -o unpack_wf
clean :
	$(RM) unpack_wf
