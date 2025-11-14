# Just change SHCI_DIR to the location of the SHCI directory
SHCI_DIR := ../..
OBJ_DIR := $(SHCI_DIR)/build #added 
FLAG := -qopenmp -qmkl -include stdio.h -include stdexcept -include cstdint
LIB_DIR := $(SHCI_DIR)/lib
PARENT_DIR := $(SHCI_DIR)/..
SRC := wfc_proj.cc
UTIL := $(SHCI_DIR)/src/util.cc
HDET := $(SHCI_DIR)/src/det/half_det.cc

all : wfc_proj.cc
	icpx -I$(LIB_DIR) -I$(PARENT_DIR) -I$(SHCI_DIR) $(UTIL) $(HDET) $(SRC) $(FLAG) -o wfc_proj
clean :
	$(RM) wfc_proj
