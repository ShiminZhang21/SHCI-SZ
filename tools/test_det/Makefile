# Just change SHCI_DIR to the location of the SHCI directory
SHCI_DIR := ../..
OBJ_DIR := $(SHCI_DIR)/build_cray #added 

LIB_DIR := $(SHCI_DIR)/lib
PARENT_DIR := $(SHCI_DIR)/..
SRC := test_det.cc
UTIL := $(SHCI_DIR)/src/util.cc
HDET := $(SHCI_DIR)/src/det/half_det.cc

all : test_det.cc
	g++ -I$(LIB_DIR) -I$(PARENT_DIR) $(UTIL) $(HDET) $(SRC) -fopenmp -o test_det
#	g++ -I$(LIB_DIR) -I$(PARENT_DIR) $(UTIL) $(SRC) $(OBJ_DIR)/det/half_det.o -fopenmp -o test_det
clean :
	$(RM) test_det
