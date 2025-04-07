
BDSIM_DIR=/usr/local/include/bdsim
G4_DIR=/usr/local/include/Geant4
BASE_DIR=${PWD}

CXX = g++
INCLUDES=-I$(BASE_DIR)/include -I$(BDSIM_DIR) -I$(BDSIM_DIR)/analysis -I$(BDSIM_DIR)/parser -I$(G4_DIR) `root-config --glibs --cflags`
CXXFLAGS = -fPIC -Wall -Wextra -g -Wno-sign-compare
LDFLAGS = -shared

BDSIM_LIBS=-lbdsim -lrebdsim

TARGET = libinterface.so

PROFILER=-lprofiler

LIB_DIR=libs
BIN_DIR=bin

INTERFACE_DIR = Interface
INTERFACE_FILES = $(wildcard $(INTERFACE_DIR)/*.cc)
OBJ_FILES = $(INTERFACE_FILES:.cc=.o)

SRC_DIR=src
SRC_FILES=$(wildcard $(SRC_DIR)/*.cc)
EXECUTABLES=$(SRC_FILES:$(SRC_DIR)/%.cc=$(BIN_DIR)/%)


all: $(LIB_DIR)/$(TARGET)  $(EXECUTABLES)

$(LIB_DIR)/$(TARGET): $(OBJ_FILES)
	$(CXX) $(LDFLAGS) $(BDSIM_LIBS) $(INCLUDES) $(PROFILER) -o $@ $^

$(INTERFACE_DIR)/%.o: $(INTERFACE_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(PROFILER) -c $< -o $@


$(BIN_DIR)/%: $(SRC_DIR)/%.cc $(LIB_DIR)/$(TARGET)
	$(CXX) $(CXXFLAGS) -L$(LIB_DIR) $(INCLUDES) $(PROFILER) $(BDSIM_LIBS) -linterface -o $@ $< -I./include


clean:
	rm libs/*





