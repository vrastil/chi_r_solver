CXX	 		:= g++
OPT 		:= -fPIC -std=c++11 -pthread -flto=jobserver
ARCH	 	:= native
OPTIMIZE	:= -Ofast -march=${ARCH} -mtune=${ARCH}
CXXFLAGS 	:= $(OPT) $(OPTIMIZE)
LDFLAGS 	:= -lboost_program_options -lboost_filesystem -lboost_system -lboost_log

BUILD   	:= ./build
OBJ_DIR		:= $(BUILD)/objects
APP_DIR		:= $(BUILD)
TARGET 		:= main.a
INCLUDE	  	:= -Iinclude/
SRC 		:= src/main.cpp src/integrator.cpp src/io.cpp src/stars_NFW.cpp

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

.PHONY: all build clean

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

