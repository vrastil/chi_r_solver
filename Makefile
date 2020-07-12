CXX	 		:= g++
OPT 		:= -fPIC -std=c++11 -pthread -flto=jobserver -MMD -MP
ARCH	 	:= native
OPTIMIZE	:= -Ofast -march=${ARCH} -mtune=${ARCH} -Wall
CXXFLAGS 	:= $(OPT) $(OPTIMIZE)
LDFLAGS 	:= -lboost_program_options -lboost_filesystem -lboost_system -lboost_log

BUILD   	:= ./build
OBJ_DIR		:= $(BUILD)/objects
APP_DIR		:= $(BUILD)
TARGET 		:= main.a
INCLUDE	  	:= -Iinclude/
SRC 		:= src/main.cpp src/integrator.cpp src/io.cpp src/stars_NFW.cpp src/units.cpp

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPS 	 := $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	+$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	+$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

.PHONY: all build clean

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

-include $(DEPS)
