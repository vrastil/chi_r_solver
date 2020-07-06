SRC := NFW_Chameleon.cpp integrator.cpp constants.cpp io.cpp
CC := g++
OPT := -Ofast -fPIC -std=c++11 -pthread
MAIN := main.a
LIB := -lboost_program_options -lboost_filesystem -lboost_system -lboost_log

all:
	$(CC) $(SRC) $(OPT) $(LIB) -o $(MAIN)