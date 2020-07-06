SRC := NFW_Chameleon.cpp integrator.cpp constants.cpp io.cpp
CC := g++
OPT := -fPIC -std=c++11 -pthread -flto=jobserver
ARCH := native
OPTIMIZE := -Ofast -march=${ARCH} -mtune=${ARCH}
MAIN := main.a
LIB := -lboost_program_options -lboost_filesystem -lboost_system -lboost_log

all:
	$(CC) $(SRC) $(OPTIMIZE) $(OPT) $(LIB) -o $(MAIN)