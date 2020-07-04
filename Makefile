SRC := NFW_Chameleon.cpp integrator.cpp constants.cpp # io.cpp
CC := g++
OPT := -Ofast -fPIC
MAIN := main.a

all:
	$(CC) $(SRC) -o $(MAIN)