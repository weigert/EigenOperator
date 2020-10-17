################################################################################
# EigenOperator Test Makefile
################################################################################

# Parameters
CC = g++ -std=c++14
FLAGS = -march=native -Wall -g -rdynamic -O3
LIBS	= -I/usr/local/include -L/usr/local/lib -pthread

# Build Targets
A_SRC = main.cpp
A_TARGET  = operator

# Make Rules
.PHONY: all
all: operator

.PHONY: operator
operator: $(A_SRC)
	$(CC) $(A_SRC) $(FLAGS) $(LIBS) -o $(A_TARGET)
