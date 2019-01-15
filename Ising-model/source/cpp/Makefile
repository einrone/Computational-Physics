CC = mpic++
CFLAGS = -g -Wall
TARGET = ising
SRC = main.cpp dumpfiles.cpp quantities.cpp solver.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

clean:
	$(RM) $(TARGET)

program: main.o dumpfiles.o execute_solve.o quantities.o solver.o
	mpic++ -o ising main.o dumpfiles.o execute_solve.o quantities.o solver.o
