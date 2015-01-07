CC = g++

INCLUDES = -I/usr/include -I~/software/boost_1_55_0/include
LFLAGS = -L~/software/boost_1_55_0/lib
OPTFLAGS = -O2
CFLAGS  = -g -Wall
TARGET = mstwrapper_rooted

#Set to -DDEBUG to compile with debug statements
DEBUG   = -DDEBUG


all: $(TARGET)

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).o $(LFLAGS)

$(TARGET).o: $(TARGET).cpp 
	$(CC) $(DEBUG) $(CFLAGS) $(INCLUDES) $(OPTFLAGS) -o $(TARGET).o -c $(TARGET).cpp

clean:
	$(RM) $(TARGET) *.o *.txt *.out
