
PROG=reptile-v1.1

CC=g++
#CFLAGS=-c -Wall -O3  
CFLAGS=-c -g -Wall 
LDFLAGS= -lz 
SRCS=$(wildcard *.cpp) 
OBJS=$(SRCS:.cpp=.o)

all: $(SRCS) $(PROG)
		
$(PROG): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $@


$(OBJS): %.o : %.cpp
	$(CC) $(CFLAGS) $<   

clean:
	rm -rf $(PROG) $(OBJS)
