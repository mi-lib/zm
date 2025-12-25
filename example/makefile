CC=gcc
CFLAGS=-ansi -Wall -O3 -funroll-loops -std=c99 -g `zm-config --cflags`
LINK=`zm-config -l`
# -std=c99 is for gamma function test.

CXX=g++
CXXFLAGS=-std=c++17 -Wall -O3 -funroll-loops -g `zm-config --cflags`
LINKCXX=`zm-config -lcpp`

TARGET=$(shell ls *.c | xargs -I{} basename {} .c | tr -s "\n" " ") $(shell ls *.cpp | xargs -I{} basename {} .cpp | tr -s "\n" " ")

all: $(TARGET)
%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LINK)
%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINKCXX)
clean :
	-@rm -f *.o *~ core *test $(TARGET)
