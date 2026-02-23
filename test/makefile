ZM_CFLAGS=`zm-config --cflags`
LINK=`zm-config -l`
LINKCXX=`zm-config -lcpp`

CC=gcc
CFLAGS=-Wall -O3 -funroll-loops -std=c99 ${ZM_CFLAGS}
# -std=c99 is for gamma function test.

CXX=g++
CXXFLAGS=-std=c++17 -Wall -O3 -funroll-loops ${ZM_CFLAGS}

%: %.c
	@$(CC) $(CFLAGS) -o $@ $< $(LINK)
%: %.cpp
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LINKCXX)
clean :
	-@rm -f *.o *~ core *test
