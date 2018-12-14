include ../config

DOX = doxygen
DOXCONF = dox.conf

TARGET = doc

all:
	@( echo "PROJECT_NUMBER = "$(VERSION);\
           echo "INCLUDE_PATH = ../include "$(PREFIX)/include )\
         | cat - $(DOXCONF) | $(DOX) -
clean:
	-@rm -fr *.o *~ \#* html man
