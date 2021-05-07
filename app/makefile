include ../config

CONFIG=zm-config
BINDIR=$(PREFIX)/bin

INCLUDE=`$(CONFIG) -I`
LIB=`$(CONFIG) -L`
DEF=`$(CONFIG) -D`
LINK=`$(CONFIG) -l`

CC=gcc
CFLAGS=-ansi -Wall -O3 $(INCLUDE) $(LIB) $(DEF)

COMPILE=$(CC) $(CFLAGS) -o $@ $< $(LINK)
INSTALL=install -m 755
UNINSTALL=rm -f

TARGET=deg2rad rad2deg le psolve pexp seqmod

#####################################################################
all: $(TARGET)

%: %.c
	@echo " CC	" $<
	@$(COMPILE)
clean:
	@echo " CLEAN"
	-@rm -f *.o *~ core $(TARGET)
install: $(TARGET)
ifneq ($(TARGET),)
	@echo " INSTALL	" $^
	@$(INSTALL) $^ $(BINDIR)
endif
uninstall:
ifneq ($(TARGET),)
	@echo " UNINSTALL	" $(TARGET)
	@$(UNINSTALL) $(BINDIR)/$(TARGET)
endif
