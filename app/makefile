include ../config

INCLUDE=-I$(PREFIX)/include
LIB=-L$(PREFIX)/lib
BINDIR=$(PREFIX)/bin

CC=gcc
CFLAGS=-ansi -Wall -O3 $(INCLUDE) $(LIB)

LINK=-lzm `zeda-config -l`

COMPILE=$(CC) $(CFLAGS) -o $@ $< $(LINK)

INSTALL=install -m 755
UNINSTALL=rm -f

CONFIGGEN=zeda-config-gen
CONFIG=zm-config
TARGET=$(CONFIG) deg2rad rad2deg le psolve pexp seqmod

#####################################################################
all: $(TARGET)

$(CONFIG):
	@echo " GENERATE" $@
	-@rm -f $@
	@$(CONFIGGEN) -I "$(INCLUDE)" -L "$(LIB)" -l "$(LINK)" -v $(VERSION) > $@
%: %.c
	@echo " CC	" $<
	@$(COMPILE)
install: $(TARGET)
	@echo " INSTALL	" $^
	@$(INSTALL) $^ $(BINDIR)
clean:
	@echo " CLEAN"
	-@rm -f *.o *~ core $(TARGET)
uninstall:
	@echo " UNINSTALL	" $(TARGET)
	@cd $(BINDIR); $(UNINSTALL) $(TARGET); cd -
