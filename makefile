CONFIG:=$(shell test -e config || cp config.org config; echo config)
include $(CONFIG)

ROOTDIR:=.
INCDIR:=$(ROOTDIR)/include/$(PROJNAME)
SRCDIR:=$(ROOTDIR)/src
LIBDIR:=$(ROOTDIR)/lib
TOOLDIR:=$(ROOTDIR)/tools
APPDIR:=$(ROOTDIR)/app
DOCDIR:=$(ROOTDIR)/doc
TESTDIR:=$(ROOTDIR)/test
SAMPLEDIR:=$(ROOTDIR)/example

CHKDEP=`which zeda-chkdep`

all: library devtools
library:
	@$(CHKDEP) $(DEPENDENCY) || exit 1
	@cd $(SRCDIR); make
devtools:
	@cd $(TOOLDIR); make
application:
	@cd $(APPDIR); make
autotest:
	@cd $(TESTDIR); ./test.sh
doc:
	@cd $(DOCDIR); make
clean:
	-@rm -f $(ROOTDIR)/*~ $(INCDIR)/*~
	@cd $(SRCDIR); make clean
	-@rm -f $(LIBDIR)/*.so
	@cd $(TOOLDIR); make clean
	@cd $(APPDIR); make clean
	@cd $(DOCDIR); make clean
	@cd $(SAMPLEDIR); ./allclean.sh
install:
	@echo " INSTALL	library"
	-@install -m 755 $(LIBDIR)/*.so $(PREFIX)/lib/
	@echo " INSTALL	header files"
	-@install -m 755 -d $(PREFIX)/include/$(PROJNAME)
	-@install -m 644 $(INCDIR)/*.h $(PREFIX)/include/$(PROJNAME)/
	@echo " INSTALL	tools"; cd $(TOOLDIR); make install
	@echo " INSTALL	applications"; cd $(APPDIR); make && make install
uninstall:
	@echo " UNINSTALL	library"
	-@rm $(PREFIX)/lib/lib$(PROJNAME).so
	@echo " UNINSTALL	header files"
	-@rm -r $(PREFIX)/include/$(PROJNAME)
	@echo " UNINSTALL	tools"; cd $(TOOLDIR); make uninstall
	@echo " UNINSTALL	applications"; cd $(APPDIR); make uninstall
