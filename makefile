MAKEFILEGEN=`which zeda-makefile-gen`

all:
ifeq ($(MAKEFILEGEN),)
	echo "ZEDA not installed."
else
	@$(MAKEFILEGEN) | make -f -
endif
autotest:
	@$(MAKEFILEGEN) | make -f - autotest
doc:
	@$(MAKEFILEGEN) | make -f - doc
clean:
	@$(MAKEFILEGEN) | make -f - clean
install:
	@$(MAKEFILEGEN) | make -f - install
uninstall:
	@$(MAKEFILEGEN) | make -f - uninstall
