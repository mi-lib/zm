MAKEFILEGEN=`which zeda-makefile-gen`
MAKEDEB=`which zeda-deb-gen`

.PHONY: doc test example

all:
ifeq ($(MAKEFILEGEN),)
	echo "ZEDA not installed."
else
	@$(MAKEFILEGEN) | make -f -
endif
test:
	@$(MAKEFILEGEN) | make -f - test
doc:
	@$(MAKEFILEGEN) | make -f - doc
example:
	@$(MAKEFILEGEN) | make -f - example
clean:
	@$(MAKEFILEGEN) | make -f - clean
install:
	@$(MAKEFILEGEN) | make -f - install
uninstall:
	@$(MAKEFILEGEN) | make -f - uninstall
deb:
ifeq ($(MAKEDEB),)
	echo "ZEDA is not installed."
else
	@$(MAKEDEB)
endif
