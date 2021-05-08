MAKEFILEGEN=`which zeda-makefile-gen`

all:
	@$(MAKEFILEGEN) tools | make -f -
install:
	@$(MAKEFILEGEN) tools | make -f - install
uninstall:
	@$(MAKEFILEGEN) tools | make -f - uninstall
clean:
	@$(MAKEFILEGEN) tools | make -f - clean
