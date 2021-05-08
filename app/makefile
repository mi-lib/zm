MAKEFILEGEN=`which zeda-makefile-gen`

all:
	@$(MAKEFILEGEN) app -ansi | make -f -
clean:
	@$(MAKEFILEGEN) app -ansi | make -f - clean
install:
	@$(MAKEFILEGEN) app -ansi | make -f - install
uninstall:
	@$(MAKEFILEGEN) app -ansi | make -f - uninstall
