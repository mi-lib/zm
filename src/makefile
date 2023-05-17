MAKEFILEGEN=`which zeda-makefile-gen`

include ../config

EXPORT_HEADER=../include/${PROJNAME}/${PROJNAME}_export.h
PROJNAME_CAPITAL=$(shell echo ${PROJNAME} | tr 'a-z' 'A-Z' )

all:
	@rm -f ${EXPORT_HEADER}
	@touch ${EXPORT_HEADER}
	@echo "/* ${PROJNAME}_export.h */" >> ${EXPORT_HEADER}
	@echo "/* This file was automatically generated. */" >> ${EXPORT_HEADER}
	@echo "/* `date` by `whoami` */" >> ${EXPORT_HEADER}
	@echo "#ifndef __${PROJNAME_CAPITAL}_EXPORT_H__" >> ${EXPORT_HEADER}
	@echo "#define __${PROJNAME_CAPITAL}_EXPORT_H__" >> ${EXPORT_HEADER}
	@echo "#if defined(__WINDOWS__) && !defined(__CYGWIN__)" >> ${EXPORT_HEADER}
	@echo "# if defined(__${PROJNAME_CAPITAL}_BUILD_DLL__)" >> ${EXPORT_HEADER}
	@echo "#  define __${PROJNAME_CAPITAL}_EXPORT __declspec(dllexport)" >> ${EXPORT_HEADER}
	@echo "# else" >> ${EXPORT_HEADER}
	@echo "#  define __${PROJNAME_CAPITAL}_EXPORT __declspec(dllimport)" >> ${EXPORT_HEADER}
	@echo "# endif" >> ${EXPORT_HEADER}
	@echo "#else" >> ${EXPORT_HEADER}
	@echo "# define __${PROJNAME_CAPITAL}_EXPORT __EXPORT" >> ${EXPORT_HEADER}
	@echo "#endif" >> ${EXPORT_HEADER}
	@echo "#endif /* __${PROJNAME_CAPITAL}_EXPORT_H__ */" >> ${EXPORT_HEADER}
	@$(MAKEFILEGEN) lib | make -f -
%.o: %.c
	@$(MAKEFILEGEN) lib | make -f - $@
clean:
	@$(MAKEFILEGEN) lib | make -f - clean
	@rm -f ${EXPORT_HEADER}
