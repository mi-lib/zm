CFLAGS=`zm-config --cflags`
LINK=`zm-config -l`

CC=gcc
CFLAGS+=-Wall -O3 $(LIB) $(INCLUDE) -funroll-loops -std=c99
# -std=c99 is for gamma function test.

%: %.c
	@$(CC) $(CFLAGS) -o $@ $< $(LINK)
clean :
	-@rm -f *.o *~ core *test
