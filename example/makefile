CFLAGS=`zm-config --cflags`
LINK=`zm-config -l`

TARGET=$(shell ls *.c | xargs -i basename {} .c | tr -s "\n" " ")

CC=gcc
CFLAGS+=-ansi -Wall -O3 -funroll-loops -std=c99 -g
# -std=c99 is for gamma function test.

all: $(TARGET)
%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LINK)
clean :
	-@rm -f *.o *~ core *test $(TARGET)
