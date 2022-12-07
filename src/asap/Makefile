CC= gcc

CFLAGS= -O3  -Wall


all:	 asap


asap:	oldfns.c asap.c asap_common.c asap_core.c gdtosvg.c  draw.c
	$(CC) $(CFLAGS)  -o asap oldfns.c asap.c asap_common.c  asap_core.c gdtosvg.c  draw.c -lm


clean:
	\rm -f asap web_asap.cgi *.o

