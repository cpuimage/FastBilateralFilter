CC = gcc -O3 $(FLAGS)

LIBS = -lm

SRCS = mt19937ar.c affinity.c arrayalloc.c noisycomputations.c imageio.c maxfilter.c young.c deriche_o3opt.c gnuplot_i.c fastbf.c
 
SRCS += fastbf_main.c

CFLAGS= -D_GNU_SOURCE

ifdef OMP
LIBS += -fopenmp
CFLAGS += -fopenmp
endif

OBJS = $(SRCS:.c=.o)

MAIN = FBF

.PHONY: depend clean rebuild

all:    $(MAIN) cleanobjs

$(MAIN): $(OBJS)
	$(CC) -o $(MAIN) $(OBJS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

cleanobjs:
	$(RM) *.o
clean:
	$(RM) *.o *~ $(MAIN)

rebuild: clean all

depend: $(SRCS)
	makedepend -I/ $^

# DO NOT DELETE THIS LINE -- make depend needs it
