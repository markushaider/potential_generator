OPT1	= -DZA=128
OPT2    = -DAllBoxSize=32000
OPT3	= -DSHIFT_X=0
OPT4	= -DSHIFT_Y=0
OPT5	= -DSHIFT_Z=0
OPT6    = -DHYDROBOX=20000

CC       =  gcc    # sets the C-compiler (default)
OPTIMIZE =   -O2 -Wall -msse2 -mmmx   # optimization and warning flags (default)

FFTW_INCL= -I/usr/local/include/
FFTW_LIBS= -L/usr/local/lib/


FFTW_LIB = $(FFTW_LIBS) -lfftw3
OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6)

EXEC   = Potcomo

OBJS   = main.o dichterechner.o reader.o potentialhaupt.o fftwdichte.o fftwkernel.o potentialverschachtler.o fitsschreiberling.o

INCL   =   proto.h variablen.h
CFLAGS =   $(OPTIONS)  $(FFTW_INCL)


LIBS   =  -lm $(FFTW_LIB)

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
 
$(OBJS): $(INCL) 

#######################################################################
.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)
