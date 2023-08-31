IDIR = include
CC = g++
CFLAGS = -I$(IDIR)

SDIR = src
ODIR = $(SDIR)/obj

_DEPS = classes.h helpers.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
_OBJ = classes.o helpers.o RDS.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

RDS.x: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~