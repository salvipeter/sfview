TARGETS=sfview sfview.ps sfview.txt

all: $(TARGETS)

OBJECTS=display.o \
	keyboard.o \
	mesh_surface.o \
	mouse.o \
	nurbs_surface.o \
	sfview.o \
	surface.o \
	utilities.o

CXXFLAGS=-g -Wall

LDFLAGS=-lglut

sfview: $(OBJECTS)

sfview.ps: sfview.1
	groff -man -Tps $< > $@

sfview.txt: sfview.1
	groff -man -Tascii $< > $@

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(TARGETS)
