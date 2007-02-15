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

# Dependencies

display.o: display.cc globals.hh utilities.hh surface.hh common.hh
keyboard.o: keyboard.cc display.hh globals.hh mesh_surface.hh \
            nurbs_surface.hh common.hh
mesh_surface.o: mesh_surface.cc display.hh globals.hh mesh_surface.hh \
                surface.hh common.hh
mouse.o: mouse.cc display.hh globals.hh common.hh
nurbs_surface.o: nurbs_surface.cc nurbs_surface.hh surface.hh globals.hh \
                 common.hh
sfview.o: sfview.cc display.hh globals.hh keyboard.hh mouse.hh common.hh
surface.o: surface.cc
utilities.o: utilities.cc common.hh
