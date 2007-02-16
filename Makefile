TARGETS=sfview sfview.ps

all: $(TARGETS)

SOURCES=display.cc \
	keyboard.cc \
	mesh_surface.cc \
	mouse.cc \
	nurbs_surface.cc \
	sfview.cc \
	surface.cc \
	utilities.cc

OBJECTS=$(subst .cc,.o,$(SOURCES))

DEPENDENCIES=$(subst .cc,.d,$(SOURCES))

CXXFLAGS=-g -Wall

LDFLAGS=-lglut

sfview: $(OBJECTS)

sfview.ps: sfview.1
	groff -man -Tps $< > $@

include $(DEPENDENCIES)

# Rule for generating dependency files, from the make manual.

%.d: %.cc
	$(CXX) -M $(CXXFLAGS) $< > $@.$$$$;                 \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	$(RM) $@.$$$$

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(TARGETS) $(DEPENDENCIES)
