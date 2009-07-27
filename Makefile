TARGETS=sfview sfview.html sfview.ps TAGS

all: $(TARGETS)

SOURCES=sfview.cc glwindow.cc utilities.cc surface.cc \
        mesh-surface.cc nurbs-surface.cc

OBJECTS=$(subst .cc,.o,$(SOURCES))

DEPENDENCIES=$(subst .cc,.d,$(SOURCES))

CXXFLAGS=-g -Wall

LDFLAGS=-lglut -lpng

sfview: $(OBJECTS)

sfview.html: sfview.1
	groff -man -Thtml $< > $@

sfview.ps: sfview.1
	groff -man -Tps $< > $@

include $(DEPENDENCIES)

# Rule for generating dependency files, from the make manual.
%.d: %.cc
	$(CXX) -M $(CXXFLAGS) $< > $@.$$$$;                 \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	$(RM) $@.$$$$

.PHONY: TAGS
TAGS:
	find . -name "*.[ch][ch]" -print | etags -

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(TARGETS) $(DEPENDENCIES)
