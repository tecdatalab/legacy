# zorro
target   = zorro
cc_objects = read.o main.o
objects  = $(cc_objects)
CXX      = g++
DEPFLAGS = -MM 
CXXFLAGS = -O -Wall
LIBS     = -lm
LDFLAGS  =

# everything below here is automatic
$(target) : $(objects)
	$(CXX) $(CXXFLAGS) -o $@ $(objects) $(LDFLAGS) $(LIBS)

depends = $(objects:.o=.d)

-include $(depends)

%.d : %.cc
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) $< \
      | sed -e 's,\($*\)\.o:,\1.o $@ : ,g' > $@


.PHONY : clean
clean :
	-rm -f $(target) $(objects) $(depends)



