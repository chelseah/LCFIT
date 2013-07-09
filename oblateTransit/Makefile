OBJECT = oblateness.o
PROJECT = oblateness
all: oblateness


$(PROJECT): $(OBJECT)
	$(CXX) -pg -o $@ $^

clean:
	$(RM) $(PROJECT)
	$(RM) *.o
	$(RM) .depend
	$(RM) testopt
	$(RM) *.out
depend:
	    $(CXX) -MM $(CXXFLAGS) *.cpp > .depend


