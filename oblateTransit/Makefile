OBJECT = oblatenessfast.o main.o
PROJECT = oblateness
OBJECT_DEBUG = oblatenessfast.o testinter.o
PROJECT_DEBUG = testinter
OBL_INTERFACES = oblateness.i
OBLFAST_INTERFACES = oblatenessfast.i
ELL_INTERFACES = elliptic.i
OBL_SOURCES= oblateness.cc oblateness.h 
OBLFAST_SOURCES= oblatenessfast.cc oblatenessfast.h
ELL_SOURCES = elliptic.cc elliptic.h
SOURCES = $(OBL_SOURCES) $(ELL_SOURCES) $(OBLFAST_INTERFACES)
INTERFACES = $(OBL_INTERFACES) $(ELL_INTERFACES) $(OBLFAST_INTERFACES)
WRAPPERS = $(INTERFACES:.i=_wrap.cxx)
PROXIES = $(INTERFACES:.i=.py)

.PHONY: all
all: $(WRAPPERS) $(SOURCES)
	./setup.py build_ext -i

%_wrap.cxx: %.i %.h ./numpy.i
	swig -c++ -python $<


.PHONY: cc
cc: $(PROJECT)

$(PROJECT) :$(OBJECT)
	$(CXX) -pg -o $@ $^

#.PHONY: test
#test: $(PROJECT_DEBUG)
#
#$(PROJECT_DEBUG) :$(OBJECT_DEBUG)
#	$(CXX) -pg -o $@ $^


.PHONY: clean
clean:
	$(RM) $(PROJECT)
	$(RM) $(PROJECT_DEBUG)
	$(RM) *.so *.pyc *_wrap.h
	$(RM) -r *.dSYM
	$(RM) -r build
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
	$(RM) *.o
	$(RM) .depend
	$(RM) *.out
	$(RM) *~
