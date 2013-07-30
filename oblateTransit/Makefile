OBJECT = oblateness.o main.o
PROJECT = oblateness
OBL_INTERFACES = oblateness.i
ELL_INTERFACES = elliptic.i
OBL_SOURCES= oblateness.cc oblateness.h
ELL_SOURCES = elliptic.cc elliptic.h
SOURCES = $(OBL_SOURCES) $(ELL_SOURCES)
INTERFACES = $(OBL_INTERFACES) $(ELL_INTERFACES)
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

.PHONY: clean
clean:
	$(RM) $(PROJECT)
	$(RM) *.so *.pyc *_wrap.h
	$(RM) -r *.dSYM
	$(RM) -r build
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
	$(RM) *.o
	$(RM) .depend
	$(RM) *.out
	$(RM) *~
