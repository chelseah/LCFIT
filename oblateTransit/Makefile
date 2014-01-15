OBJECT = oblatenessfast.o main.o
PROJECT = oblateness
OBJECT_DEBUG = oblatenessfast.o testinter.o
PROJECT_DEBUG = testinter
OBL_INTERFACES = oblateness.i
OBLFAST_INTERFACES = oblatenessfast.i
ELL_INTERFACES = elliptic.i
ZEIPEL_INTERFACES =  quadpolar.i quadsphere.i Zeipelmodel.i Zeipel.i
OBL_SOURCES= oblateness.cc oblateness.h 
OBLFAST_SOURCES= oblatenessfast.cc oblatenessfast.h
ZEIPEL_SOURCES = quadpolar.cc quadsphere.cc Zeipelmodel.cc Zeipel.cc Zeipelmodel.h Zeipel.h model.h quadpolar.h quadsphere.cc integrator.h
ELL_SOURCES = elliptic.cc elliptic.h
SOURCES = $(OBL_SOURCES) $(ELL_SOURCES) $(OBLFAST_INTERFACES) $(ZEIPEL_SOURCES)
INTERFACES = $(OBL_INTERFACES) $(ELL_INTERFACES) $(OBLFAST_INTERFACES) $(ZEIPEL_INTERFACES)
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
