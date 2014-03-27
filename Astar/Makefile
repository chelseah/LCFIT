ELL_INTERFACES = elliptic.i
ELL_SOURCES = elliptic.cc elliptic.h
PROJECT = testquad
OBJECT = testquad.o eggbox.o poly.o circle.o sphere.o quad.o quad2d.o quadpolar.o quad3d.o quadsphere.o Zeipelmodel.o Zeipel.o Laramodel.o Lara.o Laramodelfilter.o Larafilter.o 
Poly_SOURCES = testquad.cc poly.cc poly.h model.h circle.cc circle.h sphere.cc sphere.h quad.cc quad.h quad2d.cc quad2d.h quadpolar.cc quadpolar.h quad3d.cc quad3d.h quadsphere.cc quadsphere.h integrator.h 
ZEIPEL_SOURCES = quadpolar.cc quad2d.cc quadsphere.cc Zeipelmodel.cc Zeipel.cc Zeipelmodel.h Zeipel.h model.h quadpolar.h quad2d.h quadsphere.h integrator.h Laramodel.cc Laramodel.h Lara.cc Lara.h Laramodelfilter.cc Laramodelfilter.h Larafilter.cc Larafilter.h 
ZEIPEL_INTERFACES =  quadpolar.i quad2d.i quadsphere.i Zeipelmodel.i Zeipel.i Laramodel.i Lara.i Laramodelfilter.i Larafilter.i 
SOURCES = $(ELL_SOURCES) $(ZEIPEL_SOURCES) 
INTERFACES = $(ELL_INTERFACES) $(ZEIPEL_INTERFACES) 
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
	$(CXX) -Wall -pg -o $@ $^

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
