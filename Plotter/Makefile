include Makefile.arch


#------------------------------------------------------------------------------
PLOTTERO      = Plotter.$(ObjSuf)
PLOTTERS      = Plotter.$(SrcSuf)
PLOTTER       = Plotter$(ExeSuf)

OBJS          = $(PLOTTERO)
PROGRAMS      = $(PLOTTER)

#------------------------------------------------------------------------------


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

$(PLOTTER):     $(PLOTTERO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

clean:
		@rm -f $(OBJS) core


distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def so_locations
		@rm -rf cxx_repository

.SUFFIXES: .$(SrcSuf)

###

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

