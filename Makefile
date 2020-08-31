include Make.inc


all:  library 

library: libdir mlp
#cbnd

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
	(if test ! -d include ; then mkdir include; fi)
	(if test ! -d modules ; then mkdir modules; fi;)	
	($(INSTALL_DATA) Make.inc  include/Make.inc.mld2p4)
        

mlp:
	$(MAKE) -C mlprec all
cbnd: mlp
	$(MAKE) -C cbind all
install: all
	mkdir -p $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR)
	mkdir -p  $(INSTALL_INCLUDEDIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_INCLUDEDIR)/Make.inc.mld2p4
	mkdir -p $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) include/*.h $(INSTALL_INCLUDEDIR)
	mkdir -p $(INSTALL_MODULESDIR) && \
	   $(INSTALL_DATA) modules/*$(.mod) $(INSTALL_MODULESDIR)
	mkdir -p  $(INSTALL_DOCSDIR) && \
	   /bin/cp -fr docs/*pdf docs/html $(INSTALL_DOCSDIR)
	mkdir -p  $(INSTALL_DOCSDIR) && \
	   $(INSTALL_DATA) README.md LICENSE $(INSTALL_DOCSDIR)
	mkdir -p  $(INSTALL_SAMPLESDIR) && \
		 mkdir -p  $(INSTALL_SAMPLESDIR)/simple &&\
	 	 mkdir -p  $(INSTALL_SAMPLESDIR)/advanced && \
		(cd examples; /bin/cp -fr pdegen fileread $(INSTALL_SAMPLESDIR)/simple ) && \
		(cd tests; /bin/cp -fr pdegen fileread $(INSTALL_SAMPLESDIR)/advanced )
cleanlib:
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh))
	(cd include; /bin/rm -f *.a *$(.mod) *$(.fh))
	(cd modules; /bin/rm -f *.a *$(.mod) *$(.fh))

veryclean: cleanlib
	(cd mlprec; make veryclean)
	(cd examples/fileread; make clean)
	(cd examples/pdegen; make clean)
	(cd tests/fileread; make clean)
	(cd tests/pdegen; make clean)

check: all
	make check -C tests/pdegen

clean:
	(cd mlprec; make clean)
