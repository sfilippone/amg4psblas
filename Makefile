include Make.inc


all:  library 

library: libdir amgp cbnd
#cbnd

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
	(if test ! -d include ; then mkdir include; fi)
	(if test ! -d modules ; then mkdir modules; fi;)	
	($(INSTALL_DATA) Make.inc  include/Make.inc.amg4psblas)
        

amgp:
	$(MAKE) -C amgprec all
cbnd: amgp
	$(MAKE) -C cbind all
install: all
	mkdir -p $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR)
	mkdir -p  $(INSTALL_INCLUDEDIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_INCLUDEDIR)/Make.inc.amg4psblas
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
		(cd samples/simple; /bin/cp -fr pdegen fileread $(INSTALL_SAMPLESDIR)/simple ) && \
		(cd samples/advanced; /bin/cp -fr pdegen fileread $(INSTALL_SAMPLESDIR)/advanced )
cleanlib:
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh))
	(cd include; /bin/rm -f *.a *$(.mod) *$(.fh))
	(cd modules; /bin/rm -f *.a *$(.mod) *$(.fh))

veryclean: cleanlib
	(cd amgprec; make veryclean)
	(cd samples/simple/fileread; make clean)
	(cd samples/simple/pdegen; make clean)
	(cd samples/advanced/fileread; make clean)
	(cd samples/advanced/pdegen; make clean)

check: all
	make check -C samples/advanced/pdegen

clean:
	(cd amgprec; make clean)
