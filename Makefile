include Make.inc


all:  objs lib

objs: amgp cbnd

lib: libdir objs
	cd amgprec && $(MAKE) lib
	cd cbind && $(MAKE)  lib

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
	(if test ! -d include ; then mkdir include; fi)
	(if test ! -d modules ; then mkdir modules; fi;)	
	($(INSTALL_DATA) Make.inc  include/Make.inc.amg4psblas)
        

amgp:
	cd amgprec && $(MAKE) objs
cbnd: amgp
	cd cbind && $(MAKE)  objs

install: lib
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
	(cd amgprec && $(MAKE) veryclean)
	(cd samples/simple/fileread && $(MAKE) clean)
	(cd samples/simple/pdegen && $(MAKE) clean)
	(cd samples/advanced/fileread && $(MAKE) clean)
	(cd samples/advanced/pdegen && $(MAKE) clean)

check: all
	make check -C samples/advanced/pdegen

clean:
	(cd amgprec && $(MAKE) clean)
