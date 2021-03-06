# ----------------------------------------------------------------------------
#
# Makefile for DHELAS 3.0 library
# Feb. 28, 2001
#
# ----------------------------------------------------------------------------
#
# Use a TAB to precede shell commands (e.g., f90, ar, etc).
# Note: this Makefile uses features that *may not*
# be supported by make utilities other than GNU make.
#
# ----------------------------------------------------------------------------

FC            = gfortran

FFLAGS        = -w 

LD	      = ld

DEST	      = ~/lib

EXTHDRS	      =

HDRS	      =

INSTALL	      = /etc/install

LIBRARY	      = libdhelas3.gf.a

MAKEFILE      = Makefile

OBJS	      = momntx.o mom2cx.o boostx.o rotxxx.o \
		coupsm.o hdecay.o haber.o  feynhiggs.o \
		ixxxxx.o oxxxxx.o vxxxxx.o sxxxxx.o \
		iovxxx.o fvixxx.o fvoxxx.o jioxxx.o j3xxxx.o \
		iosxxx.o fsixxx.o fsoxxx.o hioxxx.o \
		vvvxxx.o jvvxxx.o gggxxx.o jggxxx.o \
		vvsxxx.o jvsxxx.o hvvxxx.o hvvxxu.o \
		vssxxx.o jssxxx.o hvsxxx.o \
		sssxxx.o hssxxx.o \
		wwwwxx.o jwwwxx.o w3w3xx.o jw3wxx.o \
		ggggxx.o jgggxx.o \
		vvssxx.o jvssxx.o hvvsxx.o \
		ssssxx.o hsssxx.o \
		eaixxx.o eaoxxx.o jeexxx.o \
		ioscxx.o fsicxx.o fsocxx.o hiocxx.o \
		iovcxx.o fvicxx.o fvocxx.o jiocxx.o \
		iovdmx.o fvidmx.o fvodmx.o jiodmx.o \
		iosgld.o fsigld.o fsogld.o hiogld.o \
		iovgld.o fvigld.o fvogld.o jiogld.o \
		txxxxx.o txxxx2.o \
		iotxkk.o iovtkk.o vvtxkk.o ftixkk.o ftoxkk.o

PRINT	      = pr

SHELL	      = /bin/sh

SYSHDRS	      =

MFLAGS        = -e

.F.o:
	$(FC) $(FFLAGS) -c $<

all:		$(LIBRARY)

helas:
		-rm -f *.o
		export FFLAGS='-O +cpp' && \
		export LIBRARY=libdhelas3.a && \
		$(MAKE) $(MFLAGS)

helas_check:
		-rm -f *.o
		export FFLAGS='-O +cpp -DHELAS_CHECK' && \
		export LIBRARY=libdhelas3_check.a && \
		$(MAKE) $(MFLAGS)

install-helas:
		export LIBRARY=libdhelas3.a && \
		$(MAKE) $(MFLAGS) install	

install-helas_check:
		export LIBRARY=libdhelas3.a && \
		$(MAKE) $(MFLAGS) install	

$(LIBRARY):	$(OBJS)
		@echo  "Loading $(LIBRARY) ... "
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags

install:	$(LIBRARY)
	        @echo Installing $(LIBRARY) in $(DEST)
	        @if [ $(DEST) != . ]; then \
	        (rm -f $(DEST)/$(LIBRARY); $(INSTALL) -f $(DEST) $(LIBRARY)); fi

$(DEST)/$(LIBRARY): $(SRCS) $(HDRS) $(EXTHDRS)
	        @$(MAKE) -f $(MAKEFILE) ROOT=$(ROOT) DEST=$(DEST) install
