#! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#! HND X
#! HND X   GAP (Gaussian Approximation Potental)
#! HND X   
#! HND X
#! HND X   Portions of GAP were written by Albert Bartok-Partay, Gabor Csanyi, 
#! HND X   Copyright 2006-2021.
#! HND X
#! HND X   Portions of GAP were written by Noam Bernstein as part of
#! HND X   his employment for the U.S. Government, and are not subject
#! HND X   to copyright in the USA.
#! HND X
#! HND X   GAP is published and distributed under the
#! HND X      Academic Software License v1.0 (ASL)
#! HND X
#! HND X   GAP is distributed in the hope that it will be useful for non-commercial
#! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied 
#! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#! HND X   ASL for more details.
#! HND X
#! HND X   You should have received a copy of the ASL along with this program
#! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original licensors,
#! HND X   Gabor Csanyi or Albert Bartok-Partay. The ASL is also published at
#! HND X   http://github.com/gabor1/ASL
#! HND X
#! HND X   When using this software, please cite the following reference:
#! HND X
#! HND X   A. P. Bartok et al Physical Review Letters vol 104 p136403 (2010)
#! HND X
#! HND X   When using the SOAP kernel or its variants, please additionally cite:
#! HND X
#! HND X   A. P. Bartok et al Physical Review B vol 87 p184115 (2013)
#! HND X
#! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

CUSTOM_F95FLAGS += -frealloc-lhs

ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif
include Makefile.inc
include Makefile.rules


ifeq (${HAVE_DESCRIPTORS_NONCOMMERCIAL},1)
  DEFINES += -DDESCRIPTORS_NONCOMMERCIAL
  GAP1_F95_FILES = make_permutations_noncommercial_v2
else
  GAP1_F95_FILES = 
endif

SOAP_TURBO_F95_FILES = soap_turbo_functions soap_turbo_radial soap_turbo_angular soap_turbo_compress soap_turbo
SOAP_TURBO_F95_SOURCES =  ${addsuffix .f90, ${SOAP_TURBO_F95_FILES}}
SOAP_TURBO_F95_OBJS =  ${addsuffix .o, ${SOAP_TURBO_F95_FILES}}

GAP1_F95_FILES += find_water_triplets_noncommercial descriptors gp_predict descriptors_wrapper clustering 
GAP1_F95_SOURCES = ${addsuffix .f95, ${GAP1_F95_FILES}}
GAP1_F95_OBJS = ${addsuffix .o, ${GAP1_F95_FILES}}

GAP2_F95_FILES = gp_fit gap_fit_module 
GAP2_F95_SOURCES = ${addsuffix .f95, ${GAP2_F95_FILES}}
GAP2_F95_OBJS = ${addsuffix .o, ${GAP2_F95_FILES}}

default: ${GAP_LIBFILE}



ifeq (${USE_MAKEDEP},1)
GAP1_F95_FPP_FILES = ${addsuffix .fpp, ${GAP1_F95_FILES}}
GAP2_F95_FPP_FILES = ${addsuffix .fpp, ${GAP2_F95_FILES}}
GAP1.depend: ${GAP1_F95_FPP_FILES}
	${SCRIPT_PATH}/${MAKEDEP} ${MAKEDEP_ARGS} -- ${addprefix ../../src/GAP/,${GAP1_F95_SOURCES}} > GAP1.depend
GAP2.depend: ${GAP2_F95_FPP_FILES} ${GAP1_F95_FPP_FILES}
	${SCRIPT_PATH}/${MAKEDEP} ${MAKEDEP_ARGS} -- ${addprefix ../../src/GAP/,${GAP2_F95_SOURCES}} > GAP2.depend

-include GAP1.depend
-include GAP2.depend
endif


PROGRAMS = gap_fit 

LIBS = -L. -lquiputils -lquip_core -lgap -latoms
ifeq (${HAVE_THIRDPARTY},1)
  LIBS += -lthirdparty
endif
LIBFILES = libatoms.a ${GAP_LIBFILE} libquip_core.a libquiputils.a

.PHONY : clean allclean depend install

Programs: ${PROGRAMS} 
	#cp ${QUIP_ROOT}/src/GAP/teach_sparse .

${PROGRAMS}: % : ${LIBFILES} ${GAP2_F95_OBJS} ${GAPFIT_LIBFILE}  %.o
	$(LINKER) $(LINKFLAGS) -o $@ ${F95OPTS} $@.o ${GAPFIT_LIBFILE} ${LIBS} ${LINKOPTS}



${GAP_LIBFILE}: ${SOAP_TURBO_F95_OBJS}  ${GAP1_F95_OBJS}
ifneq (${LIBTOOL},)
	${LIBTOOL} -o ${GAP_LIBFILE}  ${SOAP_TURBO_F95_OBJS} ${GAP1_F95_OBJS}
else
	${AR} ${AR_ADD} ${GAP_LIBFILE} $?
endif

${GAPFIT_LIBFILE}: ${GAP2_F95_OBJS}
ifneq (${LIBTOOL},)
	${LIBTOOL} -o ${GAPFIT_LIBFILE} ${GAP2_F95_OBJS}
else
	${AR} ${AR_ADD} ${GAPFIT_LIBFILE} $?
endif




install:
	@if [ ! -d ${QUIP_INSTALLDIR} ]; then \
	  echo "make install: QUIP_INSTALLDIR '${QUIP_INSTALLDIR}' doesn't exist or isn't a directory"; \
	  exit 1; \
	else	\
	  for f in ${PROGRAMS} ; do \
	    echo "Copying $$f to ${QUIP_INSTALLDIR}/$${f}${QUIP_MPI_SUFFIX}" ; \
	    cp $$f ${QUIP_INSTALLDIR}/$${f}${QUIP_MPI_SUFFIX} ; \
	  done ;\
	  #cp ${QUIP_ROOT}/src/GAP/teach_sparse ${QUIP_INSTALLDIR}; \
	fi


clean:
	rm -f *.o *.mod *.mod.save ${GAP_LIBFILE} ${GAPFIT_LIBFILE} ${PROGRAMS} GAP1.depend GAP2.depend


