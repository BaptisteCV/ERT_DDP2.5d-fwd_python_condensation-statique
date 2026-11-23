################################################################################
# Automatically-generated file. Do not edit!
################################################################################

USER_OBJS :=

#LIBS := -lCGAL -lumfpack -lcholmod -lamd -lcolamd -lcblas -lGL -lglut
#LIBS := -lumfpack -lcholmod -lamd -lcolamd -lGL -lglut  #uncomment this line to generate a DFN
LIBS = -L$(SUITESPARSE_ROOT)/lib -L$(METIS_ROOT)/lib -lumfpack -lcholmod -lamd -lcolamd -Wl,--no-as-needed -lmetis -Wl,--as-needed
