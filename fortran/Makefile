#### Define "all"
###all: dais.so simple.so brick_te.so gsic_magicc.so doeclim.so sneasy.so



# Platform detection
ifeq ($(OS),Windows_NT)
    OS_detected := Windows
else
    OS_detected := $(shell uname -s)
endif

# Fortran compiler:
ifeq ($(OS_detected),Windows)
    F90 = gfortran # g95
else ifeq ($(OS_detected),Darwin) # Mac
    F90 = /usr/local/bin/gfortran # g95
else
    F90 = ifort # g95
endif

# Shared library extension
ifeq ($(OS_detected),Windows)
    SHAREDEXT := .dll
else
    SHAREDEXT := .so
endif

# Define "all"
all: dais$(SHAREDEXT) simple$(SHAREDEXT) brick_te$(SHAREDEXT) gsic_magicc$(SHAREDEXT) doeclim$(SHAREDEXT) sneasy$(SHAREDEXT)



#### Fortran compiler:
###F90 = /usr/local/bin/gfortran # g95

# Directories: objpath, srcpath
OP = ./obj
SP = ./src

# Compiler flags:
Optimization = -O3
SOFlags = -fPIC
#Modpath = -module $(OP)
#Flags = ${SOFlags} ${Optimization} -fno-range-check #$(Modpath)
Flags = ${SOFlags} -fno-range-check #$(Modpath)

# Preprocessor options:
CPP = cpp
CPPFLAGS = -DMOC_NOISE

# ==============================================================
## The main programs/subroutines ##
# ==============================================================

dais$(SHAREDEXT): $(OP)/global.o $(OP)/dais.o $(OP)/run_dais.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

simple$(SHAREDEXT): $(OP)/global.o $(OP)/simple.o $(OP)/run_simple.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

brick_te$(SHAREDEXT): $(OP)/global.o $(OP)/brick_te.o $(OP)/run_brick_te.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

gsic_magicc$(SHAREDEXT): $(OP)/global.o $(OP)/gsic_magicc.o $(OP)/run_gsic_magicc.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

doeclim$(SHAREDEXT): $(OP)/global.o $(OP)/doeclim.o $(OP)/run_doeclim.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

sneasy$(SHAREDEXT): $(OP)/global.o $(OP)/sneasy.o $(OP)/CCM.o $(OP)/doeclim.o $(OP)/moc_boxmodel.o $(OP)/rndseed.o $(OP)/mt19937ar.o $(OP)/run_sneasy.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

# ==============================================================
## Main Object files ##
# ==============================================================

## DOECLIM ##
$(OP)/run_doeclim.o: $(OP)/global.o $(SP)/run_doeclim.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_doeclim.f90 -o $@

$(OP)/doeclim.o: $(OP)/global.o $(SP)/doeclim.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/doeclim.f90 -o $@

## DAIS ##
$(OP)/run_dais.o: $(OP)/global.o $(OP)/dais.o $(SP)/run_dais.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_dais.f90 -o $@

$(OP)/dais.o: $(OP)/global.o $(SP)/dais.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/dais.f90 -o $@

## SIMPLE ##
$(OP)/run_simple.o: $(OP)/global.o $(OP)/simple.o $(SP)/run_simple.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_simple.f90 -o $@

$(OP)/simple.o: $(OP)/global.o $(SP)/simple.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/simple.f90 -o $@

## brick_te ##
$(OP)/run_brick_te.o: $(OP)/global.o $(OP)/brick_te.o $(SP)/run_brick_te.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_brick_te.f90 -o $@

$(OP)/brick_te.o: $(OP)/global.o $(SP)/brick_te.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/brick_te.f90 -o $@

## gsic_magicc ##
$(OP)/run_gsic_magicc.o: $(OP)/global.o $(OP)/gsic_magicc.o $(SP)/run_gsic_magicc.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_gsic_magicc.f90 -o $@

$(OP)/gsic_magicc.o: $(OP)/global.o $(SP)/gsic_magicc.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/gsic_magicc.f90 -o $@

## SNEASY ##
$(OP)/run_sneasy.o: $(OP)/global.o $(OP)/sneasy.o $(OP)/doeclim.o $(OP)/CCM.o $(SP)/run_sneasy.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_sneasy.f90 -o $@

$(OP)/sneasy.o: $(OP)/global.o $(OP)/doeclim.o $(OP)/CCM.o $(OP)/moc_boxmodel.o $(OP)/mt19937ar.o $(OP)/rndseed.o $(SP)/sneasy.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/sneasy.f90 -o $@

## CCM ##
$(OP)/run_ccm.o: $(OP)/global.o $(SP)/run_ccm.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_ccm.f90 -o $@

$(OP)/CCM.o: $(OP)/global.o $(SP)/CCM.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/CCM.f90 -o $@

## MOC box model ##
$(OP)/moc_boxmodel.o: $(OP)/global.o $(OP)/rndseed.o $(OP)/moc_boxmodel.F90
	$(F90) -c $(Flags) $(Incl) $(SP)/moc_boxmodel.f90 -o $@

$(OP)/rndseed.o:  $(OP)/global.o $(OP)/mt19937ar.o $(SP)/rndseed.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/rndseed.f90 -o $@

$(OP)/mt19937ar.o:  $(SP)/mt19937ar.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/mt19937ar.f90 -o $@

## Some global variables used in all/most modules ##
$(OP)/global.o: $(SP)/global.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/global.f90 -o $@

## Preprocessing objects ##
$(OP)/moc_boxmodel.F90: $(SP)/moc_boxmodel.f90
	$(CPP) $(CPPFLAGS) $(SP)/moc_boxmodel.f90 > $@

.PHONY: clean

#Cleaning up
clean:
	rm -f $(OP)/*.o $(OP)/*.mod $(OP)/*.F90
