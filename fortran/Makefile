# Define "all"
all: dais.so simple.so brick_te.so gsic_magicc.so doeclim.so dais_fastdyn.so

# Fortran compiler:
#F90 = /usr/local/bin/gfortran # g95
F90 = /usr/bin/gfortran # g95

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
#CPP = cpp
#CPPFLAGS = -DMOC_NOISE

# ==============================================================
## The main programs/subroutines ##
# ==============================================================

dais.so: $(OP)/global.o $(OP)/dais.o $(OP)/run_dais.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

simple.so: $(OP)/global.o $(OP)/simple.o $(OP)/run_simple.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

brick_te.so: $(OP)/global.o $(OP)/brick_te.o $(OP)/run_brick_te.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

gsic_magicc.so: $(OP)/global.o $(OP)/gsic_magicc.o $(OP)/run_gsic_magicc.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

doeclim.so: $(OP)/global.o $(OP)/doeclim.o $(OP)/run_doeclim.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

dais_fastdyn.so: $(OP)/global.o $(OP)/dais_fastdyn.o $(OP)/run_dais_fastdyn.o
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

## DAIS, with fast dynamics emulator ##
$(OP)/run_dais_fastdyn.o: $(OP)/global.o $(OP)/dais.o $(SP)/run_dais_fastdyn.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_dais_fastdyn.f90 -o $@

$(OP)/dais_fastdyn.o: $(OP)/global.o $(SP)/dais_fastdyn.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/dais_fastdyn.f90 -o $@

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

## Some global variables used in all/most modules ##
$(OP)/global.o: $(SP)/global.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/global.f90 -o $@

.PHONY: clean

#Cleaning up
clean:
	rm -f $(OP)/*.o $(OP)/*.mod $(OP)/*.F90
