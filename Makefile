#!/usr/bin/make

#main building variables
DSRC    = .
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = gfortran-8
OPTSC   = -c -std=f2008 -ffree-line-length-512 -J mod
OPTSL   =  -J mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)MAIN: $(MKDIRS) $(DOBJ)main.o \
	$(DOBJ)labor2.o \
	$(DOBJ)partest.o \
	$(DOBJ)solveactivelife.o \
	$(DOBJ)solveinretirement.o \
	$(DOBJ)labor1.o \
	$(DOBJ)spline.o \
	$(DOBJ)simulation.o \
	$(DOBJ)labor3.o \
	$(DOBJ)labors.o \
	$(DOBJ)solvefirstactive.o \
	$(DOBJ)statistics.o \
	$(DOBJ)solvelastactive.o \
	$(DOBJ)lsupply.o \
	$(DOBJ)qrfac.o \
	$(DOBJ)hybrd.o \
	$(DOBJ)qform.o \
	$(DOBJ)r1mpyq.o \
	$(DOBJ)fdjac1.o \
	$(DOBJ)dogleg.o \
	$(DOBJ)dpmpar.o \
	$(DOBJ)r1updt.o \
	$(DOBJ)enorm.o
	@rm -f $(filter-out $(DOBJ)main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MAIN

#compiling rules
$(DOBJ)myinterpolation.o: ./src/MyInterpolation.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)labor2.o: ./src/labor2.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)partest.o: ./src/partest.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)myinterpolation.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)main.o: ./src/main.f90 \
	$(DOBJ)utilities.o \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)tauchen.o \
	$(DOBJ)hybrd_wrapper.o \
	$(DOBJ)solveretired_egm.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solveactivelife.o: ./src/SolveActiveLife.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solveinretirement.o: ./src/SolveInRetirement.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o \
	$(DOBJ)myinterpolation.o \
	$(DOBJ)glob0.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)labor1.o: ./src/labor1.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)pyplot_module.o: ./src/pyplot_module.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)spline.o: ./src/spline.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)simulation.o: ./src/Simulation.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)labor3.o: ./src/labor3.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)labors.o: ./src/labors.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)glob0.o: ./src/glob0.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)model_parameters.o: ./src/Model_Parameters.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solvefirstactive.o: ./src/Solvefirstactive.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)statistics.o: ./src/Statistics.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solvelastactive.o: ./src/Solvelastactive.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solveretired_egm.o: ./src/SolveRetired_EGM.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)utilities.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)pyplot_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)pppack.o: ./src/pppack.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)policyfunctions.o: ./src/PolicyFunctions.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)pppack.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)lsupply.o: ./src/lsupply.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)hybrd_wrapper.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)utilities.o: ./src/Utilities.f90 \
	$(DOBJ)model_parameters.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)tauchen.o: ./src/tauchen.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)qrfac.o: ./src/minpack-hybrd/qrfac.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hybrd.o: ./src/minpack-hybrd/hybrd.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)qform.o: ./src/minpack-hybrd/qform.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)r1mpyq.o: ./src/minpack-hybrd/r1mpyq.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)fdjac1.o: ./src/minpack-hybrd/fdjac1.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)dogleg.o: ./src/minpack-hybrd/dogleg.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)dpmpar.o: ./src/minpack-hybrd/dpmpar.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)r1updt.o: ./src/minpack-hybrd/r1updt.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hybrd_wrapper.o: ./src/minpack-hybrd/hybrd_wrapper.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)enorm.o: ./src/minpack-hybrd/enorm.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
