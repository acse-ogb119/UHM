info :	
	@echo "make one TEST=test"
	@echo " - where test is one of $(TESTS)"
	@echo "make all"
	@echo "make run"

-include ../../Make.inc

TEST  = ex1
TESTS = ex1


FC_WORK 	= $(FC) $(FFLAGS) $(EXTRA_FFLAGS) \
		  $(FCDEBUG) $(FCWARN) 

LIBS            = $(UHM_LIB) $(MUMPS_LIB) $(METIS_LIB) \
                  $(LINAL_LIB) $(FLA_LIB) $(LAPACK_LIB) \
		  $(BLAS_LIB) $(CUDA_LIB) 


%.o : %.f90
	@echo "Compiling $<"
	@$(FC_WORK) $(FCINCS) -o $@ -c $<

one : $(TEST).o
	@echo "Linking   $<"
	@$(FC_WORK) -o $(TEST) $(TEST).o  $(LIBS)  $(LDFLAGS)

all :
	for e in $(TESTS) ; do \
		make one run TEST=$$e ;\
	done

run : 
	./$(TEST)

clean :
	@/bin/rm -f *.o *~ 
	@/bin/rm -rf ./ooc_dir
	@for t in $(TESTS) ; do rm -f $$t ; done
