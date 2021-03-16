info :	
	@echo "make one TEST=test"
	@echo " - where test is one of $(TESTS)"
	@echo "make all"
	@echo "make run"

-include ../../Make.inc

TEST  = ex1
TESTS = ex1


CXX_WORK 	= $(CXX) $(CFLAGS) $(EXTRA_CFLAGS) \
		  $(CCDEBUG) $(CCWARN) 

LIBS            = $(UHM_LIB) $(MUMPS_LIB) $(METIS_LIB) \
                  $(LINAL_LIB) $(FLA_LIB) $(LAPACK_LIB) \
		  $(BLAS_LIB) $(CUDA_LIB) 


%.o : %.cxx
	@echo "Compiling $<"
	@$(CXX_WORK) $(CCINCS) -o $@ -c $<

one : $(TEST).o
	@echo "Linking   $<"
	@$(CXX_WORK) -o $(TEST) $(TEST).o  $(LIBS)  $(LDFLAGS)

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
