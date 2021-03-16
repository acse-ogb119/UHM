info :
	@echo "make one TEST=test"
	@echo " - where test is one of $(TESTS)"
	@echo "make all"
	@echo "make run"

-include ../../../Make.inc

TEST  = lu_nopiv
TESTS = chol gemm lu_incpiv lu_nopiv lu_piv syrk trsm trmm

DIRS            =

LIBS            = $(LINAL_LINK_LIBS)

CXX_WORK        = $(CXX) $(CFLGAS) $(EXTRA_CFLAGS) \
		  $(CCINCS) $(CCDEBUG) $(CCWARN) $(CCDEFS)

%.o : %.cxx
	@echo "Compiling $<"
	@$(CXX_WORK) -o $@ -c $<

one : $(TEST).o
	@echo "Linking   $<"
	@$(CXX_WORK) -o $(TEST) $(TEST).o $(LIBS) $(LDFLAGS)

all :
	@for e in $(TESTS) ; do make one TEST=$$e ; done

clean :
	@/bin/rm -f *.o *~ *.log
	@for t in $(TESTS) ; do rm -f $$t ; done
	@for d in $(DIRS); do (cd $$d; make clean ); done
