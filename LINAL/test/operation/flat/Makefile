info :
	@echo "make one TEST=test"
	@echo " - where test is one of $(TESTS)"
	@echo "make all"
	@echo "make run"

-include ../../../Make.inc

TEST  = norm
TESTS = norm

DIRS            =

LIBS            = $(LINAL_LINK_LIBS)

CXX_WORK        = $(CXX) $(CFLGAS) $(EXTRA_CFLAGS) \
		  $(CCINCS) $(CCDEBUG) $(CCWARN) $(CCDEFS)

%.o : %.cxx
	@echo "Compiling $<"
	@$(CXX_WORK) -o $@ -c $<

one : $(TEST).o
	@echo "Linking   $<"
	@$(CXX_WORK) -o $(TEST) $(TEST).o  $(LDFLAGS) $(LIBS)

all :
	@for e in $(TESTS) ; do make one TEST=$$e ; done

run :
	@echo "Run " $(TEST)
	@./$(TEST)

clean :
	@/bin/rm -f *.o *~ *.log
	@for t in $(TESTS) ; do rm -f $$t ; done
	@for d in $(DIRS); do (cd $$d; make clean ); done
