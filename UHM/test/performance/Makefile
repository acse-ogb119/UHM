info :	
	@echo "make one TEST=test"
	@echo " - where test is one of $(TESTS)"
	@echo "make all"
	@echo "make run"

-include ../../Make.inc

TEST  = uhmtest
TESTS = uhmtest mumpstest pardisotest wsmptest


CXX_WORK 	= $(CXX) $(CFLAGS) $(EXTRA_CFLAGS) \
		  $(CCDEBUG) $(CCWARN) 

LIBS            = $(UHM_LINK_LIBS)


%.o : %.cxx
	@echo "Compiling $<"
	$(CXX_WORK) $(CCINCS) -o $@ -c $<

one : $(TEST).o
	@echo "Linking   $<"
	$(CXX_WORK) -o $(TEST) $(TEST).o  $(LIBS)  $(LDFLAGS)

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
