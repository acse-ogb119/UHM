#-------------------------------------------------
# - Library path -
-include ./Make.inc

SRC_PATH = ./src
INC_PATH = ./include
OBJ_PATH = ./_obj_

# environements - already defined in Make.inc
UHM_INTERF_MUMPS_ENABLE   ?= 0
UHM_INTERF_PARDISO_ENABLE ?= 0
UHM_INTERF_WSMP_ENABLE    ?= 0

CCDEFS		= -D UHM_HIER_MATRIX_ENABLE \
		  -D UHM_MULTITHREADING_ENABLE \
		  -D UHM_QR_INC_ENABLE 

CCINCS 		= -I$(SRC_PATH) \
	  	  -I$(INC_PATH) \
		  -I$(LINAL_INC_PATH) \
	  	  -I$(FLA_INC_PATH) \
		  -I$(METIS_INC_PATH)

# includes
INC_FILE	= uhm/common.hxx \
		  uhm/const.hxx \
		  uhm/interf/mumps.hxx \
		  uhm/interf/pardiso.hxx \
		  uhm/interf/sparse.hxx \
		  uhm/interf/wsmp.hxx \
		  uhm/matrix/uhm/fla.hxx \
		  uhm/matrix/uhm/helper.hxx \
		  uhm/matrix/uhm/matrix.hxx \
		  uhm/mesh/element.hxx \
		  uhm/mesh/mesh.hxx \
		  uhm/mesh/node.hxx \
		  uhm/object.hxx \
		  uhm/operation/element.hxx \
		  uhm/operation/mesh.hxx \
		  uhm/operation/scheduler.hxx \
		  uhm/util.hxx \
		  uhm/wrapper/fort.hxx \
		  uhm/wrapper/pardiso_fort.hxx \
		  uhm/wrapper/wsmp_fort.hxx \
		  uhm.hxx

# sources
SRC_FILE 	= interf/sparse.cxx \
		  matrix/uhm/common.cxx \
		  matrix/uhm/fla/common.cxx \
		  matrix/uhm/fla/chol/check.cxx \
		  matrix/uhm/fla/chol/decompose.cxx \
		  matrix/uhm/fla/chol/solve.cxx \
		  matrix/uhm/fla/lu_nopiv/check.cxx \
	 	  matrix/uhm/fla/lu_nopiv/decompose.cxx \
		  matrix/uhm/fla/lu_nopiv/solve.cxx \
		  matrix/uhm/fla/lu_piv/check.cxx \
		  matrix/uhm/fla/lu_piv/decompose.cxx \
		  matrix/uhm/fla/lu_piv/solve.cxx \
		  matrix/uhm/fla/matrix.cxx \
		  matrix/uhm/fla/qr/check.cxx \
		  matrix/uhm/fla/qr/decompose.cxx \
		  matrix/uhm/fla/qr/solve.cxx \
		  matrix/uhm/matrix.cxx \
		  mesh/chol.cxx \
		  mesh/element.cxx \
		  mesh/graphviz.cxx \
		  mesh/interf.cxx \
		  mesh/lu_nopiv.cxx \
		  mesh/lu_piv.cxx \
		  mesh/matrix.cxx \
		  mesh/mesh.cxx \
		  mesh/node.cxx \
		  mesh/qr.cxx \
		  operation/build_tree.cxx \
		  operation/build_tree_var1.cxx \
		  operation/build_tree_var2.cxx \
		  operation/build_tree_var3.cxx \
		  operation/build_tree_var4.cxx \
		  operation/element.cxx \
		  operation/graph.cxx \
		  operation/scheduler.cxx \
		  util.cxx \
		  wrapper/fort.cxx 

## pardiso and wsmp does not require structure
## it is okay to get symbolic interfaces
SRC_FILE        += interf/pardiso/pardiso.cxx \
		   wrapper/pardiso_fort.cxx
SRC_FILE        += interf/wsmp/wsmp.cxx \
		   wrapper/wsmp_fort.cxx


## interface to other solvers
ifeq ($(UHM_INTERF_MUMPS_ENABLE),1)
CCDEFS		+= -D UHM_INTERF_MUMPS_ENABLE 
CCINCS          += -I$(MUMPS_INC_PATH) 

SRC_FILE        += interf/mumps/mumps.cxx \
		   wrapper/mumps_fort.cxx 	
endif

##
ifeq ($(UHM_INTERF_PARDISO_ENABLE),1)
CCDEFS		+= -D UHM_INTERF_PARDISO_ENABLE 
endif

##
ifeq ($(UHM_INTERF_WSMP_ENABLE),1)
CCDEFS		+= -D UHM_INTERF_WSMP_ENABLE 
endif

##
CXX_WORK 	= $(CXX) $(CFLAGS) $(EXTRA_CFLAGS) \
		  $(CCINCS) $(CCDEBUG) $(CCWARN) $(CCDEFS)

#
INC		= $(addprefix $(INC_PATH)/,$(INC_FILE))
SRC 		= $(addprefix $(SRC_PATH)/,$(SRC_FILE))
OBJ 		= $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.cxx=.o))

# make
.PHONY : lib

lib : libuhm.a
	@echo "UHM Library is created at " $(UHM_LIB_PATH)/$<

libuhm.a : all
	@echo "Creating directory " $(UHM_LIB_PATH)
	@mkdir -p $(UHM_LIB_PATH)
	@ar cr $(UHM_LIB_PATH)/$@ $(OBJ)
	@ranlib $(UHM_LIB_PATH)/$@

all : $(OBJ_PATH)/.dummy $(OBJ)

$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH) 
	@if [ -d $(OBJ_PATH) ]; then \
		touch $@; \
	else mkdir $(OBJ_PATH); touch $@; \
	fi

$(OBJ_PATH)/%.o : $(SRC_PATH)/%.cxx
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	@$(CXX_WORK) -o $@ -c $<

.PHONY : clean

clean:
	@(echo "Clean tests")
	@(cd ./test; make clean)
	@(echo "Clean object files and library")
	@(rm -rf $(OBJ_PATH))
	@(rm -f $(UHM_LIB_PATH)/*.a)
	@(echo "Clean misc")
	@(rm -rf ./docs)
	@(find . -iname "*~" -exec rm {} \;)
