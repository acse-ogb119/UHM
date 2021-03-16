#-------------------------------------------------
# - Library path -
-include ./Make.inc

SRC_PATH = ./src
INC_PATH = ./include
OBJ_PATH = ./_obj_

# compiler
CCDEFS          =
CCINCS 		= -I$(SRC_PATH) \
	  	  -I$(INC_PATH) \
	  	  -I$(FLA_INC_PATH) 

# cuBLAS addition
ifdef CUDA_PATH
CCDEFS		= -D LINAL_GPU_ENABLE
CCINCS         += -I$(CUDA_INC_PATH)
endif


CXX_WORK 	= $(CXX) $(CFLAGS) $(EXTRA_CFLAGS) \
		  $(CCINCS) $(CCDEBUG) $(CCWARN) $(CCDEFS)

# includes
INC_FILE	= linal.hxx \
		  linal/common.hxx \
		  linal/const.hxx \
		  linal/flat.hxx \
		  linal/gpu.hxx \
		  linal/hier.hxx \
		  linal/matrix.hxx \
		  linal/operation.hxx \
		  linal/util.hxx 

# sources
SRC_FILE	= dense/chol.cxx \
		  dense/gemm.cxx \
		  dense/graphviz.cxx \
		  dense/lu_incpiv.cxx \
		  dense/lu_nopiv.cxx \
		  dense/lu_piv.cxx \
		  dense/pivot.cxx \
		  dense/scal.cxx \
		  dense/syrk.cxx \
		  dense/trmm.cxx \
		  dense/trsm.cxx \
		  dense/q.cxx \
		  dense/q_inc_var1.cxx \
		  dense/q_var2.cxx \
		  dense/qr.cxx \
		  dense/qr_inc_var1.cxx \
		  dense/qr_var1.cxx \
		  dense/qr_var2.cxx \
		  flat/norm.cxx \
		  gpu/global.cxx \
		  internal/gemm.cxx \
		  internal/trmm.cxx \
		  internal/trsm.cxx \
		  util.cxx

#
INC		= $(addprefix $(INC_PATH)/,$(INC_FILE))
SRC 		= $(addprefix $(SRC_PATH)/,$(SRC_FILE))
OBJ 		= $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.cxx=.o))

# make
.PHONY : lib

lib : liblinal.a
	@echo "LINAL Library is created at " $(LINAL_LIB_PATH)/$<

liblinal.a : all
	@echo "Creating directory " $(LINAL_LIB_PATH)
	@mkdir -p $(LINAL_LIB_PATH)
	@ar cr $(LINAL_LIB_PATH)/$@ $(OBJ)
	@ranlib $(LINAL_LIB_PATH)/$@
all : $(OBJ_PATH)/.dummy $(OBJ)

$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH) 
	@if [ -d $(OBJ_PATH) ]; then \
		touch $@; \
	else mkdir $(OBJ_PATH); touch $@; \
	fi

$(OBJ_PATH)/%.o : $(SRC_PATH)/%.cxx $(INC)
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	@$(CXX_WORK) -o $@ -c $<

.PHONY : clean

clean:
	@(echo "Clean tests")
	@(cd ./test; make clean)
	@(echo "Clean object files and library")
	@(rm -rf $(OBJ_PATH))
	@(rm -f $(LINAL_LIB_PATH)/*.a)
	@(echo "Clean misc")
	@(rm -rf ./docs)
	@(find . -iname "*~" -exec rm {} \;)
