#
# LINAL AND MUMPS
#
REV = 1.0
J ?= 1

MAKE_INC     = ./MAKE_INC/dev/kyungjoo
INSTALL_PATH = ./$(REV)
TARBALL      = ./uhm.$(REV).tar.gz

.PHONY: all

slink:
	@echo "Create link file Make.inc"
	@(cd LINAL; rm -f Make.inc; ln -s ../$(MAKE_INC) Make.inc;)
	@(cd UHM; rm -f Make.inc; ln -s ../$(MAKE_INC) Make.inc;)

all: slink
	@echo "Building LINAL AND UHM"
	(cd LINAL; make lib -j$(J) )
	(cd UHM; make lib -j$(J) )

install: all
	@echo "Install Libraries" $(INSTALL_PATH)
	@mkdir -p $(INSTALL_PATH)/include
	cp -rf ./LINAL/include/* $(INSTALL_PATH)/include
	cp -rf ./UHM/include/* $(INSTALL_PATH)/include
	@mkdir -p $(INSTALL_PATH)/lib
	cp -rf ./LINAL/lib/*.a $(INSTALL_PATH)/lib
	cp -rf ./UHM/lib/*.a $(INSTALL_PATH)/lib

tar:
	@echo "Creating tar ball"
	tar -zcvf $(TARBALL) \
	--exclude='./UHM/uhmfile' \
	--exclude='./*/docs' \
	--exclude='*~' \
	--exclude='\#*' \
	./LINAL ./UHM ./Makefile ./MAKE_INC ./doxy.conf

clean:
	@(echo "Clean LINAL")
	@(cd LINAL; make clean)
	@(echo "Clean UHM")
	@(cd UHM; make clean)
	@(echo "Clean misc")
	@(find . -iname "*~" -exec rm {} \;)
	@(find . -iname "\#*" -exec rm {} \;)
