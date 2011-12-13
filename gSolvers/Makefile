
files = gPCG.cu gPBiCG.cu

default:	build_and_install

libgsolvers.so: $(files)
	nvcc -arch=sm_13 --compiler-options -O2,-Wall,-fPIC,-DWM_${WM_PRECISION_OPTION} --shared --linker-options -soname,libgsolvers.so $^ -o $@

build_and_install: libgsolvers.so
	cp libgsolvers.so ${FOAM_USER_LIBBIN}
	wmake libso

clean:
	@-rm *.so *.o
	@-wclean
	@-rm -rf lnInclude
