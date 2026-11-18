# Build the CCFULL program of this folder:  make
# On Apple silicon use a gfortran that makes arm64 code, here:  make FC=/opt/homebrew/bin/gfortran
FC = gfortran

ccfull: ccfull.f90
	$(FC) -O2 -ffree-line-length-none -o ccfull ccfull.f90
