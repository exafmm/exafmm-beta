make libserial_coulombVdW.a
cd ../../../../../tool/gpu
cp fmm/exafmm/wrapper/libserial_coulombVdW.a gcc/libgpufmmcoulomb.a
cd ../..
rm exec/gpu_M/charmm
./install.com gpu xxlarge M mpif90
cd ala3
../exec/gpu_M/charmm -i fmm.inp | tee fmm.out
