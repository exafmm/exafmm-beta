make libfmmcoulomb.a
cd ../../../../..
cd tool/gpu
cp fmm/exafmm/wrapper/libfmmcoulomb.a gcc
cd ../..
rm exec/gpu/charmm
./install.com gpu xxlarge
cd test/cbenchtest/mbco
../../../exec/gpu/charmm -i fmm.inp | tee fmm.out
