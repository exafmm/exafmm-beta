#! /bin/bash

#export LD_LIBRARY_PATH=./:$LD_LIBRARY_PATH

./emutestvg w || exit 1

md5sum r1.g80emu > r1.md5sum
md5sum rsqrt.g80emu > rsqrt.md5sum

diff r1.md5sum r1.md5sum_org
if [ $? == "0" ] 
then
  echo "succeed to generate r1.g80emu";
else
  echo "** error : failed to generate r1.g80emu **";
  /bin/rm -f r1.g80emu
fi

diff rsqrt.md5sum rsqrt.md5sum_org
if [ $? == "0" ] 
then
  echo "succeed to generate rsqrt.g80emu";
else
  echo "** error : failed to generate rsqrt.g80emu **";
  /bin/rm -f rsqrt.g80emu
fi
