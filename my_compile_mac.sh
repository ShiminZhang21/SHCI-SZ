

make -f Makefile.mac -j all 2>&1 | tee makelog_mac.txt


#move executable to build folder:
mv shci build
