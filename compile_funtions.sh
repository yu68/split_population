#cp ../functions1.py ./functions.pyx
cython functions.pyx


#for the one with numpy
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o functions.so functions.c


