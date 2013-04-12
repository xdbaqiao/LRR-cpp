 g++ -pg `pkg-config --cflags --libs eigen3` -o lrr lrr.cpp
#./lrr
#gprof lrr|gprof2dot |dot -Tpng -o result.png

#生成统计图表


#gcc with mex
#g++ -fPIC -DMATLAB_MEX_FILE `pkg-config --cflags --libs eigen3` -I/home/lvzongting/opt/matlab/extern/include -c lrralm.c

#mex with gcc
#mex -I/usr/include/eigen3 lrralm.cpp
