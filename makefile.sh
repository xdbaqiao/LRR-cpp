 g++ `pkg-config --cflags --libs eigen3` -o lrr LowRankRepresentation.cpp

#gcc with mex
#g++ -fPIC -DMATLAB_MEX_FILE `pkg-config --cflags --libs eigen3` -I/home/lvzongting/opt/matlab/extern/include -c alm_lrr_l21.cpp

#mex with gcc
#mex -I/usr/include/eigen3 alm_lrr_l21.cpp
