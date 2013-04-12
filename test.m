clear;
clc; 
Y=[];
load('/home/lvzongting/workplace/iccv/programme/ScSPM/data/Scene-15/livingroom/image_0004.mat');
Y(:,:,1)=feaSet.feaArr(:,1:1024);
load('/home/lvzongting/workplace/iccv/programme/ScSPM/data/Scene-15/livingroom/image_0014.mat');
Y(:,:,2)=feaSet.feaArr(:,1:1024);


%mex -I/usr/include/eigen3 alm_lrr_l21.cpp


% 
Y0=sortrows(1000*Y(:,:,1)',1:128);
Y0=Y0';
D=sortrows(1000*[Y(:,1:256,1) Y(:,1:256,2)  ]',1:128);
D=D';
tic;
Z = solve_lrr_c(Y0,D,0.15);
toc
tic;
Z_m = solve_lrr_m(Y0,D,0.15);
toc
mesh(Z-Z_m)