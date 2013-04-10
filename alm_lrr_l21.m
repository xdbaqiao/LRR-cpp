function [Z,E] = alm_lrr_l21(X,A,lambda)
% This routine solves the following nuclear-norm optimization problem,
% min |Z|_*+lambda*|E|_2,1
% s.t., X = AZ+E
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
%        A -- D*M matrix of a dictionary, M is the size of the dictionary
