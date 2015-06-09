function [ out ] = gab_vals( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xm = x(1);
C = x(2);
K = x(3);

a = 1.8959;
b = -2.1943;
c = 0.7969;

out(1) = (1/(xm*C*K)) - a;
out(2) = (C-2)/(xm*C) - b;
out(3) = (K*(1-C))/(xm*C) - c;

end

