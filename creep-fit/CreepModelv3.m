function [ J ] = CreepModelv3(beta, t)
J0 = 1.0767e-07;
%J1 = 0.00000001388361 - J0;
J = J0 + beta(1) * (1-exp(-t/beta(2))) + beta(3) * (1-exp(-t/beta(4)));
