function [ J ] = CreepModelv2(beta, t)
J0 = 0.00000001097289;
J1 = 0.00000001388361 - J0;
J = J0 + J1 * (1-exp(-t/beta));
