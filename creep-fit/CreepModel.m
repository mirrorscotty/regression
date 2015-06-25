function [ J ] = CreepModel(beta, t)
J=0;
J0 = beta(1);
n = (length(beta)-1)/2;

for i=1:n
    Ji(i, 1) = beta(2*i);
    taui(i, 1) = beta(2*i+1);
end

J = J0;

for i=1:n
    %Jval = Ji(i)*Ji(i);
    %tauval = taui(i)*taui(i);
    %J = J + Jval * (1-exp(-t/tauval));
    J = J + J(i) * (1-exp(-t/taui(i)));
end
%J=(1-J)*J0*J0;