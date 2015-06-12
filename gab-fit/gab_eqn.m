function [ Xdb ] = gab_eqn( beta,x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

aw = x(:,1);
T = x(:,2);
T = ones(length(T),1) + T;

for i = 1:length(aw)

C = beta(1)*exp(beta(2)/T(i));
K = beta(3)*exp(beta(4)/T(i));
Xm = beta(5)*exp(beta(6)/T(i));

Xdb(i) = C*K*Xm*aw(i)/((1-K*aw(i)).*(1-K*aw(i)+C*K*aw(i)));
Xdb = Xdb';
end

end

