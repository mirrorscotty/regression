Data = csvread('gab_data_andrieu.csv', 1, 0);
T = Data(:,1);
m = Data(:,4);
C = Data(:,2);
k = Data(:,3);

T = ones(length(T),1)*273.15 + T;

bm = regress(log(m), [ones(length(T),1) 1./T]);
bc = regress(log(C), [ones(length(T),1) 1./T]);
bk = regress(log(k), [ones(length(T),1) 1./T]);

bm(1) = exp(bm(1));
bc(1) = exp(bc(1));
bk(1) = exp(bk(1));

figure(1);
hold on;
plot(T, m, 'o');
plot(T, bm(1)*exp(bm(2)./T));
hold off;

figure(2);
hold on;
plot(T, C, 'o');
plot(T, bc(1)*exp(bc(2)./T));
hold off;

figure(3);
hold on;
plot(T,k, 'o');
plot(T, bk(1)*exp(bk(2)./T));
hold off;

bc
bk
bm
