clear;
RawData = csvread('Andrieu.csv',1,0);

for i = 1:size(RawData, 1)
    for j = 1:size(RawData, 2)
        if RawData(i,j) == 0
            RawData(i,j) = NaN;
        end
    end
end

aw = RawData(:,1);
for i = 2:size(RawData, 2);
    Xdb(:,i-1) = RawData(:,i);
end

for i = 1:size(Xdb, 2);
    y = Xdb(:,i)./aw;
    X = [ones(length(aw), 1) aw aw.*aw];

    [b, bint, r, rint, stats] = regress(y,X)
    out(i,:) = [i*10+30 b(3) b(2) b(1)];
end

