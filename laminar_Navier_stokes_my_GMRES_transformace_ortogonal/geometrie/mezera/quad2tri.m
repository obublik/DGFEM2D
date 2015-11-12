function quad2tri

load 'geometrie1'

PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3};
TTo = geometrie{4};

% tvorba triangulace
tri = [];
for i = 1:length(TP(:,1));
        tri = [tri;[TP(i,1),TP(i,2),TP(i,3)]];
        tri = [tri;[TP(i,1),TP(i,3),TP(i,4)]];
end

nt = length(TTo(:,1));
TT = zeros(2*nt,3);
for i = 1:nt
    TT(i,:) = [TTo(i,1),TTo(i,2),nt+i-1];
    TT(nt+i,:) = [i-1, TTo(i,3),TTo(i,4)];
end

geometrie{3} = tri;
geometrie{4} = TT;
geometrie{5} = 3*ones(2*nt,1);

figure
triplot(tri+1,PX,PY)
hold on;

save 'geometrie' geometrie;