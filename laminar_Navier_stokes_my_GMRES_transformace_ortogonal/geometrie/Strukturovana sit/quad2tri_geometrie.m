function quad2tri_geometrie

load geometrie;
PX = geometrie{1};
PY = geometrie{2};
QP = geometrie{3};
QQ = geometrie{4};
nq = length(QP(:,1));

% tvorba triangulace
TP = zeros(2*nq,3);
TT = zeros(2*nq,3);

for i = 1:nq
        TP(i,:) = [QP(i,1),QP(i,2),QP(i,3)];
        TP(i+nq,:) = [QP(i,1),QP(i,3),QP(i,4)];
        TT(i,:) = [QQ(i,1),QQ(i,2),i+nq-1];
        TT(i+nq,:) = [i-1,QQ(i,3),QQ(i,4)];
end

geometrie{3} = TP;
geometrie{4} = TT;
geometrie{5} = 3*ones(2*nq,1);

save 'geometrie' geometrie;

mesh{1} = [PX,PY];
mesh{2} = TP+1;

save 'mesh' mesh

figure
triplot(geometrie{3}+1,geometrie{1},geometrie{2});
axis equal