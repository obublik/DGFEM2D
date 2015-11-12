function zobraz_pom

S = load('pom.txt');
load 'geometrie';
PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3}+1;
typ = geometrie{5};
np = length(PX);
nt = length(TP(:,1));

osy = [min(PX) max(PX) min(PY) max(PY)];

Sb = prevedNaUzly(TP,S,typ,np);

% tvorba triangulace
nq = length(TP(:,1));
tri = [];
for i = 1:nq
    for j = 2:typ(i)-1
        tri = [tri;TP(i,[1,j,j+1])];
    end
end

figure('name','pom');
tricontour(tri,PX,PY,Sb,30);
axis(osy);
axis equal;
colorbar;

function Sb = prevedNaUzly(TP,S,typ,np)
    
    nt = length(S);
    Sb = zeros(np,1);
    poc = zeros(np,1);
    for i = 1:nt
        for j = 1:typ(i)
            Sb(TP(i,j)) = Sb(TP(i,j)) + S(i);
            poc(TP(i,j)) = poc(TP(i,j)) + 1;
        end
    end
    
    for i = 1:np
        Sb(i) = Sb(i)/poc(i);
    end