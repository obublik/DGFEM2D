function zobraz

load 'W.txt';
load 'geometrie';
PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3}+1;
typ = geometrie{5};
np = length(PX);
nt = length(TP(:,1));

osy = [min(PX) max(PX) min(PY) max(PY)];

I = 1:nt;
Cx = (PX(TP(I,1))+PX(TP(I,2))+PX(TP(I,3)))/3;
Cy = (PY(TP(I,1))+PY(TP(I,2))+PY(TP(I,3)))/3;

figure;
quiver(Cx,Cy,W(I,2)./W(I,1),W(I,3)./W(I,1),0.2);
axis equal;

Wb = prevedNaUzly(TP,W,typ,np);

% tvorba triangulace
nq = length(TP(:,1));
tri = [];
for i = 1:nq
    for j = 2:typ(i)-1
        tri = [tri;TP(i,[1,j,j+1])];
    end
end

kapa = 1.4;
I = 1:np;
R = Wb(I,1);
P = (kapa-1)*(Wb(I,4)-1./(2*Wb(I,1)).*(Wb(I,2).^2 + Wb(I,3).^2));
M =((Wb(I,2).^2 + Wb(I,3).^2)./(kapa*P(I).*Wb(I,1))).^(1/2);
T = P./R;

figure('name','mach');
tricontf(PX,PY,tri,M,30);
axis(osy);
axis equal;
colorbar;

figure('name','tlak');
tricontf(PX,PY,tri,P,30);
axis(osy);
axis equal;
colorbar;

figure('name','hustota');
tricontour(tri,PX,PY,T,30);
axis(osy);
axis equal;
colorbar;

figure('name','teplota');
tricontour(tri,PX,PY,T,30);
axis(osy);
axis equal;
colorbar;

display(['Max Mach: ',num2str(max(M))]);

function Wb = prevedNaUzly(TP,W,typ,np)
    
    nt = length(W(:,4));
    Wb = zeros(np,4);
    poc = zeros(np,1);
    for i = 1:nt
        for j = 1:typ(i)
            Wb(TP(i,j),:) = Wb(TP(i,j),:) + W(i,:);
            poc(TP(i,j)) = poc(TP(i,j)) + 1;
        end
    end
    
    for i = 1:np
        Wb(i,:) = Wb(i,:)/poc(i);
    end
    





