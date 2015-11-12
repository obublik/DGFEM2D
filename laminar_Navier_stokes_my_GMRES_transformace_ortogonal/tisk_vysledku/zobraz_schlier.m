function zobraz_schlier
global PX PY tri
[L_ref, u_ref, ro_ref, p_ref, eta_ref, Re, Pr, kapa] = referencni_hodnoty;

load 'W.txt';
load 'geometrie';
PX = L_ref*(geometrie{1});
PY = L_ref*(geometrie{2}-100.5);
TP = geometrie{3}+1;
TT = geometrie{4};
typ = geometrie{5};
np = length(PX);
nt = length(TP(:,1));

% I = 1:nt;
% Cx = (PX(TP(I,1))+PX(TP(I,2))+PX(TP(I,3)))/3;
% Cy = (PY(TP(I,1))+PY(TP(I,2))+PY(TP(I,3)))/3;
% 
% figure;
% quiver(Cx,Cy,W(I,2)./W(I,1),W(I,3)./W(I,1),0.2);
% axis equal;

R = prevedNaUzly(TP,W,typ,np);

% tvorba triangulace
nq = length(TP(:,1));
tri = [];
for i = 1:nq
    if(typ(i) == 3)
        tri = [tri;TP(i,1:3)];
    else
        tri = [tri;[TP(i,1),TP(i,2),TP(i,4)]];
        tri = [tri;[TP(i,2),TP(i,3),TP(i,4)]];
    end
end


tiskni(ro_ref*R,0,'Hustota');


function Wb = prevedNaUzly(TP,W,typ,np)
    
    nt = length(W(:,1));
    Wb = zeros(np,1);
    poc = zeros(np,1);
    for i = 1:nt
        for j = 1:typ(i)
            Wb(TP(i,j)) = Wb(TP(i,j)) + W(i,1);
            poc(TP(i,j)) = poc(TP(i,j)) + 1;
        end
    end
    
    for i = 1:np
        Wb(i) = Wb(i)/poc(i);
    end
    


function tiskni(Z,par,name)
global PX PY tri

figure('color','w');
trisurf(tri,PX,PY,Z);
shading interp;
osy = [min(PX) max(PX) min(PY) max(PY)];
axis equal;
axis(osy);
colorbar('SouthOutside');
xlabel('x [mm]','fontsize',14)
ylabel('y [mm]','fontsize',14)
set(gca,'fontsize',14);
title(name);

n = 20;
t = linspace(0,1,2000);
% c = (sign(sin(2*pi*n*t))/2 + 0.5)';
c = (sin(2*pi*n*t)/2 + 0.5)';
colormap([c,c,c]);













