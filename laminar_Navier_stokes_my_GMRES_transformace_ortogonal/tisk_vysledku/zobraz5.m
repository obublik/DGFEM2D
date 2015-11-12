function zobraz5
global norm
We = load('We.txt');
typ = load('rad_elementu.txt');

load 'geometrie';
PXo = geometrie{1};
PYo = geometrie{2};
TPo = geometrie{3}+1;
np = length(PXo);
nt = length(TPo);

osy = [min(PXo) max(PXo) min(PYo) max(PYo)];
norm = sqrt([1,1/6,0.5,1/15,1/9,1/3,1/4,1/12,1/20,1/28]);

n_troj = 4;
B = [3/4,1/8,1/8; 1/8,3/4,1/8; 1/8,1/8,3/4; 0.5 0.5 0.5];

% vytvoreni nove triangulace
EP = [TPo(:,[1,2]); TPo(:,[2,3]); TPo(:,[3,1])];
EP = unique(sort(EP,2),'rows');
ne = length(EP);

% vytvoreni matice stran
Er = sparse(np,np);
for i = 1:ne
    Er(EP(i,1),EP(i,2)) = i;
    Er(EP(i,2),EP(i,1)) = i;
end
TE = zeros(nt,3);
for i = 1:nt
    TE(i,1) = Er(TPo(i,1),TPo(i,2)); % index prvni strany
    TE(i,2) = Er(TPo(i,2),TPo(i,3));
    TE(i,3) = Er(TPo(i,3),TPo(i,1));
end

PX = [PXo;(PXo(EP(:,1))+PXo(EP(:,2)))/2];
PY = [PYo;(PYo(EP(:,1))+PYo(EP(:,2)))/2];

TP = zeros(4*nt,3);
W = zeros(4*nt,4);
k = 1;
for i=1:nt
    j = (4*(i-1));
    TP(j+1,:)=[TPo(i,1),TE(i,1)+np,TE(i,3)+np];
    TP(j+2,:)=[TE(i,1)+np,TPo(i,2),TE(i,2)+np];
    TP(j+3,:)=[TE(i,3)+np,TE(i,2)+np,TPo(i,3)];
    TP(j+4,:)=[TE(i,1)+np,TE(i,2)+np,TE(i,3)+np];
    
    switch typ(i)
        case 1
            nb = 1;
        case 2
            nb = 3;
        case 3
            nb = 6;
        case 4
            nb = 10;
    end
    
    for s = 1:n_troj % opakovani pres vnitrni trojuhelniky
        W_lok = [0 0 0 0];
        for m = 1:nb % opakovani pres bazove funkce
            bazf = baze(m,B(s,1),B(s,2),B(s,3));
            for p = 1:4  % opakovani pres vektor konzervativnich promennych
                W_lok(p) = W_lok(p) + We(k+m-1,p)*bazf;
            end
        end
        W(j+s,:) = W_lok;
    end
    k = k + nb;
end

Wb = prevedNaUzly(TP,W,length(PX));

np = length(PX);
kapa = 1.4;
I = 1:np;
R = Wb(I,1);
P = (kapa-1)*(Wb(I,4)-1./(2*Wb(I,1)).*(Wb(I,2).^2 + Wb(I,3).^2));
M =((Wb(I,2).^2 + Wb(I,3).^2)./(kapa*P(I).*R)).^(1/2);
% U = Wb(I,2)./Wb(I,1);
% V = Wb(I,3)./Wb(I,1);


figure('name','mach','color','w');
hold on;
triplot(TPo,PXo,PYo,'color',[0.8 0.8 0.8]);
tricontour(TP,PX,PY,M,50);
axis(osy);
axis equal;
colorbar;

figure('name','tlak','color','w');
hold on;
triplot(TPo,PXo,PYo,'color',[0.8 0.8 0.8]);
tricontour(TP,PX,PY,P,30);
axis(osy);
axis equal;
colorbar;

display(['Max Mach: ',num2str(max(M))]);


function Wb = prevedNaUzly(TP,W,np)
    
    nt = length(TP(:,1));
    Wb = zeros(np,4);
    poc = zeros(np,1);
    for i = 1:nt
        for j = 1:3
            Wb(TP(i,j),:) = Wb(TP(i,j),:) + W(i,:);
            poc(TP(i,j)) = poc(TP(i,j)) + 1;
        end
    end
    
    for i = 1:np
        Wb(i,:) = Wb(i,:)/poc(i);
    end

    
function b = baze(j,b1,b2,b3)
global norm
xi = -b1 + b2;
ni = -b1 - b2 + b3;
switch j
    case 1
        b = 1/norm(1);
    case 2
        b = xi/norm(2);
    case 3
        b = 0.5*(1+3*ni)/norm(3);
    case 4
        b = 0.125*(12*xi*xi-(ni-1)*(ni-1)/norm(4));
    case 5
        b = 0.5*xi*(3+5*ni)/norm(5);
    case 6
        b = (-0.5+ni+2.5*ni*ni)/norm(6);
    case 7
        b  = 0.125*xi*(20*xi^2 - 3*(ni-1)^2)/norm(7);
    case 8
        b = -1/16*(-12*xi^2+(-1+ni)^2)*(5+7*ni)/norm(8);
    case 9
        b = 1/4*xi*(1+18*ni+21*ni^2)/norm(9);
    case 10
        b = 0.125*(-3-15*ni+15*ni^2+35*ni^3)/norm(10);
end



