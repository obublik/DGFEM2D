function zobraz2
global PX PY TP norm
[L_ref, u_ref, ro_ref, p_ref, eta_ref] = referencni_hodnoty;

We = load('We.txt');
typ = load('rad_elementu.txt');

load 'geometrie';
PXo = geometrie{1};
PYo = geometrie{2};
TPo = geometrie{3}+1;
np = length(PXo);
nt = length(TPo);

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

kapa = 1.302;
I = 1:length(PX);
R = Wb(I,1);
P = (kapa-1)*(Wb(I,4)-1./(2*Wb(I,1)).*(Wb(I,2).^2 + Wb(I,3).^2));
M =((Wb(I,2).^2 + Wb(I,3).^2)./(kapa*P(I).*Wb(I,1))).^(1/2);
U = Wb(I,2)./Wb(I,1);
V = Wb(I,3)./Wb(I,1);
T = p_ref/ro_ref*P./R/441;
% T = P./R/(kapa-1);

% tisk vysledku
tiskni(M,1);
tiskni(p_ref*P,0);
tiskni(T,0);
tiskni(ro_ref*R,0);

% proudnice
data = load('brit_1_04.mat');
data = data.data;
h = 0.0001;
figure('color','w');
hold on;
triplot(TP,PX,PY,'color',[0.8 0.8 0.8]);
Fu = TriScatteredInterp(PX,PY,U);
Fv = TriScatteredInterp(PX,PY,V);
[X0,Y0] = meshgrid(min(PX):h:max(PX), min(PY):h:max(PY));
U0 = Fu(X0,Y0);
V0 = Fv(X0,Y0);
in = inpolygon(X0,Y0,L_ref*data(:,1),L_ref*(data(:,2)-0.5));
U0(in == 0) = NaN;
V0(in == 0) = NaN;
ns = 20;
x0 = [-0.06*ones(1,ns), -6e-3*ones(1,ns), 0.02*ones(1,ns), 0.013];
y0 = linspace(min(PY),max(PY),ns);
y0 = [y0,y0,y0,-2.5e-3];
h = streamline(X0,Y0,U0,V0,x0,y0);
set(h, 'Color', 'red');
axis equal;
box on;

% tlak na stene
y_max = max(PY);
figure('color','w');
hold on;
for i = 1:np
    if(abs(PY(i)- y_max) < 1e-8)
        plot(PX(i),p_ref*P(i),'.');
    end
end
box on
grid on
xlabel('x [mm]')
ylabel('p [Pa]')
xlim([-0.06 0.06])
title('Prubeh tlaku na stene statoru');

% prutoky
S_in = 0;
S_out = 0;
plus = [2 3 1];
for i = 1:length(W(:,1));
    for j = 1:3
        if(TT(i,j) == -2)
            jp = plus(j);
            A = [PX(TP(i,j)),PY(TP(i,j))];
            B = [PX(TP(i,jp)),PY(TP(i,jp))];
            dl = norm(A-B);
            S_in = S_in + (Wb(TP(i,j),2)+Wb(TP(i,jp),2))/2*dl;
        end
         if(TT(i,j) == -3)
            jp = plus(j);
            A = [PX(TP(i,j)),PY(TP(i,j))];
            B = [PX(TP(i,jp)),PY(TP(i,jp))];
            dl = norm(A-B);
            S_out = S_out + (Wb(TP(i,j),2)+Wb(TP(i,jp),2))/2*dl;
        end
    end
end
display(['Vstup: ',num2str(S_in),', vystup: ',num2str(S_out),', rozdil: ',num2str(S_in-S_out)])


display(['Max Mach: ',num2str(max(M))]);

% figure;
% load 'reziduum.txt'
% semilogy(reziduum(:,1),reziduum(:,2));


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


function tiskni(Z,par)
global PX PY TP

figure('color','w');
[Cs,h] = tricontf(PX,PY,TP,Z,15);
if(par)
    set(h,'edgecolor','none');
end
osy = [min(PX) max(PX) min(PY) max(PY)];
axis(osy);
axis equal;
colorbar;
xlabel('x [mm]')
ylabel('y [mm]')
title('celkovy pohled');

figure('color','w');
[Cs,h] = tricontf(PX,PY,TP,Z,15);
if(par)
    set(h,'edgecolor','none');
end
axis equal;
osy = [-0.0135 -0.0095 0.5e-3 3.5e-3];
axis(osy);
box on;
colorbar;
xlabel('x [mm]')
ylabel('y [mm]')
title('brit 1');

figure('color','w');
[Cs,h] = tricontf(PX,PY,TP,Z,15);
if(par)
    set(h,'edgecolor','none');
end
axis equal;
osy = [-0.0135+0.0243 -0.0095+0.0243 0.5e-3 3.5e-3];
axis(osy);
box on;
colorbar;
xlabel('x [mm]')
ylabel('y [mm]')
title('brit 2');

figure('color','w');
[Cs,h] = tricontf(PX,PY,TP,Z,15);
if(par)
    set(h,'edgecolor','none');
end
osy = [-0.0135 -0.0095+0.023 0 3.5e-3];
axis(osy);
axis equal;
box on;
colorbar;
xlabel('x [mm]')
ylabel('y [mm]')
title('mezera mezi brity');










