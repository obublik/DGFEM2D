function zobraz
global PX PY tri
[L_ref, u_ref, ro_ref, p_ref, eta_ref, k_ref, Re, Pr, kapa] = referencni_hodnoty;

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

Wb = prevedNaUzly(TP,W,typ,np);

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

I = 1:np;
R = Wb(I,1);
P = (kapa-1)*(Wb(I,4)-1./(2*Wb(I,1)).*(Wb(I,2).^2 + Wb(I,3).^2));
M =((Wb(I,2).^2 + Wb(I,3).^2)./(kapa*P(I).*Wb(I,1))).^(1/2);
U = Wb(I,2)./Wb(I,1);
V = Wb(I,3)./Wb(I,1);
T = p_ref/ro_ref*P./R/441;
K = Wb(I,5)./Wb(I,1);
OM = Wb(I,6)./Wb(I,1);
mut = log(R.*K./OM);

% tisk vysledku
tiskni(M,1,'Mach');
tiskni(u_ref*sqrt(U.^2+V.^2),0,'Rychlost');
tiskni(p_ref*P,0,'Tlak');
tiskni(T,0,'Teplota');
tiskni(ro_ref*R,0,'Hustota');
% tiskni(k_ref*K,0,'k');
% tiskni(ro_ref*k_ref/eta_ref*OM,0,'omega');
% tiskni(eta_ref*mut,0,'\eta_{turb}');

% % tlak na stene
% y_max = max(PY);
% figure('color','w');
% hold on;
% x = [];
% px = [];
% max_py = max(PY(PX >0.05) );
% for i = 1:np
%     if(abs(PY(i) - max_py) < 1e-8)
%         x = [x,PX(i)];
%         px = [px,p_ref*P(i)];
%     end
% end
% plot(x,px,'color','k','marker','o','markersize',6,'markerfacecolor','k','linestyle','none');
% box on
% grid on
% xlabel('x [mm]','fontsize',14)
% ylabel('p [Pa]','fontsize',14)
% set(gca,'fontsize',14);

% tlak na stene
figure('color','w');
hold on;
x_min = min(PX);
x_max = max(PX);
y_min = min(PY);
y_max = max(PY);
x = linspace(x_min,x_max,100);
y = (y_max+y_min)/2*ones(size(x));
F = TriScatteredInterp(PX,PY,P);
px = F(x,y);
plot(x,px,'color','k','marker','o','markersize',4,'markerfacecolor','k','linestyle','none');
box on
grid on
xlabel('x [mm]','fontsize',14)
ylabel('p [Pa]','fontsize',14)
set(gca,'fontsize',14);
axis([x_min x_max min(P) max(P)]);

data = [x',px'];
fid = fopen('data.txt','w');
fprintf(fid,'%12.6f  %12.6f\n',data');
fclose(fid);

% mach na stene
figure('color','w');
hold on;
x_min = min(PX);
x_max = max(PX);
y_min = min(PY);
y_max = max(PY);
x = linspace(x_min,x_max,100);
y = (y_max+y_min)/2*ones(size(x));
F = TriScatteredInterp(PX,PY,M);
px = F(x,y);
plot(x,px,'color','k','marker','o','markersize',4,'markerfacecolor','k','linestyle','none');
box on
grid on
xlabel('x [mm]','fontsize',14)
ylabel('p [Pa]','fontsize',14)
set(gca,'fontsize',14);
axis([x_min x_max min(M) max(M)]);


display(['Max Mach: ',num2str(max(M))]);

% prutoky
S_in = 0;
S_out = 0;
plus = [2 3 4 1];
for i = 1:length(W(:,1));
    for j = 1:4
        if(TT(i,j) == -2)
            jp = plus(j);
            A = [PX(TP(i,j)),PY(TP(i,j))];
            B = [PX(TP(i,jp)),PY(TP(i,jp))];
            dl = norm(A-B);
            S_in = S_in + W(i,2)*dl;
        end
         if(TT(i,j) == -3)
            jp = plus(j);
            A = [PX(TP(i,j)),PY(TP(i,j))];
            B = [PX(TP(i,jp)),PY(TP(i,jp))];
            dl = norm(A-B);
            S_out = S_out + W(i,2)*dl;
        end
    end
end
% nasobeni delkou kruznice
Q_prevod = ro_ref*u_ref*pi*(0.156+0.150);
S_in = S_in*Q_prevod;
S_out = S_out*Q_prevod;
display(['Vstup: ',num2str(S_in),', vystup: ',num2str(S_out),', rozdil: ',num2str(S_in-S_out)])

% figure;
% load 'reziduum.txt'
% semilogy(reziduum(:,1),reziduum(:,2));


function Wb = prevedNaUzly(TP,W,typ,np)
    
    nt = length(W(:,1));
    Wb = zeros(np,6);
    poc = zeros(np,1);
    I = 1:4;
    for i = 1:nt
        for j = 1:typ(i)
            Wb(TP(i,j),I) = Wb(TP(i,j),I) + W(i,:);
            poc(TP(i,j)) = poc(TP(i,j)) + 1;
        end
    end
    
    for i = 1:np
        Wb(i,:) = Wb(i,:)/poc(i);
    end
    


function tiskni(Z,par,name)
global PX PY tri

figure('color','w');
[Cs,h] = tricontf(PX,PY,tri,Z,30);
if(par)
    set(h,'edgecolor','none');
end
osy = [min(PX) max(PX) min(PY) max(PY)];
axis equal;
axis(osy);
colorbar('SouthOutside');
xlabel('x [mm]','fontsize',14)
ylabel('y [mm]','fontsize',14)
set(gca,'fontsize',14);
title(name);
saveas(gca,['data\',name,'.png'],'png');


figure('color','w');
[Cs,h] = tricontf(PX,PY,tri,Z,50);
if(par)
    set(h,'edgecolor','none');
end
osy = [0.095, 0.108, -0.1-0.01, -0.1+0.01];
axis equal;
axis(osy);
colorbar('SouthOutside');
xlabel('x [mm]','fontsize',14)
ylabel('y [mm]','fontsize',14)
set(gca,'fontsize',14);
title(name);
saveas(gca,['data\',name,'_detail.png'],'png');












