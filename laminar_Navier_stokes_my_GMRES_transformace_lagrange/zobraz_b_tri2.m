function zobraz_b_tri2

% zobrazuje vysledky
We = load('We.txt');
rad = load('rad_elementu.txt');
load 'geometrie';
PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3}+1;
typ = geometrie{5};
nt = length(TP(:,1));

% nb = length(We(:,1))/nt;

figure('color','w');
hold on;
k = 1;
M_max = 0;
for i = 1:nt
    nb = baze_poc(rad(i));
    n = typ(i);
    I = 1:n;
    X0 = PX(TP(i,I));
    Y0 = PY(TP(i,I));
    
    b = barycentr(rad(i));
    mb = length(b(:,1));
    xb = zeros(mb,1);
    yb = zeros(mb,1);
    for j = 1:mb
        xb(j) = b(j,1)*X0(1) + b(j,2)*X0(2) + b(j,3)*X0(3);
        yb(j) = b(j,1)*Y0(1) + b(j,2)*Y0(2) + b(j,3)*Y0(3);
    end
    TP_loc = delaunay(xb,yb);
    
    m = 1:nb;
    W = We(k+m-1,:);
    if(rad == 1)
        W = [W;W;W];
    end
    
    kapa = 1.4;
    P = (kapa-1)*(W(:,4)-1./(2*W(:,1)).*(W(:,2).^2 + W(:,3).^2));
    M =((W(:,2).^2 + W(:,3).^2)./(kapa*P(:).*W(:,1))).^(1/2);
    
    trisurf(TP_loc,xb,yb,M,'linestyle','none');
    
    k = k + nb;
    
    if(max(M) > M_max)
        M_max = max(M);
    end
end
box on;
colorbar;
shading interp
axis equal

J = [1,2,3,1];
for i = 1:length(TP(:,1))
    plot3(PX(TP(i,J)),PY(TP(i,J)),M_max*ones(4,1),'color','k');
end
view(2);


function nb = baze_poc(rad)
nb = rad*(rad+1)/2;


function b = barycentr(n)
n = max(n,2);
x = linspace(0,1,n);
b = [];
s = 1;
for i = 1:n
    y = linspace(0,x(i),i);
    for j = 1:i
        b(s,1) = 1-x(i);
        b(s,2) = y(j);
        b(s,3) = x(i) - y(j);
        s = s + 1;
    end
end







