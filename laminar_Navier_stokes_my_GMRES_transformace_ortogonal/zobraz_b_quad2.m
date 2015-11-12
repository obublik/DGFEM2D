function zobraz_b_quad2
global nb A_baze

% zobrazuje vysledky
A_baze = load('A_baze_quad.txt');
We = load('We.txt');
rad = load('rad_elementu.txt');
load 'geometrie';
PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3}+1;
typ = geometrie{5};
nt = length(TP(:,1));

figure('color','w');
hold on;
k = 1;
M_max = 0;
for i = 1:nt
    nb = baze_poc(rad(i));
    baze_x = zeros(1,nb);
    baze_y = zeros(1,nb);
    s = 1;
    for r = 1:rad(i)
        for j = 1:r
            baze_x(s) = r-j;
            baze_y(s) = j-1;
            s = s + 1;
        end
    end
    n = typ(i);
    I = 1:n;
    X0 = PX(TP(i,I));
    Y0 = PY(TP(i,I));
    
    b = barycentr(4);
    mb = length(b(:,1));
    Xi = zeros(2*mb,1);
    Eta = zeros(2*mb,1);
    x = [0 1 1 0];
    y = [0 0 1 1];
    for j = 1:mb
        Xi(j) = b(j,1)*x(1) + b(j,2)*x(2) + b(j,3)*x(3);
        Eta(j) = b(j,1)*y(1) + b(j,2)*y(2) + b(j,3)*y(3);
        Xi(mb+j) = b(j,1)*x(1) + b(j,2)*x(3) + b(j,3)*x(4);
        Eta(mb+j) = b(j,1)*y(1) + b(j,2)*y(3) + b(j,3)*y(4);
    end
    warning off;
    TP_loc = delaunay(Xi,Eta);
    
    V = [0,0,0,1; 0,1,0,1; 1 1 1 1; 0,0,1,1];
    A = V\X0;
    B = V\Y0;
    
    X = zeros(2*mb,1);
    Y = zeros(2*mb,1);
    for j = 1:(2*mb)
        X(j) = [Xi(j)*Eta(j),Xi(j),Eta(j),1]*A;
        Y(j) = [Xi(j)*Eta(j),Xi(j),Eta(j),1]*B;
    end
    
    W = zeros(length(X),4);
    for q = 1:4
        for j = 1:length(X)
            for m = 1:nb
                W(j,q) = W(j,q) + We(k+m-1,q)*baze(nb,m,Xi(j),Eta(j));
            end
        end
    end
    
    kapa = 1.4;
    I = 1:length(X);
    P = (kapa-1)*(W(I,4)-1./(2*W(I,1)).*(W(I,2).^2 + W(I,3).^2));
    M =((W(I,2).^2 + W(I,3).^2)./(kapa*P(I).*W(I,1))).^(1/2);
    
    trisurf(TP_loc,X,Y,M,'linestyle','none');
    
    k = k + nb;
    
    if(max(M) > M_max)
        M_max = max(M);
    end
end
box on;
colorbar;
shading interp
axis equal

J = [1,2,3,4,1];
for i = 1:length(TP(:,1))
    plot3(PX(TP(i,J)),PY(TP(i,J)),M_max*ones(5,1),'color','k');
end
view(2);


function b = baze(nb,m,x,y)
global A_baze
b = 0;
for j = 1:nb
    b = b + A_baze(m,j)*x^rx(j)*y^ry(j);
end


function nb = baze_poc(rad)
switch rad
    case 1
        nb = 1;
    case 2
        nb = 3;
    case 3
        nb = 6;
    case 4
        nb = 10;
    case 5
        nb = 15;
end


function b = barycentr(n)
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

function r = rx(n)
i = 1;
s = 1;
while 1
    for j = i:-1:1
        r = j-1;
        if(s == n)
            return;
        end
        s = s + 1;
    end
    i = i + 1;
end

function r = ry(n)
i = 1;
s = 1;
while 1
    for j = 1:i
        r = j-1;
        if(s == n)
            return;
        end
        s = s + 1;
    end
    i = i + 1;
end