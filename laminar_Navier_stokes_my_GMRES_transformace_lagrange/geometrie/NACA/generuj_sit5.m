function generuj_sit5
np = 40; % pocet elementu profilu
ns = 5; % pocet elementu symetrie
n1 = 2*np + 2*ns;
n2 = 25;

X = zeros(n1+1,n2);
Y = zeros(n1+1,n2);

data = NACA(np);
xs1 = linspace(1,0,ns+1)';
xs2 = linspace(0,1,ns+1)';

% spodni stena
X(:,1) = [xs1;data(2:2*np,1)-1;xs2];
Y(:,1) = [zeros(ns+1,1);data(2:2*np,2);zeros(ns+1,1)];

[X,Y] = Hyp_gen2(X,Y);

% posunuti krajnich bodu
min_y = min(Y(:,n2));
max_y = max(Y(:,n2));
for i = 1:length(Y(:,n2))
    if(X(i,n2) > -0.7)
        if(Y(i,n2) > 0)
            Y(i,n2) = max_y;
        end
        if(Y(i,n2) < 0)
            Y(i,n2) = min_y;
        end
    end
end

% prevod na trojuhelniky
[n1,n2] = size(X);
ind = zeros(n1,n2);
P = zeros(n1*n2,2);
s = 1;
for i = 1:n1
    for j = 1:n2
        ind(i,j) = s;
        P(s,1) = X(i,j);
        P(s,2) = Y(i,j);
        s = s + 1;
    end
end

TP = zeros(2*(n1-1)*(n2-1),3);
s = 1;
for i = 1:n1-1
    for j = 1:n2-1
        TP(s,:) = [ind(i,j),ind(i+1,j),ind(i+1,j+1)];
        TP(s+1,:) = [ind(i,j),ind(i+1,j+1),ind(i,j+1)];
        s = s + 2;
    end
end

% unifikace
[P,I,J] = uniquetol(P+1,1e-9,'rows');
P = P-1;
TP = [J(TP(:,1)),J(TP(:,2)),J(TP(:,3))];

mesh{1} = P;
mesh{2} = TP;

save mesh mesh;

figure;
triplot(TP,P(:,1),P(:,2),'k');
axis equal;

% PX = P(:,1);
% PY = P(:,2);
% ind = 1:length(PX);
% ind = ind(PY == 0 & PX > 0);
% P(ind,:)

%__________________________________________________________________________
function [Xp,Yp] = Hyp_gen2(X,Y)
n1 = length(X(:,1));
n2 = length(X(1,:));
K = 10;
lambda = 4;

de = 1/n1;
dn = 1/n2;

a = 0.5;
for j = 1:n2-1
    for i = 2:n1-1
        dx = (X(i+1,j)-X(i-1,j))/(2*de);
        dy = (Y(i+1,j)-Y(i-1,j))/(2*de);
        
        g11 = dx.^2 + dy.^2;
        
        V = K*sqrt(g11)*exp(-lambda*(1-j/n2));

        X(i,j+1) =  a*X(i,j) + (1-a)/2*(X(i+1,j) + X(i-1,j)) - dn*(V/g11)*dy;
        Y(i,j+1) =  a*Y(i,j) + (1-a)/2*(Y(i+1,j) + Y(i-1,j)) + dn*(V/g11)*dx;
    end
    g11 = (((X(2,j)-X(1,j))/de).^2 + ((Y(2,j)-Y(1,j))/de).^2);
    V = K*sqrt(g11)*exp(-lambda*(1-j/n2));
    X(1,j+1) = X(1,j) - dn*(V/g11)*(Y(2,j) - Y(1,j))/de;
    Y(1,j+1) = Y(1,j) + dn*(V/g11)*(X(2,j) - X(1,j))/de;
        
    g11 = (((X(n1,j)-X(n1-1,j))/de).^2 + ((Y(n1,j)-Y(n1-1,j))/de).^2);
    V = K*sqrt(g11)*exp(-lambda*(1-j/n2));
    X(n1,j+1) = X(n1,j) - dn*(V/g11)*(Y(n1,j) - Y(n1-1,j))/de;
    Y(n1,j+1) = Y(n1,j) + dn*(V/g11)*(X(n1,j) - X(n1-1,j))/de;
end

Xp = X;
Yp = Y;
