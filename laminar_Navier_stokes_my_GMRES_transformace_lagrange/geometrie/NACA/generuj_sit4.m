function generuj_sit4
np = 80; % pocet elementu profilu
n1 = 2*np;
n2 = 30;

X = zeros(n1,n2);
Y = zeros(n1,n2);

data = NACA2(np);

I = 1:n1;
X(I,1) = data(I,1);
Y(I,1) = data(I,2);

[X,Y] = Hyp_gen4(X,Y);

alfa = -1.25/180*pi;
X = cos(alfa)*X - sin(alfa)*Y;
Y = sin(alfa)*X + cos(alfa)*Y;

X = [X;X(1,:)];
Y = [Y;Y(1,:)];

figure;
hold on;
for i = 1:n1
    plot(X(i,:),Y(i,:),'k');
end
for j = 1:n2
    plot([X(:,j);X(1,j)],[Y(:,j);Y(1,j)],'k');
end
axis equal;

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


%__________________________________________________________________________
function [Xp,Yp] = Hyp_gen4(X,Y)
n1 = length(X(:,1));
n2 = length(X(1,:));
K = 10;
lambda = 4;

de = 1/n1;
dn = 1/n2;

xp = X(:,1);
yp = Y(:,1);
xn = xp;
yn = yp;
for op = 1:100
    for i = 1:n1
        ip = mod(i+1-1,n1)+1;
        im = mod(i-1-1,n1)+1;
        xn(i) = (xp(ip) + xp(im))/2;
        yn(i) = (yp(ip) + yp(im))/2;
    end
    xp = xn;
    yp = yn;
end


% explicitni
for j = 1:2
    for i = 1:n1
        ip = mod(i+1-1,n1)+1;
        im = mod(i-1-1,n1)+1;
        
        if(j > 1)
            nx = X(ip,j)-X(im,j);
            ny = Y(ip,j)-Y(im,j);
        else
            nx = xp(ip)-xp(im);
            ny = yp(ip)-yp(im);
        end
        
        g11 = ((nx/(2*de)).^2 + (ny/(2*de)).^2);
        
        V = K*sqrt(g11)*exp(-lambda*(1-j/n2));

        X(i,j+1) = X(i,j) - dn*(V/g11)*ny/(2*de);
        Y(i,j+1) = Y(i,j) + dn*(V/g11)*nx/(2*de);
    end
end

% implicitni
A = sparse(2*n1,2*n1);
b = zeros(2*n1,1);
koef = dn/(2*de);
I1 = 1:2:2*n1;
I2 = 2:2:2*n1;
for j = 3:n2-1
    for i = 1:n1
        if(i == n1)
            ip = 1;
        else
            ip = i+1;
        end
        if(i == 1)
            im = n1;
        else
            im = i-1;
        end
        
        dxde = (X(ip,j)-X(im,j))/(2*de);
        dyde = (Y(ip,j)-Y(im,j))/(2*de);
        
        g11 = dxde^2 + dyde^2;
        
        V0 = K*sqrt(g11)*exp(-lambda*(1-j/n2));
        
        dxdn = -dyde*V0/g11;
        dydn = dxde*V0/g11;
        
        ap = dxde*dxdn - dyde*dydn;
        bp = dxde*dydn + dxdn*dyde;
        cp = g11;
        
        C = 1/cp*[ap bp; bp -ap];
        tlum = -sqrt((ap^2 + bp^2)/cp^2);
        
        V = K*sqrt(g11)*exp(-lambda*(1-(j+1)/n2));
        
        S = (V+V0)/cp*[-dyde; dxde];
        
        I = 2*(i-1)+1:2*i;
        Ip = 2*(ip-1)+1:2*ip;
        Im = 2*(im-1)+1:2*im;
        
        A(I,I) = (1-2*tlum)*eye(2);
        A(I,Ip) = koef*C + tlum*eye(2);
        A(I,Im) = -koef*C + tlum*eye(2);
        b(I) = dn*S + [X(i,j);Y(i,j)];
    end
    R = A\b;

    X(:,j+1) = R(I1);
    Y(:,j+1) = R(I2);
end

Xp = X;
Yp = Y;


%__________________________________________________________________________
function data = NACA2(n)

designation = '0012';

t = str2double(designation(3:4))/100;
m = str2double(designation(1))/100;
p = str2double(designation(2))/10;

a0 = 0.2969;
a1 =-0.1260;
a2 =-0.3516;
a3 = 0.2843;
a4 = -0.1036; %-0.1015;

% neekvidistantni deleni
c = 1;
x = linspace(0,1,n+2)';
x2 = zeros(n+1,1);
for i = 2:n+1
    x2(i) = x2(i-1) + (x(i)*(1-x(i)))^2;
end
x = c*x2/max(x2);

yt = (t/0.2)*(a0*sqrt(x)+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);

xc1 = x(find(x<=p));
xc2 = x(find(x>p));
xc=[xc1 ; xc2];

if(p == 0)
    xu = x;
    yu = yt;

    xl = x;
    yl = -yt;
else
    yc1 = (m/p^2)*(2*p*xc1-xc1.^2);
    yc2 = (m/(1-p)^2)*((1-2*p)+2*p*xc2-xc2.^2);
    yc = [yc1 ; yc2];

    dyc1_dx = (m/p^2)*(2*p-2*xc1);
    dyc2_dx = (m/(1-p)^2)*(2*p-2*xc2);
    dyc_dx = [dyc1_dx ; dyc2_dx];
    theta = atan(dyc_dx);

    xu = x - yt.*sin(theta);
    yu = yc + yt.*cos(theta);

    xl = x + yt.*sin(theta);
    yl = yc - yt.*cos(theta);
end

X = [flipud(xl) ; xu(2:end)];
Y = [flipud(yl) ; yu(2:end)];

% aby byly koncove body na ose x
Y(1) = 0;
Y(length(Y)) = 0;

data = [X,Y];
