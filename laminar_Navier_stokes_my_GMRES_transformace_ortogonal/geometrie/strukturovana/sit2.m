function sit2
global geom poc;

geom = {};
poc = 0;

n1 = 50;
n2 = 30;

dy = 0.0005;

obdelnik([0,0],[1,0],[1,0.05],[0,0.05],n1,n2,1,3,1,3,dy,dy);
% obdelnik([0,0.05],[1,0.05],[1,0.1],[0,0.1],n1,n2,1,0,1,0,dy,dy);

[P,Q] = vytvor_sit;
kresli;
save 'P' P;
save 'Q' Q;

display(['Pocet bodu site je: ',num2str(length(P(:,1)))]);


function obdelnik(A,B,C,D,n1,n2,z1,z2,z3,z4,dx0,dy0)
% A, B, C, D jsou krajni body obdelnika [x,y]
% n1, n2 jsou pocty bunek v x-ove a y-ove soureadnici
% z1, z2 jsou zahusteni, 0 - doleva, 1 - nic, 2 - doprava, 3 - obe strany
global geom poc;
    X = zeros(n1+1,n2+1);
    Y = zeros(n1+1,n2+1);
    
    a = deleni(dx0,n1,norm(B-A),z1);
    b = deleni(dy0,n2,norm(C-B),z2);
    c = deleni(dx0,n1,norm(D-C),z3);
    d = deleni(dy0,n2,norm(A-D),z4);
    
    for i = 1:n1+1
        for j = 1:n2+1
            f1 = (1-a(i))*(1-d(j));
            f2 = a(i)*(1-b(j));
            f3 = c(i)*b(j);
            f4 = (1-c(i))*d(j);
            X(i,j) = A(1)*f1 + B(1)*f2 + C(1)*f3 + D(1)*f4;
            Y(i,j) = A(2)*f1 + B(2)*f2 + C(2)*f3 + D(2)*f4;
        end
    end
    
    obd{1} = n1;
    obd{2} = n2;
    obd{3} = X;
    obd{4} = Y;
    poc = poc + 1;
    geom{poc} = obd;

    
function kresli
global geom poc;
figure;
hold on;
for i = 1:poc
    obd = geom{i};
    n1 = obd{1};
    n2 = obd{2};
    X = obd{3};
    Y = obd{4};
    I = 1:n1+1;
    J = 1:n2+1;
    for i = 1:n1+1
        plot(X(i,J),Y(i,J));
    end
    for j = 1:n2+1
        plot(X(I,j),Y(I,j));
    end
end
axis equal


function [P,Q] = vytvor_sit
global geom poc;

% zjednoznacneni bodu
Pp = [];
s = 1;
for k = 1:poc
    obd = geom{k};
    n1 = obd{1};
    n2 = obd{2};
    X = obd{3};
    Y = obd{4};
    np = (n1+1)*(n2+1);
    Ppom = zeros(np,5);
    p = 1;
    for i = 1:(n1+1)
        for j = 1:(n2+1)
            Ppom(p,:) = [X(i,j),Y(i,j),s,0,0];
            p = p + 1;
            s = s + 1;
        end
    end
    Pp = [Pp;Ppom];
end
Pp = myqs(Pp);

Pp(1,4) = 1;
Pp(1,5) = 1;
for i = 2:length(Pp(:,1))
    if(norm(Pp(i,1:2)-Pp(i-1,1:2)) < 1e-8)
        Pp(i,4) = Pp(i-1,4);
    else
        Pp(i,4) = Pp(i-1,4) + 1;
        Pp(i,5) = 1;
    end
end

P = Pp(Pp(:,5) == 1, 1:2);
I = Pp(:,3);
slov(I) = Pp(:,4);

% vytvoreni matice indexu
s = 1;
for k = 1:poc
    obd = geom{k};
    n1 = obd{1};
    n2 = obd{2};
    ind = zeros(n1+1,n2+1);
    
    for i = 1:(n1+1)
        for j = 1:(n2+1)
            ind(i,j) = slov(s);
            s = s + 1;
        end
    end
    obd{5} = ind;
    geom{k} = obd;
end

% vytvoreni matice Q
Q = [];
for k = 1:poc
    obd = geom{k};
    n1 = obd{1};
    n2 = obd{2};
    ind = obd{5};
    
    nq = n1*n2;
    Qpom = zeros(nq,4);
    p = 1;
    for i = 1:n1
        for j = 1:n2
            Qpom(p,:) = [ind(i,j),ind(i+1,j),ind(i+1,j+1),ind(i,j+1)];
            p = p + 1;
        end
    end
    Q = [Q;Qpom];
end

function y = deleni(dy0,n,L,par)
dy0 = dy0/L;
dy = zeros(n,1);
dy(1) = dy0;
y = zeros(n+1,1);
y(2) = dy0;

if(par == 1)
    y = linspace(0,1,n+1)';
    return;
end

if(par == 3)
    lich = mod(n,2);
    n = fix(n/2);
    dy0 = 2*dy0;
    dy = zeros(n,1);
    dy(1) = dy0;
    y = zeros(n+1,1);
    y(2) = dy0;
    y = zeros(n+1,1);
    y(2) = dy0;
end

a = 0;
b = 10;
c = 0.3;
s = 1;
while 1
    A = linspace(c,0,n);
    for i = 2:n
        dy(i) = dy(i-1)*(1+A(i));
        y(i+1) = y(i) + dy(i);
    end
    if(y(n+1) < 1)
        a = c;
    else
        b = c;
    end
    c = (a+b)/2;
    if((b-a)/2 < 1e-8)
        break
    end
    s = s+1;
    if(s > 10000)
        display('Chyba')
        break;
    end
end
y = y/y(n+1);

if(par == 2)
    y = 1-y(end:-1:1);
end

if(par == 3)
    if(lich == 0)
        y = [y;2-y((end-1):-1:1)]/2;
    else
        dy = y(end)-y(end-1);
        y1 = [y;y(end)+dy];
        y2 = 2+dy-y((end-1):-1:1);
        y = [y1;y2];
        y = y/y(end);
    end
end


function out = myqs(data)
if(isempty(data) ~= 1)
    n = length(data(:,1));
    if(n < 2)
        out = data;
    else
        ind = floor(n/2);
        p = data(ind,:);
        L = zeros(size(data));
        iL = 1;
        R = zeros(size(data));
        iR = 1;
        if(n > 1)
            for i = 1:n
                if(i ~= ind)
                    if(je_mensi(data(i,:),p))
                        L(iL,:) = data(i,:);
                        iL = iL + 1;
                    else
                        R(iR,:) = data(i,:);
                        iR = iR + 1;
                    end
                end
            end
            L = myqs(L(1:iL-1,:));
            R = myqs(R(1:iR-1,:));
            out = [L; p; R];
        end
    end
else
    out = [];
end

function c = je_mensi(a,b)
c = 0;
if(b(1) - a(1) > 1e-5)
    c = 1;
elseif(abs(a(1) -b(1)) < 1e-5)
    if(b(2) - a(2) > 1e-5)
        c = 1;
    end
end



