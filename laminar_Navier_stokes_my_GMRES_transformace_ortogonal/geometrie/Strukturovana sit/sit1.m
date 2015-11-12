function sit1
n1 = 38;
n2 = 12;

X = zeros(n1,n2);
Y = zeros(n1,n2);

for i = 1:n1
    x = 3*(i-1)/(n1-1);
    X(i,1) = x;
    X(i,n2) = x;
    
    Y(i,1) = f(x);
    Y(i,n2) = 1;
end
    
for j = 1:n2
    X(1,j)  = 0;
    X(n1,j) = 3;
    
    y = (j-1)/(n2-1);
    Y(1,j)  = y;
    Y(n1,j) = y;
end

[X,Y] = AlgebGen(X,Y);

PX = zeros(n1*n2,1);
PY = zeros(n1*n2,1);
Ip = zeros(n1,n2);
k = 1;
for i = 1:n1
    for j = 1:n2
        PX(k) = X(i,j);
        PY(k) = Y(i,j);
        Ip(i,j) = k;
        k = k + 1; 
    end
end
Ip = Ip - 1;

k = 1;
TP = zeros((n1-1)*(n2-1),4);
Iq = zeros(n1,n2);
for i = 1:n1-1
    for j = 1:n2-1
        TP(k,:) = [Ip(i+1,j),Ip(i+1,j+1),Ip(i,j+1),Ip(i,j)];
        Iq(i,j) = k;
        k = k + 1; 
    end
end
Iq = Iq - 1;

k = 1;
TT = zeros((n1-1)*(n2-1),4);
for i = 1:n1-1
    for j = 1:n2-1
        if(i > 1)
            TT(k,3) = Iq(i-1,j);
        else
            TT(k,3) = -2;
        end
        
        if(i < n1-1)
            TT(k,1) = Iq(i+1,j);
        else
            TT(k,1) = -3;
        end
        
        if(j > 1)
            TT(k,4) = Iq(i,j-1);
        else
            TT(k,4) = -1;
        end
        
        if(j < n2-1)
            TT(k,2) = Iq(i,j+1);
        else
            TT(k,2) = -1;
        end
        
        k = k + 1; 
    end
end

geometrie{1} = PX;
geometrie{2} = PY;
geometrie{3} = TP;
geometrie{4} = TT;
geometrie{5} = 4*ones(length(TP(:,1)),1);

save 'geometrie' geometrie

figure
hold on;
for i = 1:length(X(:,1))
    plot(X(i,:),Y(i,:));
end
for j = 1:length(X(1,:))
    plot(X(:,j),Y(:,j));
end
axis equal;

function y = f(x)
h = 0.1;
r = h/2 + 1/(8*h);
ys = 0.1-r;
    if(x < 1 || x > 2)
        y = 0;
    else
        y = sqrt(r^2-(x-1.5)^2) + ys;
    end