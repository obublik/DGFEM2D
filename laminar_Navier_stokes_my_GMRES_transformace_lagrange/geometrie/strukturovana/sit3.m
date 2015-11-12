function sit3
n1 = 50;
n2 = 20;

X = zeros(n1,n2);
Y = zeros(n1,n2);

x = nlinspace(0,3,n1);
for i = 1:n1
    X(i,1) = x(i);
    X(i,n2) = x(i);
    
    Y(i,1) = -f(x(i));
    Y(i,n2) = f(x(i));
    
    t = linspace(0,1,n2);
    for j = 2:n2-1
        X(i,j) = X(i,1) + (X(i,n2)-X(i,1))*t(j);
        Y(i,j) = Y(i,1) + (Y(i,n2)-Y(i,1))*t(j);
    end
end

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

k = 1;
TP = zeros((n1-1)*(n2-1),4);
for i = 1:n1-1
    for j = 1:n2-1
        TP(k,:) = [Ip(i,j),Ip(i+1,j),Ip(i+1,j+1),Ip(i,j+1)];
        k = k + 1; 
    end
end

P = [PX,PY];
Q = TP;

save 'P' P
save 'Q' Q

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
y = 0.5 - 0.2*exp(-(x-1)^2/0.1);

function y = nlinspace(a,b,n)
y = zeros(n,1);
for i = 1:n
    x = (i-1)/(n-1);
    y(i+1) = y(i) + 1-0.999*exp(-(3*x-1.31)^2/100);
end
y = y - min(y);
y = y/(max(y));
y = a + y*(b-a);



