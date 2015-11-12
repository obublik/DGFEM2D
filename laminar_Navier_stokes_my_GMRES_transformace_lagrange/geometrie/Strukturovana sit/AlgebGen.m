function [Xp,Yp] = AlgebGen(X,Y)
n1 = length(X(:,1));
n2 = length(X(1,:));
F = zeros(n1,1);
N = zeros(n2,1);

for i = 1:n1
    F(i) = (i-1)/(n1-1);
end

for j = 1:n2
    N(j) = (j-1)/(n2-1);
end

Xp = zeros(n1,n2);
Yp = zeros(n1,n2);

for i = 1:n1
    for j = 1:n2
        Xp(i,j) = (1-F(i))*X(1,j)+F(i)*X(n1,j)+(1-N(j))*X(i,1)+N(j)*X(i,n2)-(1-F(i))*(1-N(j))*X(1,1)-F(i)*(1-N(j))*X(n1,1)-(1-F(i))*N(j)*X(1,n2)-F(i)*N(j)*X(n1,n2);
        Yp(i,j) = (1-F(i))*Y(1,j)+F(i)*Y(n1,j)+(1-N(j))*Y(i,1)+N(j)*Y(i,n2)-(1-F(i))*(1-N(j))*Y(1,1)-F(i)*(1-N(j))*Y(n1,1)-(1-F(i))*N(j)*Y(1,n2)-F(i)*N(j)*Y(n1,n2);
    end
end