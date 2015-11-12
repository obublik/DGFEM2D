function [Xp,Yp] = ElipticGen(X,Y)
Xp = X;
Yp = Y;
n1 = length(X(:,1));
n2 = length(X(1,:));
g11 = zeros(n1,n2);
g22 = zeros(n1,n2);
g12 = zeros(n1,n2);
xtemp = zeros(n1,n2);
ytemp = zeros(n1,n2);

[X,Y] = AlgebGen(X,Y);
I = 2:n1-1;
J = 2:n2-1;
n = 300;
for op = 1:n
    g11(I,J) = ((X(I+1,J)-X(I-1,J)).^2 + (Y(I+1,J)-Y(I-1,J)).^2)/4;
    g22(I,J) = ((X(I,J+1)-X(I,J-1)).^2 + (Y(I,J+1)-Y(I,J-1)).^2)/4;
    g12(I,J) = ((X(I+1,J)-X(I-1,J)).*(X(I,J+1)-X(I,J-1)) + (Y(I+1,J)-Y(I-1,J)).*(Y(I,J+1)-Y(I,J-1)))/4;
    
    xtemp(I,J) = 1./(2*(g11(I,J)+g22(I,J))).*(g22(I,J).*X(I+1,J)-0.5*g12(I,J).*X(I+1,J+1)+0.5*g12(I,J).*X(I+1,J-1)+g11(I,J).*X(I,J+1)+g11(I,J).*X(I,J-1)+g22(I,J).*X(I-1,J)-0.5*g12(I,J).*X(I-1,J-1)+0.5*g12(I,J).*X(I-1,J+1));
    ytemp(I,J) = 1./(2*(g11(I,J)+g22(I,J))).*(g22(I,J).*Y(I+1,J)-0.5*g12(I,J).*Y(I+1,J+1)+0.5*g12(I,J).*Y(I+1,J-1)+g11(I,J).*Y(I,J+1)+g11(I,J).*Y(I,J-1)+g22(I,J).*Y(I-1,J)-0.5*g12(I,J).*Y(I-1,J-1)+0.5*g12(I,J).*Y(I-1,J+1));
    
    err = sum(sum((X(I,J)-xtemp(I,J)).^2 + (Y(I,J)-ytemp(I,J)).^2));
    
    X(I,J) = xtemp(I,J);
    Y(I,J) = ytemp(I,J);
end
Xp = X;
Yp = Y;