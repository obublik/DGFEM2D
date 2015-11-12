function A = GS_tri(n)
global integ
load integ;

A = eye(n,n);
E = eye(n,n);
for i = 2:n
    for j = 1:i-1
        c = -proj(A(j,:),E(i,:));
        for k = 1:n
            A(i,k) = A(i,k) + c*A(j,k);
        end
    end
    A(i,:) = A(i,:)/norma(A(i,:));
end

fid = fopen('A_baze_tri.txt','w');
for i = 1:n
    for j = 1:n
        fprintf(fid,'%6.16f ',A(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

% K = zeros(n,n);
% for i = 2:n
%     for j = 1:i
%         K(i,j) = soucin(A(i,:),A(j,:));
%     end
% end
% K


function p = proj(a,b)
p = soucin(a,b)/soucin(a,a);

function n = norma(a)
n = sqrt(soucin(a,a));

function s = soucin(a,b)
global integ
m = length(a);
s = 0;
for i = 1:integ.n
    x = integ.b(i,1);
    y = integ.b(i,2);
    fa = 0;
    fb = 0;
    for j = 1:m
        fa = fa + a(j)*x^rx(j)*y^ry(j);
        fb = fb + b(j)*x^rx(j)*y^ry(j);
    end
    s = s + integ.w(i)*fa*fb/2;
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
