function vyhlad_sit
load geometrie

PX = geometrie{1};
PY = geometrie{2};
QP = geometrie{3};
QQ = geometrie{4};

figure
hold on
triplot(QP+1,PX,PY,'g')
axis equal

np = length(PX);
nt = length(QQ);
B = zeros(np,1);
j_plus = [2 3 1];
for i = 1:nt
    for j = 1:3
        jp = j_plus(j);
        if(QQ(i,j) < 0)
            B(QP(i,j)+1) = 1;
            B(QP(i,jp)+1) = 1;
        end
    end
end
for op = 1:2
    B1 = B;
    for i = 1:nt
        for j = 1:3
            jp = j_plus(j);
            if(B1(QP(i,j)+1) == 1 || B1(QP(i,jp)+1) == 1)
                B(QP(i,j)+1) = 1;
                B(QP(i,jp)+1) = 1;
            end
        end
    end
end
plot(PX(B==1),PY(B==1),'.r')

Ih = zeros(nt*3,1);
Jh = zeros(nt*3,1);
Hh = zeros(nt*3,1);
poc = zeros(np);
for i = 1:nt
    for j = 1:3
        jp = j_plus(j);
        poc(QP(i,jp)+1) = poc(QP(i,jp)+1) + 1;
    end
end
s = 1;
for i = 1:nt
    for j = 1:3
        jp = j_plus(j);
        Ih(s) = QP(i,jp) + 1;
        Jh(s) = QP(i,j) + 1;
        Hh(s) = 1/poc(QP(i,jp) + 1);
        s = s + 1;
    end
end
L = sparse(Ih,Jh,Hh,np,np);

I = 1:np;
I = I(B == 0);
for it = 1:10
    PXn = L*PX;
    PYn = L*PY;
    PX(I) = PXn(I);
    PY(I) = PYn(I);
end

triplot(QP+1,PX,PY,'k')
axis equal


% geometrie{1} = PX;
% geometrie{2} = PY;
% save 'geometrie' geometrie;