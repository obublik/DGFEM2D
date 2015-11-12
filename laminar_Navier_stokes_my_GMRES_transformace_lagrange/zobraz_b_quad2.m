function zobraz_b_quad2

% zobrazuje vysledky
We = load('We.txt');
rad = load('rad_elementu.txt');
load 'geometrie';
PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3}+1;
nt = length(TP(:,1));

figure('color','w');
hold on;
kwe = 1;
M_max = 0;
for i = 1:nt
    rad2 = max(rad(i),2);
    
    nb = baze_poc(rad2);
    b = linspace(0,1,rad2);
    X = zeros(nb,1);
    Y = zeros(nb,1);
    s = 1;
    for j = 1:rad2
        for k = 1:rad2
            a4 = b(k)*(1-b(j));
            a1 = (1-b(k))*(1-b(j));
            a2 = (1-b(k))*b(j);
            a3 = b(k)*b(j);
            
            X(s) = a1*PX(TP(i,1)) + a2*PX(TP(i,2)) + a3*PX(TP(i,3)) + a4*PX(TP(i,4));
            Y(s) = a1*PY(TP(i,1)) + a2*PY(TP(i,2)) + a3*PY(TP(i,3)) + a4*PY(TP(i,4));
            s = s + 1;
        end
    end
    warning off;
    TP_loc = delaunay(X,Y);
    
    m = 1:rad(i)^2;
    W = We(kwe+m-1,:);
    if(rad == 1)
        W = [W;W;W;W];
    end
    
    kapa = 1.4;
    P = (kapa-1)*(W(:,4)-1./(2*W(:,1)).*(W(:,2).^2 + W(:,3).^2));
    M =((W(:,2).^2 + W(:,3).^2)./(kapa*P(:).*W(:,1))).^(1/2);
    trisurf(TP_loc,X,Y,M,'linestyle','none');
    
    kwe = kwe + nb;
    
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


function nb = baze_poc(rad)
nb = rad^2;

