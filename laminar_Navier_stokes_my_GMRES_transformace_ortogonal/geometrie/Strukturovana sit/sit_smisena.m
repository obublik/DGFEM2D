function sit_smisena
n1 = 100;
n2 = 40;
L = 5;

X = zeros(n1,n2);
Y = zeros(n1,n2);

for i = 1:n1
    x = L*(i-1)/(n1-1);
    X(i,1) = x;
    X(i,n2) = x;
    
    Y(i,1) = fd(x);
    Y(i,n2) = fh(x);
end
    
for j = 1:n2
    X(1,j)  = 0;
    X(n1,j) = L;
    
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
typ = 4*ones(length(TP(:,1)),1);

nq = length(TT(:,1));
k = nq;
for i = 1:nq
    if(PY(TP(i,1)+1) > 0.5)
        if(TT(i,1) > -1)
            TT(TT(i,1)+1,3) = k;
        end
        if(TT(i,4) > -1)
            TT(TT(i,4)+1,2) = k;
        end
        TTpom = [TT(i,1),i-1,TT(i,4),0];
        TT(i,:) = [k,TT(i,2),TT(i,3),0];
        TT = [TT;TTpom];
        typ(i) = 3;
        typ(k+1) = 3;
        k = k + 1;
        
        TP = [TP; [TP(i,1),TP(i,2),TP(i,4),0]];
        TP(i,:) = [TP(i,4),TP(i,2),TP(i,3),0];
    end
end

geometrie{1} = PX;
geometrie{2} = PY;
geometrie{3} = TP;
geometrie{4} = TT;
geometrie{5} = typ;

save 'geometrie' geometrie

figure
hold on;
TP = TP + 1;
for i = 1:length(TP(:,1))
    if(typ(i) == 3) 
        plot([PX(TP(i,1)),PX(TP(i,2)),PX(TP(i,3)),PX(TP(i,1))],[PY(TP(i,1)),PY(TP(i,2)),PY(TP(i,3)),PY(TP(i,1))]);
%         text(sum(PX(TP(i,:)))/3,sum(PY(TP(i,:)))/3,num2str(i-1));
    else
        plot([PX(TP(i,1)),PX(TP(i,2)),PX(TP(i,3)),PX(TP(i,4)),PX(TP(i,1))],[PY(TP(i,1)),PY(TP(i,2)),PY(TP(i,3)),PY(TP(i,4)),PY(TP(i,1))]);
%         text(sum(PX(TP(i,:)))/4,sum(PY(TP(i,:)))/4,num2str(i-1));
    end
end
axis equal;

% k = 1;
% for i = 1:length(TP(:,1))
%     if(typ(i) == 3)
%         k = k + 1;
%     end
%     
%     if(typ(i) == 3 && k == 131) 
%         plot(PX(TP(i,1)),PY(TP(i,1)),'ob');
%         plot(PX(TP(i,2)),PY(TP(i,2)),'or');
%         plot(PX(TP(i,3)),PY(TP(i,3)),'og');
%         
%         plot(sum(PX(TP(TT(i,1)+1,:)))/3,sum(PY(TP(TT(i,1)+1,:)))/3,'.b');
%         plot(sum(PX(TP(TT(i,2)+1,:)))/3,sum(PY(TP(TT(i,2)+1,:)))/3,'.r');
%         plot(sum(PX(TP(TT(i,3)+1,:)))/3,sum(PY(TP(TT(i,3)+1,:)))/3,'.g');
%         break;
%     end
% end


function y = fd(x)
    h = 0.3;
    r = h/2 + 1/(8*h);
    ys = h-r;
    if(x < 1 || x > 2)
        y = 0;
    else
        y = sqrt(r^2-(x-1.5)^2) + ys;
    end
    
    
function y = fh(x)
    h = 0.3;
    r = h/2 + 1/(8*h);
    ys = 1 + r - h;
    if(x < 1 || x > 2)
        y = 1;
    else
        y = -sqrt(r^2-(x-1.5)^2) + ys;
    end
    
    
    