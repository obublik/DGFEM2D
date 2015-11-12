function dual_mesh

load 'geometrie'

PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3}+1;
TT = geometrie{4};

nt = length(TT(:,1));
np = length(PX);
Tdual{np} = [];
for i = 1:np
    Tdual{i} = zeros(15,8);
end
poc = zeros(np,1);
for i = 1:nt
    xs = sum(PX(TP(i,:)))/3;
    ys = sum(PY(TP(i,:)))/3;
    for j = 1:3
        jp = mod(j,3)+1;
        if(TT(i,j) > -1)
            xs1 = sum(PX(TP(TT(i,j)+1,:)))/3;
            ys1 = sum(PY(TP(TT(i,j)+1,:)))/3;
            poc(TP(i,jp)) = poc(TP(i,jp)) + 1;
            Tdual{TP(i,jp)}(poc(TP(i,jp)),:) = [TP(i,j),xs1,ys1,xs,ys,1,TT(i,j)+1,0];
        else
            xs1 = sum(PX(TP(i,[j,jp])))/2;
            ys1 = sum(PY(TP(i,[j,jp])))/2;
            
            poc(TP(i,j)) = poc(TP(i,j)) + 1;
            Tdual{TP(i,j)}(poc(TP(i,j)),:) =  [TP(i,jp),xs,ys,xs1,ys1,2,i,0];
            poc(TP(i,j)) = poc(TP(i,j)) + 1;
            Tdual{TP(i,j)}(poc(TP(i,j)),:) =  [TT(i,j),xs1,ys1,PX(TP(i,j)),PY(TP(i,j)),3,TP(i,j),TP(i,jp)];
            
            poc(TP(i,jp)) = poc(TP(i,jp)) + 1;
            Tdual{TP(i,jp)}(poc(TP(i,jp)),:) =  [TT(i,j),PX(TP(i,jp)),PY(TP(i,jp)),xs1,ys1,4,TP(i,jp),TP(i,jp)];
            poc(TP(i,jp)) = poc(TP(i,jp)) + 1;
            Tdual{TP(i,jp)}(poc(TP(i,jp)),:) = [TP(i,j),xs1,ys1,xs,ys,5,TP(i,j),TP(i,jp)];
        end
    end
end

% serazeni
for i = 1:np
    xs = sum(Tdual{i}(:,2))/poc(i);
    ys = sum(Tdual{i}(:,3))/poc(i);
    for j = 1:poc(i)-1
        for k = 1:poc(i)-1
            alfa = uhel(Tdual{i}(k,2)-xs, Tdual{i}(k,3)-ys);
            beta = uhel(Tdual{i}(k+1,2)-xs, Tdual{i}(k+1,3)-ys);
            if(alfa > beta)
                pom = Tdual{i}(k,:);
                Tdual{i}(k,:) = Tdual{i}(k+1,:);
                Tdual{i}(k+1,:) = pom;
            end
        end
    end
end

% vytvoreni bodu
PP = sparse(np,np);
s = nt+1;
for i = 1:nt
    for j = 1:3
        jp = mod(j,3)+1;
        if(TT(i,j) < 0)
            PP(TP(i,j),TP(i,j)) = s;
            s = s + 1;
            PP(TP(i,j),TP(i,jp)) = s;
            PP(TP(i,jp),TP(i,j)) = s;
            s = s + 1;
        end
    end
end
nb = s-1;
PX = zeros(nb,1);
PY = zeros(nb,1);
TP = zeros(np,12);
for i = 1:np
    for j = 1:poc(i)
        switch Tdual{i}(j,6);
            case 1
                TP(i,j) = Tdual{i}(j,7);
                PX(TP(i,j)) = Tdual{i}(j,2);
                PY(TP(i,j)) = Tdual{i}(j,3);
            case 2
                TP(i,j) = Tdual{i}(j,7);
                PX(TP(i,j)) = Tdual{i}(j,2);
                PY(TP(i,j)) = Tdual{i}(j,3);
            otherwise
                TP(i,j) = PP(Tdual{i}(j,7),Tdual{i}(j,8));
                PX(TP(i,j)) = Tdual{i}(j,2);
                PY(TP(i,j)) = Tdual{i}(j,3);
        end
    end
end

% vytvoreni matice TT
TT = zeros(np,12);
for i = 1:np
    for j = 1:poc(i)
        t = Tdual{i}(j,1);
        if(t > 0)
            TT(i,j) = t-1;
        else
            TT(i,j) = t;
        end
    end
end

% tisk
figure
hold on
for i = 1:np
    for j = 1:poc(i)
        jp = mod(j,poc(i))+1;
        X = [PX(TP(i,j)),PX(TP(i,jp))];
        Y = [PY(TP(i,j)),PY(TP(i,jp))];
        switch(TT(i,j))
            case -1
                plot(X,Y,'k');
            case -2
                plot(X,Y,'g');
            case -3
                plot(X,Y,'r');
            case -4
                plot(X,Y,'m');
            otherwise
                plot(X,Y,'b');
        end
    end
end
axis equal

I = 1:5:np;
% for i = I
%     xs = 0;
%     ys = 0;
%     for j = 1:poc(i)
%         xs = xs + PX(TP(i,j));
%         ys = ys + PY(TP(i,j));
%     end
%     xs = xs/poc(i);
%     ys = ys/poc(i);
%     text(xs,ys,num2str(i-1),'color','r');
%     for j = 1:poc(i)
%         jp = mod(j,poc(i)) + 1;
%         xs = 0.75*PX(TP(i,j)) + 0.25*PX(TP(i,jp));
%         ys = 0.75*PY(TP(i,j)) + 0.25*PY(TP(i,jp));
%         text(xs,ys,num2str(TT(i,j)));
%     end
% end

for i = I
    plot(PX(TP(i,1)),PY(TP(i,1)),'or');
    for j = 1:poc(i)
        text(PX(TP(i,j)),PY(TP(i,j)),num2str(j),'color','m');
    end
end
axis equal

geometrie{1} = PX;
geometrie{2} = PY;
geometrie{3} = TP-1;
geometrie{4} = TT;
geometrie{5} = poc;

save 'geometrie_dual' geometrie


function alfa = uhel(b1,b2)
alfa = acos(b2/sqrt(b1^2+b2^2));
if(b1 < 0)
    alfa = 2*pi - alfa;
end



