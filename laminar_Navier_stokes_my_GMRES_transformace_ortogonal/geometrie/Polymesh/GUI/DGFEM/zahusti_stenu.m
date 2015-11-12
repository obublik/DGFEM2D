function zahusti_stenu

load 'geometrie';
PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3} + 1;
TT = geometrie{4} + 1;
typ = geometrie{5};

nt = length(TT(:,1));
np = length(PX);

TP2 = zeros(2*nt,6);
TT2 = zeros(2*nt,6);
typ2 = zeros(2*nt,1);
for i = 1:nt
    TP2(i,1:typ(i)) = TP(i,1:typ(i));
    TT2(i,1:typ(i)) = TT(i,1:typ(i));
    typ2(i) = typ(i);
end

list = [];
sp = np + 1;
st = nt + 1;
PP = sparse(2*np,2*np);
for i = 1:nt
    if(typ(i) == 4)
        for j = 1:4
            if(TT(i,j) == 0)
                i1 = j;
                i2 = mod(j,4) + 1;
                i3 = mod(j+1,4) + 1;
                i4 = mod(j+2,4) + 1;

                j1 = TP(i,j);
                j2 = TP(i,mod(j,4) + 1);
                j3 = TP(i,mod(j+1,4) + 1);
                j4 = TP(i,mod(j+2,4) + 1);

                if(PP(j1,j4) == 0)
                    PX(sp) = (PX(j1) + PX(j4))/2;
                    PY(sp) = (PY(j1) + PY(j4))/2;
                    PP(j1,j4) = sp;
                    PP(j4,j1) = sp;
                    j5 = sp;
                    sp = sp + 1;
                else
                    j5 = PP(j1,j4);
                end
                if(PP(j2,j3) == 0)
                    PX(sp) = (PX(j2) + PX(j3))/2;
                    PY(sp) = (PY(j2) + PY(j3))/2;
                    PP(j2,j3) = sp;
                    PP(j3,j2) = sp;
                    j6 = sp;
                    sp = sp + 1;
                else
                    j6 = PP(j2,j3);
                end
                PX(sp) = (PX(j1) + PX(j2))/2;
                PY(sp) = (PY(j1) + PY(j2))/2;
                PP(j1,j2) = sp;
                PP(j1,j2) = sp;
                j7 = sp;
                sp = sp + 1;

                PX(sp) = (PX(j5) + PX(j6))/2;
                PY(sp) = (PY(j5) + PY(j6))/2;
                PP(j5,j6) = sp;
                PP(j6,j5) = sp;
                j8 = sp;
                sp = sp + 1;

                TP2(i,1:5) = [j5,j8,j6,j3,j4];
                TT2(i,1:5) = [st,st+1,TT(i,i2),TT(i,i3),TT(i,i4)];
                typ2(i) = 5;
                TP2(st,1:4) = [j1,j7,j8,j5];
                TT2(st,1:4) = [TT(i,i1),st+1,i,TT(i,i4)];
                typ2(st) = 4;
                st = st + 1;
                TP2(st,1:4) = [j7,j2,j6,j8];
                TT2(st,1:4) = [TT(i,i1),TT(i,i2),i,st-1];
                typ2(st) = 4;
                st = st + 1;

                list = [list,TT(i,i2),TT(i,i4)];
                break;
            end
        end
    end
end

% rohove bunky
list = unique(list);
for i = 1:length(list)
    for j = 1:4
        j1 = TP(i,j);
        j2 = TP(i,mod(j,4) + 1);
        if(PP(j1,j2) == 0)
            PX(sp) = (PX(j1) + PX(j2))/2;
            PY(sp) = (PY(j1) + PY(j2))/2;
            PP(j1,j2) = sp;
            PP(j2,j1) = sp;
            sp = sp + 1;
        end
    end
    xs = sum(PX(TP(i,:)))/4;
    ys = sum(PY(TP(i,:)))/4;
    PX(sp) = (PX(j1) + PX(j2))/2;
    PY(sp) = (PY(j1) + PY(j2))/2;
    PP(j1,j2) = sp;
    PP(j2,j1) = sp;
    sp = sp + 1;
end

for i = 1:length(list)
     j1 = TP(i,j);
     j2 = TP(i,mod(j,4) + 1);
     j3 = TP(i,mod(j+1,4) + 1);
     j4 = TP(i,mod(j+2,4) + 1);
    
end

nt2 = st-1;
TP2 = TP2(1:nt2,:);
TT2 = TT2(1:nt2,:);
typ2 = typ2(1:nt2);

% hledani sousedu
Te = sparse(nt2,nt2);
for i = 1:nt2
    for j = 1:typ2(i)
        jp = mod(j,typ2(i)) + 1;
        if(TT2(i,j) > 0)
            Te(TP2(i,j),TP2(i,jp)) = i;
        end
    end
end

for i = 1:nt2
    for j = 1:typ2(i)
        if(TT2(i,j) > 0)
            jp = mod(j,typ2(i)) + 1;
            TT2(i,j) = Te(TP2(i,jp),TP2(i,j));
        end
    end
end
TT2 = TT2-1;

geometrie{1} = PX;
geometrie{2} = PY;
geometrie{3} = TP2-1;
geometrie{4} = TT2;
geometrie{5} = typ2;

save 'geometrie' geometrie

figure
hold on
for i = 1:nt2
    for j = 1:typ2(i)
        jp = mod(j,typ2(i)) + 1;
        plot(PX(TP2(i,[j,jp])), PY(TP2(i,[j,jp])));
        switch TT2(i,j)
            case -1
                plot(PX(TP2(i,[j,jp])), PY(TP2(i,[j,jp])),'color','k','linewidth',2);
            case -2
                plot(PX(TP2(i,[j,jp])), PY(TP2(i,[j,jp])),'color','g','linewidth',2);
            case -3
                plot(PX(TP2(i,[j,jp])), PY(TP2(i,[j,jp])),'color','r','linewidth',2);
            otherwise
                plot(PX(TP2(i,[j,jp])), PY(TP2(i,[j,jp])));
         end
    end
end
axis equal

% % tisk
% poc = typ2;
% np = nt2;
% TP = TP2;
% TT = TT2;
% figure
% hold on
% for i = 1:np
%     for j = 1:poc(i)
%         jp = mod(j,poc(i))+1;
%         X = [PX(TP(i,j)),PX(TP(i,jp))];
%         Y = [PY(TP(i,j)),PY(TP(i,jp))];
%         switch(TT(i,j))
%             case -1
%                 plot(X,Y,'k');
%             case -2
%                 plot(X,Y,'g');
%             case -3
%                 plot(X,Y,'r');
%             case -4
%                 plot(X,Y,'m');
%             otherwise
%                 plot(X,Y,'b');
%         end
%     end
% end
% axis equal
% 
% I = 1:5:np;
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
% 
% for i = I
%     plot(PX(TP(i,1)),PY(TP(i,1)),'or');
%     for j = 1:poc(i)
%         text(PX(TP(i,j)),PY(TP(i,j)),num2str(j),'color','m');
%     end
% end
% axis equal
















