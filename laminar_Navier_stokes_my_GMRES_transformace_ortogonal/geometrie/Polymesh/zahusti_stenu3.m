function zahusti_stenu3

load 'geometrie';
PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3} + 1;
TT = geometrie{4};
typ = geometrie{5};

nt = length(TT(:,1));
np = length(PX);

zahusti = zeros(np,1);
for i = 1:nt
    for j = 1:typ(i)
        if(TT(i,j) == -1)
            jp = mod(j,typ(i)) + 1;
            zahusti(TP(i,j)) = 1;
            zahusti(TP(i,jp)) = 1;
        end
    end
end

zt = zeros(nt,1);
for i = 1:nt
    if(typ(i) == 4)
        for j = 1:typ(i)
            if(zahusti(TP(i,j)) == 1)
               zt(i) = 1;
               break;
            end
        end
    end
end

nt2 = nt;
s = np + 1;
PP = sparse(2*np,2*np);
for i = 1:nt
    for j = 1:typ(i)
        PP(TP(i,j),TP(i,j)) = TP(i,j);
    end
    if(zt(i) == 1)
        for j = 1:typ(i)
            j1 = TP(i,j);
            j2 = TP(i,mod(j,typ(i)) + 1);
            if(PP(j1,j2) == 0)
                PX = [PX; (PX(j1)+PX(j2))/2];
                PY = [PY; (PY(j1)+PY(j2))/2];
                PP(j1,j2) = s;
                PP(j2,j1) = s;
                s = s + 1;
            end
        end
        j1 = TP(i,1);
        j2 = TP(i,2);
        j3 = TP(i,3);
        j4 = TP(i,4);
        n = typ(i);
        PX = [PX; sum(PX(TP(i,1:n)))/n];
        PY = [PY; sum(PY(TP(i,1:n)))/n];
        PP(j1,j3) = s;
        PP(j3,j1) = s;
        PP(j2,j4) = s;
        PP(j4,j2) = s;
        s = s + 1;
        nt2 = nt2 + 3;
    end
end

TP2 = zeros(nt2,6);
TT2 = zeros(nt2,6);
typ2 = zeros(nt2,1);
s = 1;
for i = 1:nt
    if(PP(TP(i,1),TP(i,3)) > 0 && typ(i) == 4)
        i1 = TP(i,1);
        i2 = TP(i,2);
        i3 = TP(i,3);
        i4 = TP(i,4);

        TP2(s,1:4) = [PP(i1,i1),PP(i1,i2),PP(i1,i3),PP(i1,i4)];
        TT2(s,1:4) = [TT(i,1),1,1,TT(i,4)];
        typ2(s) = 4;
        s = s + 1;
        TP2(s,1:4) = [PP(i1,i2),PP(i2,i2),PP(i2,i3),PP(i1,i3)];
        TT2(s,1:4) = [TT(i,1),TT(i,2),1,1];
        typ2(s) = 4;
        s = s + 1;
        TP2(s,1:4) = [PP(i1,i3),PP(i2,i3),PP(i3,i3),PP(i3,i4)];
        TT2(s,1:4) = [1,TT(i,2),TT(i,3),1];
        typ2(s) = 4;
        s = s + 1;
        TP2(s,1:4) = [PP(i1,i4),PP(i1,i3),PP(i3,i4),PP(i4,i4)];
        TT2(s,1:4) = [1,1,TT(i,3),TT(i,4)];
        typ2(s) = 4;
        s = s + 1;
    else
        v = zeros(2*typ(i),1);
        v2 = zeros(2*typ(i),1);
        p = 1;
        for j = 1:typ(i)
            j1 = TP(i,j);
            j2 = TP(i,mod(j,typ(i)) + 1);
            v(p) = PP(j1,j1);
            v2(p) = TT(i,j);
            p = p + 1;
            v(p) = PP(j1,j2);
            v2(p) = TT(i,j);
            p = p + 1;
        end
        
        I = 1:length(v);
        I = I(v ~= 0);
        TP2(s,1:length(I)) = v(I);
        TT2(s,1:length(I)) = v2(I);
        typ2(s) = length(I);
        s = s + 1;
    end
end

% hledani sousedu
np2 = length(PX);
Te = sparse(np2,np2);
for i = 1:nt2
    for j = 1:typ2(i)
        jp = mod(j,typ2(i)) + 1;
        if(TT2(i,j) > 0)
            Te(TP2(i,j),TP2(i,jp)) = i-1;
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
                plot(PX(TP2(i,[j,jp])), PY(TP2(i,[j,jp])),'color','k');
            case -2
                plot(PX(TP2(i,[j,jp])), PY(TP2(i,[j,jp])),'color','g');
            case -3
                plot(PX(TP2(i,[j,jp])), PY(TP2(i,[j,jp])),'color','r');
            otherwise
                plot(PX(TP2(i,[j,jp])), PY(TP2(i,[j,jp])));
         end
    end
end
axis equal
















