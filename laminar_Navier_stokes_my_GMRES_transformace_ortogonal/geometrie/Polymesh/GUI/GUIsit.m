function GUIsit(action)
global P Q typ IH obd Xobd Yobd;
% graficke uzivatelske prostredi
if (nargin < 1)
    action ='initialize';
end;

if strcmp(action,'initialize')
    load P;
    load Q;
    
    figure(...
        'Name','Program pro generovani nestrukturovanych siti', ...
        'Position',[170 50 1000 600], ...
        'Color',[0.8,0.8,0.8]);
    
    axes(...
        'Units','normalized', ...
        'Position',[0.04 0.05 0.6 0.9], ...
        'Tag', 'fig');
  
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.84 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Tlacitka pro zadavani okrajovych podminek');
    
    % tlacitka pro zadavani okrajovych podminek
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.70 0.78 0.1 0.05], ...
        'visible', 'on', ...
        'String','Vyber body', ...
        'Callback','GUIsit(''vyberbody'')');

    uicontrol(...
        'Style','popupmenu', ...
        'Units','normalized', ...
        'Position',[0.8 0.78 0.1 0.05], ...
        'BackGroundColor','w', ...
        'String',{'stena','vstup','vystup','nevazka stena'}, ...
        'tag','druhZobraz');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.78 0.1 0.05], ...
        'visible', 'on', ...
        'String','Prirad typ hranice', ...
        'Callback','GUIsit(''priradtyp'')');
    
    
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.72 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Vypocteni geometrie');
    
    % tlacitko pro celkovy tisk, vyhlazeni a vypocteni geometrie
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.75 0.66 0.1 0.05], ...
        'visible', 'on', ...
        'String','Vyhlad sit', ...
        'Callback','GUIsit(''vyhladsit'')');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.66 0.1 0.05], ...
        'visible', 'on', ...
        'background','red', ...
        'String','Vypocti sit', ...
        'Callback','GUIsit(''vypoctisit'')');
    
    % Tlacitko pro konec
    uicontrol('Style','push', 'Units','normalized', 'Position',[0.9 0 0.1 0.05], 'String', 'Konec', 'Callback','GUIsit(''konec'')');
    
    % vykresluje sit
    maxX = max(P(:,1));
    minX = min(P(:,1));
    maxY = max(P(:,2));
    minY = min(P(:,2));
    axis([(minX-0.1) (maxX+0.1) (minY-0.1) (maxY-0.1)]);
    axis('equal');
    kresliSit;
    IH = najdiHranicniBody(P,Q); % vraci indexy bodu na hranici
    % vykresluje body hranice
    hold on;
    for i = IH
        plot(P(i,1),P(i,2),'.','Color','b');
    end
    hold off;
    
    % inicializacni funkce
    typ = zeros(length(P(:,1)),1);
    obd = [];
    
    
elseif strcmp(action,'vyberbody') % vybira body hranice
    if(isempty(obd) == 0)
        delete(obd);
    end
    [Xobd,Yobd] = ginput(2);
    hold on;
    obd = plot([Xobd(1),Xobd(2),Xobd(2),Xobd(1),Xobd(1)],[Yobd(1),Yobd(1),Yobd(2),Yobd(2),Yobd(1)],'r');
    hold off;

    
elseif strcmp(action,'priradtyp') % prirazuje zvoleny typ hranici
    if(Xobd(2) < Xobd(1))
        pom = Xobd(1);
        Xobd(1) = Xobd(2);
        Xobd(2) = pom;
    end
    if(Yobd(2) < Yobd(1))
        pom = Yobd(1);
        Yobd(1) = Yobd(2);
        Yobd(2) = pom;
    end
    IB = IH(P(IH,1) > Xobd(1) & P(IH,1) < Xobd(2) & P(IH,2) > Yobd(1) & P(IH,2) < Yobd(2));
    hran = get(findobj('tag','druhZobraz'), 'Value');
    hold on;
    switch hran
        case 1
            typ(IB) = 1;
            plot(P(IB,1),P(IB,2),'.','Color','k');
        case 2
            typ(IB) = -1;
            plot(P(IB,1),P(IB,2),'.','Color','g');
        case 3
            typ(IB) = -2;
            plot(P(IB,1),P(IB,2),'.','Color','r');
        case 4
            typ(IB) = 2;
            plot(P(IB,1),P(IB,2),'.','Color','m');
    end 
    hold off;

elseif strcmp(action,'zadejfunkci') % vypocitava vsechny hodnoty site potrebne pro vypocet
    open('fun.m');
 
elseif strcmp(action,'vypoctisit') % vypocitava vsechny hodnoty site potrebne pro vypocet
    % generovani stran
    E = [Q(:,[1,2]); Q(:,[2,3]); Q(:,[3,4]); Q(:,[4,1])];
    E = unique(sort(E,2),'rows');
    ne = length(E(:,1));
   
   % vypocet zakladnich geometrickych vztahu potrebnych pro vypocet
    P = [P,typ];
    
    Etyp = vypoctiGeometriiM(P,Q,E);
    
    cla;
    hold on;
    for i = 1:ne
        switch Etyp(i)
            case 0
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)]);
            case 1
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'k','linewidth',2);
            case (-1)
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'g','linewidth',2);
            case (-2)
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'r','linewidth',2);
            case (2)
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'m','linewidth',2);
        end
    end
    hold off;
    

elseif strcmp(action,'vyhlad sit') % vyhlazuje sit
    np = length(P(:,1));
    nq = length(Q(:,1));
    E = [Q(:,[1,2]); Q(:,[2,3]); Q(:,[3,4]); Q(:,[4,1])];
    E = unique(sort(E,2),'rows');
    ne = length(E(:,1));
    
    PPpoc = zeros(np,1);
    PPpom{np} = 0;
    for i = 1:ne
        PPpoc(E(i,1)) = PPpoc(E(i,1)) + 1;
        PPpoc(E(i,2)) = PPpoc(E(i,2)) + 1;  
        PPpom{E(i,1)} = [PPpom, E(i,2)];
        PPpom{E(i,2)} = [PPpom, E(i,1)];
    end
    PPpocmax = max(Ppoc);
    PP = zeros(np,PPpocmax);
    for i = 1:np
        for j = 1:PPpoc(i)
            PP(i,j) = PPpom{i}(j);
        end
    end
    
    % vlastni vyhlazeni
    I = 1:np;
    I = I(I ~= IH);
    for k = 1:3
        for i = I
            sx = 0;
            sy = 0;
            for j = 1:PPpoc(i)
                sx = sx + P(PP(i,j),1);
                sy = sy + P(PP(i,j),2);
            end
            sx = sx/PPpoc(i);
            sy = sy/PPpoc(i);
            P(i,1) = P(i,1)+sx;
            P(i,2) = P(i,2)+sy;
        end
    end
    kresliSit;
    
    
elseif strcmp(action,'konec') % ukoncuje program
    close(clf);
end
 
    
function Etyp = vypoctiGeometriiM(P,Q,E)
h = waitbar(0,'Please wait...');
PX = P(:,1);
PY = P(:,2);
Ptyp = P(:,3);

np = length(P(:,1));
nq = length(Q(:,1));
ne = length(E(:,1));

% matice obsahujici vsechny potrebne hodnoty pro vypocet
TP = zeros(nq,4);
TT = zeros(nq,4);
TE = zeros(nq,4);

% vytvoreni matice trojuhelniku/ctvercu
Qr = sparse(np,np);
for i = 1:nq
    Qr(Q(i,1),Q(i,2)) = i;
    Qr(Q(i,2),Q(i,3)) = i;
    Qr(Q(i,3),Q(i,4)) = i;
    Qr(Q(i,4),Q(i,1)) = i;
end

waitbar(0.25,h);

% vytvoreni matice stran
Er = sparse(np,np);
for i = 1:ne
    Er(E(i,1),E(i,2)) = i;
    Er(E(i,2),E(i,1)) = i;
end

waitbar(0.5,h);

% plneni matice Tg
for i = 1:nq
    TP(i,1) = Q(i,1)-1; % indexy bodu
    TP(i,2) = Q(i,2)-1;
    TP(i,3) = Q(i,3)-1;
    TP(i,4) = Q(i,4)-1;
    TT(i,1) = Qr(Q(i,2),Q(i,1))-1; % indexy sousednich ctvercu
    TT(i,2) = Qr(Q(i,3),Q(i,2))-1;
    TT(i,3) = Qr(Q(i,4),Q(i,3))-1;
    TT(i,4) = Qr(Q(i,1),Q(i,4))-1;
    TE(i,1) = Er(Q(i,1),Q(i,2)); % index prvni strany
    TE(i,2) = Er(Q(i,2),Q(i,3));
    TE(i,3) = Er(Q(i,3),Q(i,4));
    TE(i,4) = Er(Q(i,4),Q(i,1));
end

waitbar(0.75,h);
plus = [2,3,4,1];
Etyp = zeros(ne,1);
for i = 1:nq
    for j = 1:4
        jp = plus(j);
        if(TT(i,j) < 0)
            if(Ptyp(TP(i,j)+1) == -1 && Ptyp(TP(i,jp)+1) == -1)
                Etyp(TE(i,j)) = -1;
            elseif(Ptyp(TP(i,j)+1) == -2 && Ptyp(TP(i,jp)+1) == -2)
                Etyp(TE(i,j)) = -2;
            elseif(Ptyp(TP(i,j)+1) == 2 && Ptyp(TP(i,jp)+1) == 2)
                Etyp(TE(i,j)) = 2;
            else
                Etyp(TE(i,j)) = 1;
            end
        end
    end
end

for i = 1:nq
    for j = 1:4
        switch(Etyp(TE(i,j)))
            case 1 % vazka stena
                TT(i,j) = -1;
            case -1 % vstup
                TT(i,j) = -2;
            case -2 % vystup
                TT(i,j) = -3;
            case 2 % nevazka stena
                TT(i,j) = -4;
        end
    end
end

geometrie{1} = PX;
geometrie{2} = PY;
geometrie{3} = TP;
geometrie{4} = TT;
geometrie{5} = 4*ones(nq,1);

save 'DGFEM\geometrie' geometrie;

waitbar(1,h);
close(h);


function kresliSit
global P Q;
    hold on;
    for i = 1:length(Q(:,1))
        plot([P(Q(i,1),1), P(Q(i,2),1), P(Q(i,3),1), P(Q(i,4),1), P(Q(i,1),1)], [P(Q(i,1),2), P(Q(i,2),2), P(Q(i,3),2), P(Q(i,4),2), P(Q(i,1),2)],'g');
    end
    hold off;

function IH = najdiHranicniBody(P,Q)
m = length(P(:,1));
n = length(Q(:,1));
I = 1:m;
S = zeros(m,1);
for j = 1:n;
    v1 = [P(Q(j,2),1) - P(Q(j,1),1),P(Q(j,2),2) - P(Q(j,1),2)];
    v2 = [P(Q(j,3),1) - P(Q(j,2),1),P(Q(j,3),2) - P(Q(j,2),2)];
    v3 = [P(Q(j,4),1) - P(Q(j,3),1),P(Q(j,4),2) - P(Q(j,3),2)];
    v4 = [P(Q(j,1),1) - P(Q(j,4),1),P(Q(j,1),2) - P(Q(j,4),2)];
    vv1 = (v1(1)^2 + v1(2)^2)^(1/2);
    vv2 = (v2(1)^2 + v2(2)^2)^(1/2);
    vv3 = (v3(1)^2 + v3(2)^2)^(1/2);
    vv4 = (v4(1)^2 + v4(2)^2)^(1/2);
    S(Q(j,1)) = S(Q(j,1)) + abs(acos(-(v1(1)*v4(1) + v1(2)*v4(2))/(vv1*vv4)));
    S(Q(j,2)) = S(Q(j,2)) + abs(acos(-(v2(1)*v1(1) + v2(2)*v1(2))/(vv2*vv1)));
    S(Q(j,3)) = S(Q(j,3)) + abs(acos(-(v3(1)*v2(1) + v3(2)*v2(2))/(vv3*vv2)));
    S(Q(j,4)) = S(Q(j,4)) + abs(acos(-(v4(1)*v3(1) + v4(2)*v3(2))/(vv4*vv3)));
end

IH = I(S(I) <= (2*pi-0.1));











