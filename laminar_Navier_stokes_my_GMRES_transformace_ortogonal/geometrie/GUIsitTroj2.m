function GUIsitTroj2(action)
global geom node cnect prvek P T typ IH obd Xobd Yobd;
warning off;

% graficke uzivatelske prostredi
if (nargin < 1)
    action ='initialize';
end;

if strcmp(action,'initialize')
    clc;
    load 'geom';
    figure(...
        'Name','Program pro generovani nestrukturovanych siti', ...
        'Position',[170 150 1000 600], ...
        'Color',[0.8,0.8,0.8]);
    
    axes(...
        'Units','normalized', ...
        'Position',[0.04 0.05 0.6 0.9], ...
        'Tag', 'fig');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.95 0.4 0.04],'visible', 'on', 'background','blue', 'foregroundcolor', 'white','String', 'Tlacitka pro zadavani geometrie');
    
    % tlacitka pro tvorbu geometrie
    % obdelnik
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.90 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','xld');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.67 0.90 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','xld');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.71 0.90 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','yld');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.73 0.90 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','yld');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.77 0.90 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','d');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.79 0.90 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','d');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.83 0.90 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','v');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.85 0.90 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','v');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.95 0.90 0.05 0.04], ...
        'visible', 'on', ...
        'String','obdelnik', ...
        'Callback','GUIsitTroj2(''obdelnik'')');
    
    % elipsa
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','sx');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.67 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','sx');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.70 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','sy');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.72 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','sy');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','a');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.77 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','a');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.80 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','b');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.82 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','b');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.85 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','uhel');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.87 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','uhel');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.90 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','n');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.92 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','10','tag','n');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.95 0.85 0.05 0.04], ...
        'visible', 'on', ...
        'String','elipsa', ...
        'Callback','GUIsitTroj2(''elipsa'')');
    
    % nacteni geometrie
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.80 0.12 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','Zadej jmeno souboru');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.77 0.80 0.06 0.04],'BackGroundColor','w','visible', 'on', 'String','...','tag','soub');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.83 0.80 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','px');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.85 0.80 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','px');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.89 0.80 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','py');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.91 0.80 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','py');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.95 0.80 0.05 0.04], ...
        'visible', 'on', ...
        'String','Nacti', ...
        'Callback','GUIsitTroj2(''nacti'')');
    
    % tlacitko pro prijmuti prvku
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.92 0.73 0.08 0.05], ...
        'visible', 'on', ...
        'String','Prijmi prvek', ...
        'Callback','GUIsitTroj2(''prijmiPrvek'')');
    
    % tlacitko pro smazani geometrie
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.84 0.73 0.08 0.05], ...
        'visible', 'on', ...
        'String','Smaz geometrii', ...
        'Callback','GUIsitTroj2(''smazgeometrii'')');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.68 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Zahusteni site');
    
    % tlacitka pro zahusteni site
    uicontrol(...
        'Style','radiobutton', ...
        'Units','normalized', ...
        'Position',[0.78 0.62 0.07 0.05], ...
        'BackGroundColor',[0.8,0.8,0.8], ...
        'visible', 'on', ...
        'HitTest', 'off', ...
        'String','Zahustit?', ...
        'Tag','zahustit');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.85 0.62 0.1 0.05], ...
        'visible', 'on', ...
        'String','Zadej funkci', ...
        'Callback','GUIsitTroj2(''zadejfunkci'')');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.57 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Tlacitka pro vytvoreni site');
    
    % tlacitko pro generovani site
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.76 0.51 0.04 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','hmax');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.8 0.51 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0.1','tag','hmax');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.85 0.51 0.1 0.05], ...
        'visible', 'on', ...
        'String','Generuj sit', ...
        'Callback','GUIsitTroj2(''generujsit'')');
    
    % tlacitka pro nacteni site
    uicontrol(...
        'Style','radiobutton', ...
        'Units','normalized', ...
        'Position',[0.7 0.51 0.06 0.05], ...
        'BackGroundColor',[0.8,0.8,0.8], ...
        'visible', 'on', ...
        'HitTest', 'off', ...
        'String','Nacist?', ...
        'Tag','nacist');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.46 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Tlacitka pro zadavani okrajovych podminek');
    
    % tlacitka pro zadavani okrajovych podminek
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.70 0.4 0.1 0.05], ...
        'visible', 'on', ...
        'String','Vyber body', ...
        'Callback','GUIsitTroj2(''vyberbody'')');

    uicontrol(...
        'Style','popupmenu', ...
        'Units','normalized', ...
        'Position',[0.8 0.4 0.1 0.05], ...
        'BackGroundColor','w', ...
        'String',{'vstup','vystup','stena','nevazka stena'}, ...
        'tag','druhZobraz');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.4 0.1 0.05], ...
        'visible', 'on', ...
        'String','Prirad typ hranice', ...
        'Callback','GUIsitTroj2(''priradtyp'')');
    
    
    %====================================
    % tlacitka pro upravu site
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.35 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Uprava site');
    
    % tlacitko pro zahusteni u steny
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.75 0.29 0.1 0.05], ...
        'visible', 'on', ...
        'String','Zahusti u steny', ...
        'Callback','GUIsitTroj2(''zahusti_u_steny'')');
    
    % tlacitko pro vyhlazeni site
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.29 0.1 0.05], ...
        'visible', 'on', ...
        'String','Vyhlad sit', ...
        'Callback','GUIsitTroj2(''vyhladsit'')');
    
    %====================================
    % tlacitka pro celkovy tisk a vypocteni geometrie
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.24 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Vypocteni geometrie');
    
    % textova pole pro informace o siti
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.17 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'Pocet trojuhelniku:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.17 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'pocTroj');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.12 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'Pocet bodu:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.12 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'pocBod');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.07 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'Prumerna kvalita:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.07 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'prumKval');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.02 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'Minimalni kvalita:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.02 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'minKval');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.18 0.1 0.05], ...
        'visible', 'on', ...
        'background','red', ...
        'String','Vypocti sit', ...
        'Callback','GUIsitTroj2(''vypoctisit'')');
    
    % Tlacitko pro konec
    uicontrol('Style','push', 'Units','normalized', 'Position',[0.9 0 0.1 0.05], 'String', 'Konec', 'Callback','GUIsitTroj2(''konec'')');
    
    % inicializacni funkce
    zobrazHranici;
    obd = [];

%--------------------------------------------------------------------------    
elseif strcmp(action,'obdelnik') % pridava obdelnik do dane geometrie
    zobrazHranici;
    x = str2double(get(findobj('tag','xld'),'string'));
    y = str2double(get(findobj('tag','yld'),'string'));
    d = str2double(get(findobj('tag','d'),'string'));
    v = str2double(get(findobj('tag','v'),'string'));
    
    prvek = [x,y; x+d,y; x+d, y+v; x, y+v];
    tiskni(prvek);
 
%--------------------------------------------------------------------------     
elseif strcmp(action,'elipsa') % pridava elipsu do dane geometrie
    zobrazHranici;
    sx = str2double(get(findobj('tag','sx'),'string'));
    sy = str2double(get(findobj('tag','sy'),'string'));
    a = str2double(get(findobj('tag','a'),'string'));
    b = str2double(get(findobj('tag','b'),'string'));
    alfa = str2double(get(findobj('tag','uhel'),'string'));
    n = str2double(get(findobj('tag','n'),'string'));
    
    I = 1:n;
    Xp = a*cos(2*pi*I/n)';
    Yp = b*sin(2*pi*I/n)';
    X = cos(alfa)*Xp + sin(alfa)*Yp + sx;
    Y = -sin(alfa)*Xp + cos(alfa)*Yp + sy;
    prvek = [X,Y];
    tiskni(prvek);
    
%--------------------------------------------------------------------------     
elseif strcmp(action,'nacti') % nacita libovolnou polygonialni geometrii
    zobrazHranici;
    soub = get(findobj('tag','soub'),'string');
    px = str2double(get(findobj('tag','px'),'string'));
    py = str2double(get(findobj('tag','py'),'string'));
    eval(['load ',soub]);
    
    prvek = [data(:,1)+px, data(:,2)+py];
    tiskni(prvek);
    
%--------------------------------------------------------------------------     
elseif strcmp(action,'prijmiPrvek') % prijima prvek geometrie
    if(isempty(prvek) ~= 1) 
        if(isempty(geom) == 1)    
            geom{1} = prvek;
        else
            geom{length(geom)+1} = prvek;
        end
    end
   prvek = [];
   zobrazHranici;
    
%--------------------------------------------------------------------------    
elseif strcmp(action,'smazgeometrii') % maze geometrii
    geom = [];
    obd = [];
    cla;
   
%--------------------------------------------------------------------------     
elseif strcmp(action,'generujsit') % generujesit
    obd = [];
    n = length(geom);
    
    % generovani vektoru uzlu
    node = [];
    for i = 1:n
        node = [node;geom{i}];
    end
    
    % generovani vektoru propojeni
    cnect = [];
    for i = 1:n
        n1 = length(geom{i}(:,1));
        n2 = length(cnect);
        cnect = [cnect; (1:n1-1)' + n2,(2:n1)' + n2;  n1 + n2, n2+1];
    end
    
    hdata.hmax = str2double(get(findobj('tag','hmax'),'string'));
    zahustit = get(findobj('tag','zahustit'),'Value');
    if(zahustit == 1)
        hdata.fun = @fun;
    end
    
    nacti = get(findobj('tag','nacist'),'Value');
    if(nacti == 0)
        [P,T,junk,stat] = meshfaces(node,cnect,[],hdata,[]); % vypujceny kod
        
        % nastaveni statistiky
        set(findobj('tag','pocTroj'),'String',num2str(stat.Triangles));
        set(findobj('tag','pocBod'),'String',num2str(stat.Nodes));
        set(findobj('tag','prumKval'),'String',num2str(stat.Mean_quality));
        set(findobj('tag','minKval'),'String',num2str(stat.Min_quality));
    else
        load mesh;
        P = mesh{1};
        T = mesh{2};
        
        % nastaveni statistiky
        set(findobj('tag','pocTroj'),'String',num2str(length(T(:,1))));
        set(findobj('tag','pocBod'),'String',num2str(length(P(:,1))));
        set(findobj('tag','prumKval'),'String',num2str(0));
        set(findobj('tag','minKval'),'String',num2str(0));
    end

    typ = zeros(length(P(:,1)),1);
    %P je matice, ktera ma v prvnim radku souradnice bodu X, ve druhem sloupci je souradnice Y, ve
    %tretim sloupci je okajova podminka bodu

    % urceni vnitrnich bodu
    nt = length(T(:,1));
    np = length(P(:,1));
    uhel = zeros(nt,3);
    I = 1:nt;
    uhel(I,1) = abs(acos(((P(T(I,2),1)-P(T(I,1),1)).*(P(T(I,3),1)-P(T(I,1),1)) + (P(T(I,2),2)-P(T(I,1),2)).*(P(T(I,3),2)-P(T(I,1),2)))./(((P(T(I,2),1)-P(T(I,1),1)).^2 + (P(T(I,2),2)-P(T(I,1),2)).^2).^(1/2).*((P(T(I,3),1)-P(T(I,1),1)).^2 + (P(T(I,3),2)-P(T(I,1),2)).^2).^(1/2))));
    uhel(I,2) = abs(acos(((P(T(I,3),1)-P(T(I,2),1)).*(P(T(I,1),1)-P(T(I,2),1)) + (P(T(I,3),2)-P(T(I,2),2)).*(P(T(I,1),2)-P(T(I,2),2)))./(((P(T(I,3),1)-P(T(I,2),1)).^2 + (P(T(I,3),2)-P(T(I,2),2)).^2).^(1/2).*((P(T(I,1),1)-P(T(I,2),1)).^2 + (P(T(I,1),2)-P(T(I,2),2)).^2).^(1/2))));
    uhel(I,3) = abs(acos(((P(T(I,1),1)-P(T(I,3),1)).*(P(T(I,2),1)-P(T(I,3),1)) + (P(T(I,1),2)-P(T(I,3),2)).*(P(T(I,2),2)-P(T(I,3),2)))./(((P(T(I,1),1)-P(T(I,3),1)).^2 + (P(T(I,1),2)-P(T(I,3),2)).^2).^(1/2).*((P(T(I,2),1)-P(T(I,3),1)).^2 + (P(T(I,2),2)-P(T(I,3),2)).^2).^(1/2))));

    % uhel prirazeny jednotlivym bodum
    up = zeros(np,1);
    for i = 1:nt
        up(T(i,1)) = up(T(i,1)) + uhel(i,1);
        up(T(i,2)) = up(T(i,2)) + uhel(i,2);
        up(T(i,3)) = up(T(i,3)) + uhel(i,3);
    end

    % IH jsou indexy bodu na hranici
    I = 1:length(P(:,1));
    eps = 0.1;
    IH = I(up(I) < 2*pi-eps);
    
    % vykresleni site
    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'g');
    plot(P(IH,1),P(IH,2),'.','Color','b');
    xmax = max(P(:,1));
    xmin = min(P(:,1));
    ymax = max(P(:,2));
    ymin = min(P(:,2));
    axis([xmin-0.3 xmax+0.3 ymin-0.3 ymax+0.3]);
    axis('equal');
    hold off;

%--------------------------------------------------------------------------  
elseif strcmp(action,'vyhladsit') % vyhlazuje sit
    % generovani stran
    P = vyhlad_sit(T,P,typ);
    
    % vykresleni site
    I = 1:length(P(:,1));
    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'g');
    IB = I(typ == -1);
    plot(P(IB,1),P(IB,2),'.','Color','g');
    IB = I(typ == -2);
    plot(P(IB,1),P(IB,2),'.','Color','r');
    IB = I(typ == 1);
    plot(P(IB,1),P(IB,2),'.','Color','k');
    IB = I(typ == 2);
    plot(P(IB,1),P(IB,2),'.','Color','m');
    xmax = max(P(:,1));
    xmin = min(P(:,1));
    ymax = max(P(:,2));
    ymin = min(P(:,2));
    axis([xmin-0.3 xmax+0.3 ymin-0.3 ymax+0.3]);
    axis('equal');
    hold off;
    
%--------------------------------------------------------------------------     
elseif strcmp(action,'vyberbody') % vybira body hranice
    if(isempty(obd) == 0)
        delete(obd);
    end
    [Xobd,Yobd] = ginput(2);
    hold on;
    obd = plot([Xobd(1),Xobd(2),Xobd(2),Xobd(1),Xobd(1)],[Yobd(1),Yobd(1),Yobd(2),Yobd(2),Yobd(1)],'r');
    hold off;

%--------------------------------------------------------------------------     
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
            typ(IB) = -1;
            plot(P(IB,1),P(IB,2),'.','Color','g');
        case 2
            typ(IB) = -2;
            plot(P(IB,1),P(IB,2),'.','Color','r');
       case 3
            typ(IB) = 1;
            plot(P(IB,1),P(IB,2),'.','Color','k');
       case 4
            typ(IB) = 2;
            plot(P(IB,1),P(IB,2),'.','Color','m');
    end 
    hold off;

%--------------------------------------------------------------------------     
elseif strcmp(action,'zadejfunkci') % otevira m-file pro zadavani funkce pro zahusteni
    open('fun.m');
    
%--------------------------------------------------------------------------     
elseif strcmp(action,'zahusti_u_steny')  % zahusteni u sten
    [T,P,typ] = adaptuj(T,P,typ);
    
    % vykresleni site
    I = 1:length(P(:,1));
    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'g');
    IB = I(typ == -1);
    plot(P(IB,1),P(IB,2),'.','Color','g');
    IB = I(typ == -2);
    plot(P(IB,1),P(IB,2),'.','Color','r');
    IB = I(typ == 1);
    plot(P(IB,1),P(IB,2),'.','Color','k');
    IB = I(typ == 2);
    plot(P(IB,1),P(IB,2),'.','Color','m');
    xmax = max(P(:,1));
    xmin = min(P(:,1));
    ymax = max(P(:,2));
    ymin = min(P(:,2));
    axis([xmin-0.3 xmax+0.3 ymin-0.3 ymax+0.3]);
    axis('equal');
    hold off;
    
    % nastaveni statistiky
    set(findobj('tag','pocTroj'),'String',num2str(length(T(:,1))));
    set(findobj('tag','pocBod'),'String',num2str(length(P(:,1))));
    
%--------------------------------------------------------------------------  
elseif strcmp(action,'vypoctisit') % vypocitava vsechny hodnoty site potrebne pro vypocet
    % generovani stran
    E = [T(:,[1,2]); T(:,[2,3]); T(:,[3,1])];
    E = unique(sort(E,2),'rows');
    ne = length(E(:,1));
    typE = zeros(ne,1);
    np = length(P(:,1));
    
    % vypocet zakladnich geometrickych vztahu potrebnych pro vypocet
    nq = length(T(:,1));
    Q = zeros(nq,3);
    I = 1:nq;
    J = 1:3;
    Q(I,J) = T(I,J);
    
    P = [P,typ];
    
    Etyp = vypoctiGeometrii_impl(P,Q,E);
    vypoctiGeometrii_expl(P,Q,[E,Etyp]);

    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'b');
    for i = 1:ne
        switch Etyp(i)
            case 1
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'k','linewidth',2);
            case 2
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'m','linewidth',2);
            case (-1)
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'g','linewidth',2);
            case (-2)
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'r','linewidth',2);
        end
    end
    hold off;
    
%--------------------------------------------------------------------------     
elseif strcmp(action,'konec') % ukoncuje program
    save 'geom' geom;
    close(clf);
end


% _________________________________________________________________________
function zobrazHranici
    global geom;
    cla;
    if(isempty(geom) == 0)
        n = length(geom);
        xmax = -1000;
        xmin = 1000;
        ymax = -1000;
        ymin = 1000;
        hold on;
        for i = 1:n
            plot([geom{i}(:,1);geom{i}(1,1)],[geom{i}(:,2);geom{i}(1,2)],'k','linewidth',2);
            pxmax = max(geom{i}(:,1));
            pxmin = min(geom{i}(:,1));
            pymax = max(geom{i}(:,2));
            pymin = min(geom{i}(:,2));
            if(pxmax > xmax)
                xmax = pxmax;
            end
            if(pxmin < xmin)
                xmin = pxmin;
            end
            if(pymax > ymax)
                ymax = pymax;
            end
            if(pymin < ymin)
                ymin = pymin;
            end
        end
        axis([xmin-0.3 xmax+0.3 ymin-0.3 ymax+0.3]);
        axis('equal');
        hold off;
    else
        cla;
    end
 
    
% _________________________________________________________________________    
function tiskni(prvek)
    prvek = [prvek; prvek(1,1) prvek(1,2)];
    hold on;
    plot(prvek(:,1),prvek(:,2),'r');
    hold off;
    
    
function Etyp = vypoctiGeometrii_impl(P,Q,E)
h = waitbar(0,'Please wait...');
X = P(:,1);
Y = P(:,2);
Ptyp = P(:,3);

np = length(P(:,1));
nq = length(Q(:,1));
ne = length(E(:,1));

% matice obsahujici vsechny potrebne hodnoty pro vypocet
QP = zeros(nq,3);
QQ = zeros(nq,3);
QE = zeros(nq,3);
QSx = zeros(nq,3);
QSy = zeros(nq,3);
QTx = zeros(nq,1);
QTy = zeros(nq,1);
QO = zeros(nq,1);
Qpricka = zeros(nq,1);
QLR = zeros(nq,3);

EP = zeros(ne,2);
EQ = zeros(ne,2);
Enx = zeros(ne,1);
Eny = zeros(ne,1);
Etyp = zeros(ne,1);

% vytvoreni matice trojuhelniku/ctvercu
Qr = sparse(np,np);
for i = 1:nq
    Qr(Q(i,1),Q(i,2)) = i;
    Qr(Q(i,2),Q(i,3)) = i;
    Qr(Q(i,3),Q(i,1)) = i;
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
    QP(i,1) = Q(i,1); % indexy bodu
    QP(i,2) = Q(i,2);
    QP(i,3) = Q(i,3);
    QQ(i,1) = Qr(Q(i,2),Q(i,1)); % indexy sousednich ctvercu
    QQ(i,2) = Qr(Q(i,3),Q(i,2));
    QQ(i,3) = Qr(Q(i,1),Q(i,3));
    QE(i,1) = Er(Q(i,1),Q(i,2)); % index prvni strany
    QE(i,2) = Er(Q(i,2),Q(i,3));
    QE(i,3) = Er(Q(i,3),Q(i,1));
    QSx(i,1) = Y(Q(i,2)) - Y(Q(i,1)); % Sx1
    QSx(i,2) = Y(Q(i,3)) - Y(Q(i,2));
    QSx(i,3) = Y(Q(i,1)) - Y(Q(i,3));
    QSy(i,1) = -(X(Q(i,2)) - X(Q(i,1))); % Sy1
    QSy(i,2) = -(X(Q(i,3)) - X(Q(i,2)));
    QSy(i,3) = -(X(Q(i,1)) - X(Q(i,3)));
    QTx(i) = (X(Q(i,1)) + X(Q(i,2)) + X(Q(i,3)))/3; % x-ova souradnice teziste
    QTy(i) = (Y(Q(i,1)) + Y(Q(i,2)) + Y(Q(i,3)))/3;
    QO(i) = abs(sum(cross([X(Q(i,2))-X(Q(i,1)),Y(Q(i,2))-Y(Q(i,1)),0],[X(Q(i,3))-X(Q(i,1)),Y(Q(i,3))-Y(Q(i,1)),0])))/2; % obsah trojuhelniku
end

waitbar(0.75,h);

I = 1:ne;
vel(I) = ((X(E(I,2))-X(E(I,1))).^2 + (Y(E(I,2))-Y(E(I,1))).^2).^(1/2); % velikost strany
for i = 1:ne;
    EP(i,1) = E(i,1); % index bodu strany
    EP(i,2) = E(i,2);
    EQ(i,1) = Qr(E(i,1),E(i,2)); % index prvniho sousedniho ctverce
    EQ(i,2) = Qr(E(i,2),E(i,1));
    Enx(i) = (Y(E(i,2))-Y(E(i,1)))/vel(i); % nx
    Eny(i) =  -(X(E(i,2))-X(E(i,1)))/vel(i); % ny
end

% pokud neni strana vnitrni, sousedni ctverec se presune na prvni misto
for i = 1:ne;
    if(EQ(i,1) == 0 || EQ(i,2) == 0)
        if(EQ(i,1) == 0)
            EQ(i,1) = EQ(i,2);
            EQ(i,2) = 0;
        end
    end
end

% typ strany
for i = 1:ne;
    if(EQ(i,2) == 0)
        if((Ptyp(EP(i,1)) == 2 && Ptyp(EP(i,2)) == 2) || (Ptyp(EP(i,1)) == 2 && Ptyp(EP(i,2)) < 0) || (Ptyp(EP(i,1)) < 0 && Ptyp(EP(i,2)) == 2))
            Etyp(i) = 2;
        elseif((Ptyp(EP(i,1)) == -1 && Ptyp(EP(i,2)) == -1) || (Ptyp(EP(i,1)) == -1 && Ptyp(EP(i,2)) == -2) || (Ptyp(EP(i,1)) == -2 && Ptyp(EP(i,2)) == -1))
            Etyp(i) = -1;
        elseif(Ptyp(EP(i,1)) == -2 && Ptyp(EP(i,2)) == -2)
            Etyp(i) = -2;
        else
            Etyp(i) = 1;
        end
    end
end

% Nasledujici tri indexy prirazuji stranam ve ctverci cislo 1 nebo 2
% podle toho, zdali je ctverec pro danou stranu levy nebo pravy. Toho
% bude vyuzito pri linearni rekonstrukci
for i = 1:nq
    if(i == EQ(QE(i,1),1))
        QLR(i,1) = 1;
    else
        QLR(i,1) = 2;
    end
    
    if(i == EQ(QE(i,2),1))
        QLR(i,2) = 1;
    else
        QLR(i,2) = 2;
    end
    
    if(i == EQ(QE(i,3),1))
        QLR(i,3) = 1;
    else
        QLR(i,3) = 2;
    end
end

% vypocteni nejmensi strany ve ctverci
QS = sqrt(QSx.^2 + QSy.^2);
for i = 1:nq
      Qpricka(i) = min([QS(i,1),QS(i,2),QS(i,3)]);
end


% provazani bodu a ctvercu
PQpoc = zeros(np,1);
PQpom{np} = [];
for i = 1:nq
    for j = 1:3
        PQpoc(QP(i,j)) = PQpoc(QP(i,j)) + 1;
        PQpom{QP(i,j)} = [PQpom{QP(i,j)},i];
    end
end

Pmax = max(PQpoc);
PQ = zeros(np,Pmax);
for i = 1:np
    for j = 1:length(PQpom{i})
        PQ(i,j) = PQpom{i}(j);
    end
end 

Qg{1} = nq;
Qg{2} = 0;
Qg{3} = QP;
Qg{4} = QQ;
Qg{5} = QE;
Qg{6} = QSx;
Qg{7} = QSy;
Qg{8} = QTx;
Qg{9} = QTy;
Qg{10} = QO;
Qg{11} = Qpricka;
Qg{12} = QLR;

Eg{1} = ne;
Eg{2} = EP;
Eg{3} = EQ;
Eg{4} = Enx;
Eg{5} = Eny;
Eg{6} = Etyp;

Pg{1} = np;
Pg{2} = X;
Pg{3} = Y;
Pg{4} = PQpoc;
Pg{5} = PQ;
Pg{6} = P(:,3); % typ bodu

geometrie{1} = Pg;
geometrie{2} = Qg;
geometrie{3} = Eg;

save 'implicitne\geometrie' geometrie;

waitbar(1,h);
close(h);


function vypoctiGeometrii_expl(P,Q,E)
h = waitbar(0,'Please wait...');
PX = P(:,1);
PY = P(:,2);
Etyp = E(:,3);

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
    Qr(Q(i,3),Q(i,1)) = i;
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
    TT(i,1) = Qr(Q(i,2),Q(i,1))-1; % indexy sousednich ctvercu
    TT(i,2) = Qr(Q(i,3),Q(i,2))-1;
    TT(i,3) = Qr(Q(i,1),Q(i,3))-1;
    TE(i,1) = Er(Q(i,1),Q(i,2)); % index prvni strany
    TE(i,2) = Er(Q(i,2),Q(i,3));
    TE(i,3) = Er(Q(i,3),Q(i,1));
end

waitbar(0.75,h);

for i = 1:nq
    for j = 1:3
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
geometrie{5} = 3*ones(nq,1);

save 'explicitne\geometrie' geometrie;

waitbar(1,h);
close(h);


% _________________________________________________________________________
function s = obsah(A,B,C)
v1 = B-A;
v2 = C-A;
s = abs(v1(1)*v2(2) - v1(2)*v2(1))/2;


% _________________________________________________________________________
function [T,P,Ptyp] = adaptuj(T,P,Ptyp)

np = length(P(:,1)); % pocet bodu
nt = length(T(:,1)); % pocet trojuhelniku

% vytvoreni vektoru sten
E = [T(:,[1,2]); T(:,[2,3]); T(:,[3,1])];
E = unique(sort(E,2),'rows');
ne = length(E(:,1));

% vytvoreni matice stran
Er = sparse(np,np);
for i = 1:ne
    Er(E(i,1),E(i,2)) = i;
    Er(E(i,2),E(i,1)) = i;
end

QE = zeros(nt,3);
% plneni matice QE
for i = 1:nt
    QE(i,1) = Er(T(i,1),T(i,2)); % index prvni strany
    QE(i,2) = Er(T(i,2),T(i,3));
    QE(i,3) = Er(T(i,3),T(i,1));
end

BolAd = zeros(nt,1);
% oznaceni trojuhelniku u sten_____________________________________________
for i = 1:nt
    for k = 1:3
        if(Ptyp(T(i,k)) == 2)
            BolAd(i) = 1;
        end
    end
end

% oznaceni vedlejsich trojuhelniku
Ppom = zeros(np,1);
I = 1:nt;
Ih = I(BolAd(I) == 1);
Ppom(T(Ih,1)) = 1; % oznaceni bodu oznacenych trojuhelniku
Ppom(T(Ih,2)) = 1;
Ppom(T(Ih,3)) = 1;

for i = 1:nt
    for k = 1:3
        if(Ppom(T(i,k)) == 1)
            BolAd(i) = 1;
        end
    end
end
% _________________________________________________________________________

% vlastni adaptace
Ead = zeros(ne,1);

for i = 1:nt % oznaceni trojuhelniku, ktere se budou adaptovat
    if(BolAd(i) == 1)
        for k = 1:3
            Ead(QE(i,k)) = 1;
        end
    end
end

while(1 == 1)
    bol = 1;
    for i = 1:nt % kontrola spravnosti deleni trojuhelniku, spravne jsou 1 nebo 3 delici strany
        s = 0;
        for k = 1:3 % soucet delicich stran
            s = s + Ead(QE(i,k));
        end
        if(s == 2) % spatny pocet delicich stran
            for k = 1:3
                Ead(QE(i,k)) = 1;
            end
            bol = 0;
        end
    end
    
    if(bol == 1)
        break;
    end
end

% vytvoreni novych bodu
ES = zeros(ne,1);
for i = 1:ne
    if(Ead(i) == 1)
        Pn = (P(E(i,1),:) + P(E(i,2),:))/2;
        P = [P; Pn];
        if(Ptyp(E(i,1)) == 2 && Ptyp(E(i,2)) == 2)
            Ptyp = [Ptyp; 2];
        elseif(Ptyp(E(i,1)) == 1 && Ptyp(E(i,2)) == 1)
            Ptyp = [Ptyp; 1];
        elseif(Ptyp(E(i,1)) == -1 && Ptyp(E(i,2)) == -1)
            Ptyp = [Ptyp; -1];
        elseif(Ptyp(E(i,1)) == -2 && Ptyp(E(i,2)) == -2)
            Ptyp = [Ptyp; -2];
        elseif((Ptyp(E(i,1)) == 1 && Ptyp(E(i,2)) == -1) || (Ptyp(E(i,1)) == -1 && Ptyp(E(i,2)) == 1) || (Ptyp(E(i,1)) == 1 && Ptyp(E(i,2)) == -2) || (Ptyp(E(i,1)) == -2 && Ptyp(E(i,2)) == 1))
            Ptyp = [Ptyp; 1];
        elseif((Ptyp(E(i,1)) == 2 && Ptyp(E(i,2)) == -1) || (Ptyp(E(i,1)) == -1 && Ptyp(E(i,2)) == 2) || (Ptyp(E(i,1)) == 2 && Ptyp(E(i,2)) == -2) || (Ptyp(E(i,1)) == -2 && Ptyp(E(i,2)) == 2))
            Ptyp = [Ptyp; 2];
        else
            Ptyp = [Ptyp; 0];
        end
        np = np+1;
        ES(i) = np;
    end
end

Tpom = T;          
% vlastni deleni
for i = 1:nt
   v = [Ead(QE(i,1)),Ead(QE(i,2)),Ead(QE(i,3))];
   if(sum(v) ~= 0)
       switch sum(v)
           case 1
               for k = 1:3
                   if(v(k) == 1)
                       p = k;
                   end
               end
               T(i,:) = [ES(QE(i,p)),Tpom(i,circ(p+2)),Tpom(i,p)];
               T = [T; [ES(QE(i,p)),Tpom(i,circ(p+1)),Tpom(i,circ(p+2))]];
           case 3
               T(i,:) = [ES(QE(i,1)),ES(QE(i,2)),ES(QE(i,3))];
               T = [T; [Tpom(i,1),ES(QE(i,1)),ES(QE(i,3))]];
               T = [T; [Tpom(i,2),ES(QE(i,2)),ES(QE(i,1))]];
               T = [T; [Tpom(i,3),ES(QE(i,3)),ES(QE(i,2))]];
       end
   end
end

% odstraneni mozneho spatneho oznaceni vnitrnich uzlu
% urceni vnitrnich bodu
    nt = length(T(:,1));
    np = length(P(:,1));
    uhel = zeros(nt,3);
    I = 1:nt;
    uhel(I,1) = abs(acos(((P(T(I,2),1)-P(T(I,1),1)).*(P(T(I,3),1)-P(T(I,1),1)) + (P(T(I,2),2)-P(T(I,1),2)).*(P(T(I,3),2)-P(T(I,1),2)))./(((P(T(I,2),1)-P(T(I,1),1)).^2 + (P(T(I,2),2)-P(T(I,1),2)).^2).^(1/2).*((P(T(I,3),1)-P(T(I,1),1)).^2 + (P(T(I,3),2)-P(T(I,1),2)).^2).^(1/2))));
    uhel(I,2) = abs(acos(((P(T(I,3),1)-P(T(I,2),1)).*(P(T(I,1),1)-P(T(I,2),1)) + (P(T(I,3),2)-P(T(I,2),2)).*(P(T(I,1),2)-P(T(I,2),2)))./(((P(T(I,3),1)-P(T(I,2),1)).^2 + (P(T(I,3),2)-P(T(I,2),2)).^2).^(1/2).*((P(T(I,1),1)-P(T(I,2),1)).^2 + (P(T(I,1),2)-P(T(I,2),2)).^2).^(1/2))));
    uhel(I,3) = abs(acos(((P(T(I,1),1)-P(T(I,3),1)).*(P(T(I,2),1)-P(T(I,3),1)) + (P(T(I,1),2)-P(T(I,3),2)).*(P(T(I,2),2)-P(T(I,3),2)))./(((P(T(I,1),1)-P(T(I,3),1)).^2 + (P(T(I,1),2)-P(T(I,3),2)).^2).^(1/2).*((P(T(I,2),1)-P(T(I,3),1)).^2 + (P(T(I,2),2)-P(T(I,3),2)).^2).^(1/2))));

    % uhel prirazeny jednotlivym bodum
    up = zeros(np,1);
    for i = 1:nt
        up(T(i,1)) = up(T(i,1)) + uhel(i,1);
        up(T(i,2)) = up(T(i,2)) + uhel(i,2);
        up(T(i,3)) = up(T(i,3)) + uhel(i,3);
    end

    % IH jsou indexy bodu na hranici
    I = 1:length(P(:,1));
    eps = 0.1;
    IH = I(up(I) >= 2*pi-eps);
    Ptyp(IH) = 0;

% _________________________________________________________________________
function j = circ(i)
    j = mod(i-1,3)+1;

    
    
% _________________________________________________________________________    
function P = vyhlad_sit(T,P,Ptyp)
    nt = length(T(:,1));
    np = length(P(:,1));
    DP = zeros(np,2);
    poc = zeros(np,1);
    
    for i = 1:nt
        DP(T(i,2),:) = DP(T(i,2),:) + P(T(i,1),:);
        poc(T(i,2),:) = poc(T(i,2),:) + 1;
        DP(T(i,3),:) = DP(T(i,3),:) + P(T(i,2),:);
        poc(T(i,3),:) = poc(T(i,3),:) + 1;
        DP(T(i,1),:) = DP(T(i,1),:) + P(T(i,3),:);
        poc(T(i,1),:) = poc(T(i,1),:) + 1;
    end
    
    I = 1:np;
    I = I(Ptyp == 0);
    P(I,1) = DP(I,1)./poc(I);
    P(I,2) = DP(I,2)./poc(I);
    