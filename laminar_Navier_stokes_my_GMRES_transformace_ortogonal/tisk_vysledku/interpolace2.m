function interpolace2(a,b)
% tato funkce interpoluje data ze site s oznacenim a do site s oznacenim b

g = load(['geometrie',num2str(a)]);
geometrie = g.geometrie;
PX1 = geometrie{1};
PY1 = geometrie{2};
TP1 = geometrie{3}+1;

g = load(['geometrie',num2str(b)]);
geometrie = g.geometrie;
PX2 = geometrie{1};
PY2 = geometrie{2};
TP2 = geometrie{3}+1;

W1 = load('W.txt');

% prevod hodnot do uzlu
nt = length(TP1);
I = 1:nt;
TTx1 = (PX1(TP1(I,1)) + PX1(TP1(I,2)) + PX1(TP1(I,3)))/3;
TTy1 = (PY1(TP1(I,1)) + PY1(TP1(I,2)) + PY1(TP1(I,3)))/3;

% pridani vnejsi ho ohraniceni
h = 1e4;
x_min = min(PX1)-h;
x_max = max(PX1)+h;
y_min = min(PY1)-h;
y_max = max(PY1)+h;
PXo = [x_min; x_max; x_max; x_min];
PYo = [y_min; y_min; y_max; y_max];
TTx1 = [TTx1; PXo];
TTy1 = [TTy1; PYo];
Wo = [ones(4,1), zeros(4,2), ones(4,1)];
W1 = [W1;Wo];

% vlastni interpolace do stredu bunek nove site
nt2 = length(TP2(:,1));
I = 1:nt2;
W = zeros(nt2,4);
TTx2 = (PX2(TP2(I,1)) + PX2(TP2(I,2)) + PX2(TP2(I,3)))/3;
TTy2 = (PY2(TP2(I,1)) + PY2(TP2(I,2)) + PY2(TP2(I,3)))/3;
for k = 1:4
    F = TriScatteredInterp(TTx1,TTy1,W1(:,k),'linear');
    W(I,k) = F(TTx2,TTy2);
end

% ulozeni dat
fid = fopen('W.txt','w');
for i = 1:nt2
    if(isnan(W(i,1)))
        W(i,:) = [3.1, 0, 0, 3.1];
    end
    fprintf(fid,'%15.12f %15.12f %15.12f %15.12f\n',W(i,:));
end
fclose(fid);






