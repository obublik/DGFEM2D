function zobraz_sit

load 'geometrie';
PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3}+1;
typ = geometrie{5};
nt = length(TP(:,1));
typ
% tvorba triangulace
tri = [];
for i = 1:nt
    if(typ(i) == 3)
        tri = [tri;TP(i,1:3)];
    else
        tri = [tri;[TP(i,1),TP(i,2),TP(i,4)]];
        tri = [tri;[TP(i,2),TP(i,3),TP(i,4)]];
    end
end
length(tri(:,1))
figure('color','w');
hold on;
triplot(tri,PX,PY,'k');
axis equal;
box on;

% load 'rad_elementu.txt';
% for i = 1:nt
%     if(typ(i) == 3)
%         xc = sum(PX(TP(i,1:3)))/3;
%         yc = sum(PY(TP(i,1:3)))/3;
%     else
%         xc = sum(PX(TP(i,1:4)))/3;
%         yc = sum(PY(TP(i,1:4)))/3;
%     end
%     switch rad_elementu(i)
%         case 1
%             text(xc,yc,'1')
%         case 2
%             text(xc,yc,'2')
%         case 3
%             text(xc,yc,'3')
%     end
% end