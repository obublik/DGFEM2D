function quad2tri

load 'P';
load 'Q'

% tvorba triangulace
tri = [];
for i = 1:length(Q(:,1));
        tri = [tri;[Q(i,1),Q(i,2),Q(i,4)]];
        tri = [tri;[Q(i,2),Q(i,3),Q(i,4)]];
end

mesh{1} = P;
mesh{2} = tri;

save 'mesh' mesh;