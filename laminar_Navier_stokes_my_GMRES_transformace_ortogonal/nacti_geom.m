function nacti_geom

load geometrie;
PX = geometrie{1};
PY = geometrie{2};
TP = geometrie{3};
TT = geometrie{4};
typ = geometrie{5};
np = length(PX);
nt = length(TT(:,1));

% zapisuje matici TT do TT.txt
fid = fopen('TT.txt','w');
for i = 1:nt
    for j = 1:typ(i)
        fprintf(fid,'%i ',TT(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

% zapisuje matici TP do TP.txt
fid = fopen('TP.txt','w');
for i = 1:nt
    for j = 1:typ(i)
        fprintf(fid,'%i ',TP(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

% zapisuje vektor typ do typ.txt
fid = fopen('typ.txt','w');
for i = 1:nt
    fprintf(fid,'%i\n',typ(i));
end
fclose(fid);

% zapisuje body souboru PXY.txt
fid = fopen('PXY.txt','w');
for i = 1:np
    fprintf(fid,'%15.12f %15.12f\n',PX(i),PY(i));
end
fclose(fid);
















