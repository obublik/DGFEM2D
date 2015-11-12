function zobraz_sit

fid = fopen('We.txt','r');
radka = fgetl(fid);

figure
hold on
s = 0;
while(radka ~= -1)
    We = str2num(radka);
    nk = We(1);
    x = We(2:(nk+1));
    y = We((nk+2):(2*nk+1));
    radka = fgetl(fid);
    plot([x,x(1)],[y,y(1)])
    s = s + 1;
end
axis equal
display(['nt = ',num2str(s)]);