function data = NACA(n)

designation = '0012';

t = str2num(designation(3:4))/100;
m = str2num(designation(1))/100;
p = str2num(designation(2))/10;

a0= 0.2969;
a1=-0.1260;
a2=-0.3516;
a3= 0.2843;
a4=-0.1015;

c = 1;
x = linspace(0,c,n+1)';

% neekvidistantni deleni
x = c*(x/c).^2;

yt = (t/0.2)*(a0*sqrt(x)+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);

xc1 = x(find(x<=p));
xc2 = x(find(x>p));
xc=[xc1 ; xc2];

if(p == 0)
    xu = x;
    yu = yt;

    xl = x;
    yl = -yt;
else
    yc1 = (m/p^2)*(2*p*xc1-xc1.^2);
    yc2 = (m/(1-p)^2)*((1-2*p)+2*p*xc2-xc2.^2);
    yc = [yc1 ; yc2];

    dyc1_dx = (m/p^2)*(2*p-2*xc1);
    dyc2_dx = (m/(1-p)^2)*(2*p-2*xc2);
    dyc_dx = [dyc1_dx ; dyc2_dx];
    theta = atan(dyc_dx);

    xu = x - yt.*sin(theta);
    yu = yc + yt.*cos(theta);

    xl = x + yt.*sin(theta);
    yl = yc - yt.*cos(theta);
end

X = [flipud(xl) ; xu(2:end)];
Y = [flipud(yl) ; yu(2:end)];

% aby byly koncove body na ose x
Y(1) = 0;
Y(length(Y)) = 0;

alfa = -5*pi/180;
xs = 1;
data = zeros(length(X),2);
data(:,1) = (X-xs)*cos(alfa) - Y*sin(alfa) + xs;
data(:,2) = (X-xs)*sin(alfa) + Y*cos(alfa);







