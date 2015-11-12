function solver2(rad,iter,pokrac,CFL)

[L_ref, u_ref, ro_ref, p_ref, eta_ref, k_ref, Re, Pr, kapa] = referencni_hodnoty;

Re = 1000;
if(rad == 1)
    Re = -1;
end

% okrajove podminky
pi = 1;     % vstupni tlak
ri = 1;     % vstupni teplota T0!!!!!!!!
alfa = 0;   % vstupni uhel nabehu
% M = 0.95;
% po = pi/((1+(kapa-1)/2*M^2)^(kapa/(kapa-1)));   % vystupni tlak
po  = 0.737;

vstup_typ = 0; % vstupni okrajova podminka
ri2 = 1.4;
ui2 = 3;
vi2 = 0;
pi2 = 1;
ei2 = pi2/(kapa-1)+0.5*ri2*(ui2^2 + vi2^2);
n_vl = 6; % pocet vlaken
int_iter = 1; %rad; % pocet vnitrnich iteraci
tol = 1e-3;
eps_A = 1e-4; % tlumeni razu
eps_0 = 0; % tlumeni vsude
c_IP = 0;
if(eps_0 > 0 || Re > -1)
    c_IP = 1;
end

% load geometrie_dual;
load geometrie;
PX = geometrie{1};
TP = geometrie{3};
typ = geometrie{5};

np = length(PX);
nt = length(TP(:,1));

% zapisuje okrajove podminky a ridici hodnoty do souboru hodnoty.txt
fid = fopen('parameters_solver.txt', 'w');
fprintf(fid,'%6.8f\n', np);
fprintf(fid,'%6.8f\n', nt);
fprintf(fid,'%6.8f\n', CFL);
fprintf(fid,'%6.8f\n', rad);
fprintf(fid,'%6.8f\n', iter);
fprintf(fid,'%6.8f\n', n_vl);
fprintf(fid,'%6.8f\n', int_iter);
fprintf(fid,'%6.8f\n', tol);
fprintf(fid,'%6.8f\n', c_IP);
fprintf(fid,'%6.8f\n', eps_A);
fprintf(fid,'%6.8f\n', eps_0);
fclose(fid);

fid = fopen('parameters_equations.txt', 'w');
fprintf(fid,'%6.8f\n', pi);
fprintf(fid,'%6.8f\n', ri);
fprintf(fid,'%6.8f\n', alfa);
fprintf(fid,'%6.8f\n', po);
fprintf(fid,'%6.8f\n', kapa);
fprintf(fid,'%6.8f\n', Re);
fprintf(fid,'%6.8f\n', Pr);

fprintf(fid,'%6.8f\n', vstup_typ);
fprintf(fid,'%6.8f\n', ri2);
fprintf(fid,'%6.8f\n', ri2*ui2);
fprintf(fid,'%6.8f\n', ri2*vi2);
fprintf(fid,'%6.8f\n', ei2);
fclose(fid);

% generuje matici s pocatecnimi podminkami
if(pokrac == 0)
    if(vstup_typ == 0)
%         po = pi;
        Minl = sqrt((2/(kapa-1))*(-1 + (pi/po)^((kapa-1)/kapa)));
        Rinl = ri*(1+((kapa-1)/2)*Minl*Minl)^(1/(1-kapa));
        Vinl = Minl*sqrt(((kapa*po)/Rinl));
        uinl = Vinl*cos(alfa);
        vinl = Vinl*sin(alfa);
        Einl = po/(kapa-1)+0.5*Rinl*Vinl*Vinl;
        W = [Rinl, Rinl*uinl, Rinl*vinl, Einl];
    else
        W = [ri2, ri2*ui2, ri2*vi2, ei2];
    end
    fid = fopen('W.txt','w');
    for i = 1:nt
        fprintf(fid,'%15.12f %15.12f %15.12f %15.12f\n',W);
    end
    fclose(fid);
end

display(['Pocet bunek site je nt = ',num2str(nt),', pocet stupnu volnosti je ',num2str(stupen(rad,typ)),'.']);
tic
!java  -Xmx8000m Resic_implicit
toc

function n = stupen(r,typ)
n = 0;
for i = 1:length(typ)
    if(typ == 3)
        n = n + r*(r+1)/2;
    else
        n = n + r^2;
    end
end













