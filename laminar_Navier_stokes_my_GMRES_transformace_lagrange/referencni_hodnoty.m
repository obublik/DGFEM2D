function [L_ref, u_ref, ro_ref, p_ref, eta_ref, k_ref, Re, Pr, kapa] = referencni_hodnoty

R = 8314;
cp = 1005;
cv = cp-287;
kapa = cp/cv;
eta = 1.57e-5;
k = 0.0347;
Pr = cp*eta/k;

L_ref = 0.001;
p_ref = 978064;
T_ref = 292.76;
u_ref = sqrt(cv*T_ref);
ro_ref = p_ref/u_ref^2;
eta_ref = eta;
k_ref = 4;
Re = ro_ref*u_ref*L_ref/eta_ref;