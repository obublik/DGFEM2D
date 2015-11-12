import java.io.*;
import java.util.StringTokenizer;

class Element{
// implementace DGFEM s ortogonalnimi polynomi
	// geometrie
	double[] PX;
	double[] PY;
	double[] nx;
	double[] ny;
	double[] S;
    double O;
	double tez;
    double xs, ys;

	// sit
	Element[] Sit;
    Stena[] Se;
	int[] TT;
    int[] TP;
	int nk; //pocet sten
    int ns; // globalni index (zacatek)
    int index;

	// konstanty
    double h = 1e-8;
    double c_IP_max;
    double c_IP;
    double eps_max;
    double eps_0;
    double eps;
    double tol;
    double[] Wpoc;
    int dif;
    int zdroj;
	
    // tolerance pro hodnotu hustoty
	double R_eps = 1e-2;
    
	// rovnice
	int nr; // pocet rovnic
    Rovnice Eq;
	
	// bazove funkce, koeficienty bazovych funkci
	int rad;      // rad presnosti elementu
	int ne; 	  // pocet bazovych funkci
    int[] baze_xi; // mocniny baze
    int[] baze_eta; // mocniny baze
    double[][] A_baze; //koeficieny bazovych polynomu
	int nip;	  // pocet integracnich bodu na stene
	double[] xip; // poloha integracnich bodu na stene
	double[] w;   // vahy Legendreova-Gaussova integracniho vzorce
    double[] Jac; // Jacobian transformace

	int ni_int;       // pocet integracnich bodu uvnitr elementu
    double[][] bp; // baricentricke souradnice pomocne
	double[][] b_int; // skutecne souradnice integracnich bodu
	double[] w_int;   // vahy
	double[][] Vb; // matice prechodu mezi dvema systemy bazi
	double[][] dbdx; // hodnota derivace bazove funkce ve vnitrnim integracnim bode
	double[][] dbdy;
    double k_el; // kvalita elementu
    double[] Is; // integral bazove funkce

	double[][] M;	 // matice hmotnosti
	double[][] Minv; // inverzni matice hmotnosti
    double[][] M_in;
    double[][] M_in_inv; // inverze M_in
    double[] RHS_loc;
    Soused[] M_sous; //matice sousedu
    Transformace Tr;
    
    // matice globalnich indexu a globalni matice s pravou stranou
    int [] gi_U;

	double[][] W;	  // hodnota W v n+1 te casove hladine
    double[][] Wo;   // hodnota v n te casove hladine
    double[][] Wo2; // hodnota v n-1 casove hladine
	double[][] K;
	double[][] G;
	double[] Wp;   // pracovni hodnota W
    double[][] E = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}; // jednotkova matice

	// hodnoty potrebne pro limiter
	double pom; // pomocna hodnota
	
	// casove hodnoty
	double dt;
    double dto;
	double t;
	
	// pomocne hodnoty
	int[] k_plus;


	Element(Element[] Sit, Rovnice Eq, double[] PX, double[] PY, int[] TP, int[] TT, int nk, int rad, double[] Wpoc, double tol, double c_IP_max, double eps_max, double eps_0) throws IOException{
		this.Sit = Sit;
		this.TT = TT;
        this.TP = TP;
		this.PX = PX;
		this.PY = PY;
		this.nk = nk;
		this.rad = rad;
        this.Wpoc = Wpoc;
        this.tol = tol;
        this.c_IP_max = c_IP_max;
        this.eps_max = eps_max;
        this.eps_0 = eps_0;
        
        // alokace objektu rovnic
        this.Eq = Eq;
        nr = Eq.nr;
        dif = Eq.difuze;
        zdroj = Eq.zdroj;

        k_plus = new int[nk];
		for(int k = 0; k < nk; k++){
			k_plus[k] = k+1;
		}
		k_plus[nk-1] = 0;
        
        // vypocet transformace mezi referencnim a realnym elementem
        Tr = new Transformace(nk,PX,PY);
        
        // vypocet geometrie
	    spocti_geom();
        
        // definuje bazove funkce
        bazove_funkce(rad);
        
        // definice vah integracnich vzorcu a vypocet matice hmotnosti
        generuj_rad_elementu();
        
        // generovani matice hmotnosti
	    generuj_Matice();
        
		// alokace poli
	    W = new double[ne][nr];
        Wo = new double[ne][nr];
        Wo2 = new double[ne][nr];
		for(int j = 0; j < nr; j++){
            for(int i = 0; i < ne; i++){
                W[i][j] = Wpoc[j];
                Wo[i][j] = Wpoc[j];
                Wo2[i][j] = Wpoc[j];
            }
		}
        M_in = new double[nr*ne][nr*ne];
        RHS_loc = new double[nr*ne];
        gi_U = new int[nr*ne];
	}
    
    void inicializuj(){
        // tvoreni objektu sten
        Se = new Stena[nk];
        int k_TT = -1;
        int kp_TT = -1;
        for(int k = 0; k < nk; k++){
            int kp = k_plus[k];
            for(int j = 0; j < nk; j++){ // hledani odpovidajicich si bodu na hranici
                if(TT[k] > -1){
                    if(TP[k] == Sit[TT[k]].TP[j])
                        k_TT = j;
                    if(TP[kp] == Sit[TT[k]].TP[j])
                        kp_TT = j;
                }
            }         
            Se[k] = new Stena(TT[k],k,kp,k_TT,kp_TT);
        }
    }
    
	// Aplikace limiteru
	void limiter(){
        if(rad > 1){
            double g_shock = shock_senzor();
            pom = g_shock;
            c_IP = c_IP_max;
            if(g_shock > tol/2 && tol > 1e-10){
                eps = eps_max*(g_shock - tol/2)*(2/tol);
            }
            else{
                eps = 0;
            }
        }
        else{
            eps = 0;
            c_IP = 0;
        }
	}
    
	// Generovani radku globalni matice a vektoru prave strany
	void sestav_lokalni_matice(){
        
		vynuluj_matice();
        
        // naplneni matic
//         double a1 = (2*dt + dto)/(dt*(dt+dto)); //3/(2*dt); 
//         double a2 = -(dt+dto)/(dt*dto); //-2/dt;
//         double a3 = dt/(dto*(dt+dto)); //1/(2*dt);
        
        double a1 = 1/dt; 
        double a2 = -1/dt;
        double a3 = 0;
        
        // vnitrni element - krivkovy i objemovy integral
        double[][] V = new double[ne][nr];
        
        double[][] Rw = R(V,M_sous);
        for(int m = 0; m < nr; m++){
            for(int i = 0; i < ne; i++){
                V[i][m] = h;
                double[][] Rwh = R(V,M_sous);
                RHS_loc[ne*m+i] = Rw[i][m];
                for(int j = 0; j < ne; j++){
                    RHS_loc[ne*m+i] = RHS_loc[ne*m+i] - M[i][j]*(a1*W[j][m] + a2*Wo[j][m] + a3*Wo2[j][m]);
                }
                for(int p = 0; p < nr; p++){
                    for(int j = 0; j < ne; j++){
                        M_in[ne*p+j][ne*m+i] = -(Rwh[j][p] - Rw[j][p])/h;
                    }
                }
                V[i][m] = 0;
            }
        }
        
        // pricteni matice hmotnosti
        for(int m = 0; m < nr; m++){
            for(int i = 0; i < ne; i++){
                for(int j = 0; j < ne; j++){
                     M_in[ne*m+i][ne*m+j] = M_in[ne*m+i][ne*m+j] + M[i][j]*a1;
                }
            }
        }
        
        // sousedni elementy - pouze krivkovy integral pres k-tou stenu
        for(int k = 0; k < nk; k++){
            Rw = R_ste(k,V,M_sous,null);
            for(int m = 0; m < nr; m++){
                for(int i = 0; i < M_sous[k].neR; i++){
                    M_sous[k].V[i][m] = h;
                    double[][] Rwh = R_ste(k,V,M_sous,null);
                    for(int p = 0; p < nr; p++){
                        for(int j = 0; j < ne; j++){
                            M_sous[k].MR[ne*p+j][M_sous[k].neR*m+i] = -(Rwh[j][p] - Rw[j][p])/h;
                        }
                    }
                    M_sous[k].V[i][m] = 0;
                }
             }
        }
	}
    
     // tato funkce vypocitava reziduum__________________________________________
	double[][] R(double[][] V, Soused[] Sous){
		double[][] K = new double[ne][nr];

		// vypocet toku hranici
		for(int k = 0; k < nk; k++){ // opakovani pres jednotlive steny
            K = R_ste(k, V, Sous, K);
		}
		
		if(rad > 1){
			// pripocteni objemovych integralu
            double[] W_int = new double[nr];
            double[] dW_int = new double[2*nr];
			double[] f = null;  // nevazky tok
			double[] g = null;
			double[] fv = null; // vazky tok
			double[] gv = null;
            double[] product = null;
            
			for(int p = 0; p < ni_int; p++){ // hodnoty funkci f a g v integracnich bodech
                for(int j = 0; j < nr; j++){
                    W_int[j] = 0;
                    dW_int[j] = 0;
                    dW_int[nr+j] = 0;
                    for(int m = 0; m < ne; m++){
                        W_int[j] = W_int[j] + (W[m][j] + V[m][j])*Vb[p][m];
                        dW_int[j] = dW_int[j] + (W[m][j] + V[m][j])*dbdx[p][m];
                        dW_int[nr+j] = dW_int[nr+j] + (W[m][j] + V[m][j])*dbdy[p][m];
                    }
                }
                
				f = Eq.fp(W_int, 1.0, 0.0);
				g = Eq.fp(W_int, 0.0, 1.0);
	
                if(dif == 1){
					fv = Eq.fpv(W_int,dW_int, 1.0, 0.0);
					gv = Eq.fpv(W_int,dW_int, 0.0, 1.0);
				}
                
                // produkce
                if(zdroj == 1){
                    product = Eq.P(W_int,dW_int);
                }
                
				for(int j = 0; j < nr; j++){
					for(int m = 0; m < ne; m++){
						K[m][j]  = K[m][j] + Jac[p]*w_int[p]*(f[j]*dbdx[p][m] + g[j]*dbdy[p][m]) - (eps+eps_0)*Jac[p]*w_int[p]*(dW_int[j]*dbdx[p][m] + dW_int[nr+j]*dbdy[p][m]);
                        if(dif == 1)
							K[m][j] = K[m][j] - Jac[p]*w_int[p]*(fv[j]*dbdx[p][m] + gv[j]*dbdy[p][m]);
                        if(zdroj == 1)
							K[m][j] = K[m][j] + Jac[p]*w_int[p]*product[j]*Vb[p][m];
					}
				}
			}
		}
        return K;
    }
    
    // tato funkce vypocitava reziduum__________________________________________
	double[][] R_ste(int k, double[][] V, Soused[] Sous, double[][] K){
        if(K == null){
            K = new double[ne][nr];
        }
		double[] fn = null;
        double[] fvn = null;
        double[] WL = new double[nr];
        double[] dWL =  new double[2*nr];
		double[] WR = new double[nr];
        double[] dWR =  new double[2*nr];
        
		// vypocet toku hranici
        for(int p = 0; p < nip; p++){ // opakovani pres integracni body
            // hodnoty v integracnim bode
            for(int j = 0; j < nr; j++){
                WL[j] = 0;
                dWL[j] = 0;
                dWL[nr+j] = 0;
                for(int m = 0; m < ne; m++){
                    WL[j] = WL[j] + (W[m][j] + V[m][j])*Se[k].VpL[p][m];
                    dWL[j] = dWL[j] + (W[m][j] + V[m][j])*Se[k].dxVpL[p][m];
                    dWL[nr+j] = dWL[nr+j] + (W[m][j] + V[m][j])*Se[k].dyVpL[p][m];
                }
            }
            if(TT[k] > -1){
                double[][] WRp = Sit[TT[k]].W;
                double[][] Vs = Sous[k].V;
                for(int j = 0; j < nr; j++){
                    WR[j] = 0;
                    dWR[j] = 0;
                    dWR[nr+j] = 0;
                    for(int m = 0; m < Sit[TT[k]].ne; m++){
                        WR[j] = WR[j] + (WRp[m][j] + Vs[m][j])*Se[k].VpR[p][m];
                        dWR[j] = dWR[j] + (WRp[m][j] + Vs[m][j])*Se[k].dxVpR[p][m];
                        dWR[nr+j] = dWR[nr+j] + (WRp[m][j] + Vs[m][j])*Se[k].dyVpR[p][m];
                    }
                }
            }
            else{
                WR = Eq.hodnota_WR(WL,nx[k],ny[k],TT[k]);
                for(int j = 0; j < 2*nr; j++)
                    dWR[j] = dWL[j];
            }
            
            // nevazky tok v integracnim bode
            fn = Eq.Fn(WL, WR, nx[k], ny[k], TT[k]);	// nevazky tok
            
            // vazky tok v integracnim bode
            if(dif == 1){
                fvn = Eq. Fvn(WL, WR, dWL, dWR, nx[k], ny[k], TT[k]); // vazky tok
            }

            for(int m = 0; m < ne; m++){
                double b = Se[k].VpL[p][m];
                double dbx = Se[k].dxVpL[p][m];
                double dby = Se[k].dyVpL[p][m];
                for(int j = 0; j < nr; j++){
                    K[m][j]  = K[m][j] - (w[p]*fn[j]*S[k]/2)*b  - c_IP*w[p]*(WL[j]-WR[j])*S[k]/2*b;
                    if(TT[k] > -1)
                        K[m][j]  = K[m][j]  + eps_0*(w[p]*((dWL[j]+dWR[j])/2*nx[k] + (dWL[nr+j]+dWR[nr+j])/2*ny[k])*S[k]/2)*b;
                    if(dif == 1){
                        K[m][j] = K[m][j] + (w[p]*fvn[j]*S[k]/2)*b;
                    }
                }
            }
        }
        
        return K;
    }
    
	// tato funkce vraci hodnotu m-te bazove funkce v bode x, y
	double baze(int m, double xi, double eta){
        double b = 0;
        for(int i = 0; i < ne; i++){
            b = b + A_baze[i][m]*Math.pow(xi,baze_xi[i])*Math.pow(eta,baze_eta[i]);
        }
        return b;
	}

	// tato funkce vraci hodnotu x-ove derivace m-te bazove funkce v bode x, y
	double D_baze(int m, double xi, double eta, char s){
        double dxi = 0;
		double deta = 0;
        for(int i = 0; i < ne; i++){
            if(baze_xi[i] > 0)
                dxi = dxi + A_baze[i][m]*baze_xi[i]*Math.pow(xi,baze_xi[i]-1)*Math.pow(eta,baze_eta[i]);
            if(baze_eta[i] > 0)
                deta = deta + A_baze[i][m]*baze_eta[i]*Math.pow(xi,baze_xi[i])*Math.pow(eta,baze_eta[i]-1);
        }
		if(s == 'x')
			return Tr.itdx(xi,eta,dxi,deta);
		else
			return Tr.itdy(xi,eta,dxi,deta);
	}
	
	// limiter zapornych hodnot ________________________________________________
	void lim_zap_hodnot(){ // limituje zaporne hodnoty hustoty
	}
	
	//__________________________________________________________________________
	// reziduum
	double rezid(){
        double rez = 0;
        for(int m = 0; m < nr; m++){
        	for(int j = 0; j < ne; j++){
            	rez = rez + Math.abs(W[j][m] - Wo[j][m])/dt;
          	}
        }
        
        return rez;
	}

    // ulozeni W do Wo
	void uloz_Wo(){
         for(int j = 0; j < nr; j++){
        	for(int m = 0; m < ne; m++){
                Wo2[m][j] = Wo[m][j];
            	Wo[m][j] = W[m][j];
          	}
        }
        dto = dt;
	}
    
    // ulozeni Wo do W (pri potizich s resenim)
    void Wo_to_W(){
        for(int j = 0; j < nr; j++){
        	for(int m = 0; m < ne; m++){
            	W[m][j] = Wo[m][j];
          	}
        }
	}
    
	//__________________________________________________________________________
	double delta_t(double CFL){ //vypocet maximalniho casoveho kroku
		double[] Ws = spocti_Ws();
    	double p = Eq.tlak(Ws);
    	double a = Math.sqrt(Eq.kapa*p/Ws[0]);
    	double u = Math.abs(Ws[1]/Ws[0]);
    	double v = Math.abs(Ws[2]/Ws[0]);
    	double lamA = u + a;
    	double lamB = v + a;
    	dt = CFL/((lamA + lamB)/tez);
        
        return dt;
	}

    //__________________________________________________________________________
    double shock_senzor(){
        double g_shock = 0;
        double pod = 0;
        for(int p = 0; p < ni_int; p++){ // hodnoty funkci f a g v integracnich bodech
				double[] W_int = new double[nr];
                for(int j = 0; j < nr; j++){
                    for(int m = 0; m < ne; m++){
                        W_int[j] = W_int[j] + W[m][j]*Vb[p][m];
                    }
                }
                double[] Ws = spocti_Ws();
                g_shock = g_shock + Jac[p]*w_int[p]*(W_int[0]-Ws[0])*(W_int[0]-Ws[0]);
                pod = pod + Jac[p]*w_int[p]*W_int[0]*W_int[0];
			}
        
		return (g_shock/pod);
    }

    
    // trida stena ============================================
    class Soused{
        double[][] MR;
        double[][] V;
        int nr, ne, neR, typ;
        
        Soused(int typ, int ne, int neR, int nr){
            if(typ > -1){
                MR = new double[nr*ne][nr*neR];
            }
            V = new double[neR][nr];
            this.nr = nr;
            this.ne = ne;
            this.neR = neR;
            this.typ = typ;
        }
        
        void vynuluj(){
            if(typ > -1)
                for(int i = 0; i < ne*nr; i++){
                    for(int j = 0; j < neR*nr; j++){
                        MR[i][j] = 0;
                    }
                }    
         }
    }
	// trida stena ============================================
    
    
    void vynuluj_matice(){
         for(int i = 0; i < ne*nr; i++){
             RHS_loc[i] = 0;
            for(int j = 0; j < ne*nr; j++){
                M_in[i][j] = 0;
            }
         }
        for(int k = 0; k < nk; k++){
             M_sous[k].vynuluj();
        }
    }
    
    int nastav_globalni_index_U(int s){
        for(int i = 0; i < ne*nr; i++){
            gi_U[i] = s;
            s = s + 1;
        }
        return s;
    }
    
    void alokuj_sousedy(){
        M_sous = new Soused[nk];
        for(int k = 0; k < nk; k++){
            if(TT[k] > -1)
                M_sous[k] = new Soused(TT[k],ne,Sit[TT[k]].ne,nr);
            else
                M_sous[k] = new Soused(TT[k],ne,0,nr);
		}
    }
    
    void vypocti_inverzi_predpodminovace(){
        // vypocitava inverzi diagonaly, ktera se pouzije pro predpodminovac
        M_in_inv = inverze(M_in);
    }
    
    void GM_resid(double[] x, double[] r, int par){
        int n = nr*ne;
        double[] p = new double[n];
        for(int i = 0; i < n; i++){
            if(par == 1){
                p[i]  = RHS_loc[i];
                for(int j = 0; j < n; j++){
                    p[i] = p[i] - M_in[i][j]*x[gi_U[j]];
                }
            }
            else{
                p[i]  = 0;
                for(int j = 0; j < n; j++){
                    p[i] = p[i] + M_in[i][j]*x[gi_U[j]];
                }
            }
            for(int k = 0; k < nk; k++){
                if(TT[k] > -1){
                    for(int j = 0; j < Sit[TT[k]].nr*Sit[TT[k]].ne; j++){
                        p[i] = p[i] + M_sous[k].MR[i][j]*x[Sit[TT[k]].gi_U[j]];
                    }
                }
            }
        }
        for(int i = 0; i < n; i++){
            int ig = gi_U[i];
            r[ig] = 0;
            for(int j = 0; j < n; j++){
                r[ig] = r[ig] + M_in_inv[i][j]*p[j];
            }
        }
    }
    
    double sqr(){
        double nrm = 0;
         for(int i = 0; i < nr*ne; i++){
            nrm = nrm + RHS_loc[i]*RHS_loc[i];
        }
        return nrm;
    }
    
    void uloz_reseni(double[] x){
        for(int m = 0; m < nr; m++){
            for(int j = 0; j < ne; j++){
                int s = gi_U[ne*m+j];
                W[j][m] = W[j][m] + x[s];
            }
        }
    }
    
   // vypocet geometrie _______________________________________________________
	void spocti_geom(){
        double[] Sx = new double[nk];
        for(int k = 0; k < nk; k++){
            int kp = k_plus[k];
        	Sx[k] = (PY[kp] - PY[k]);
        }
        
        double[] Sy = new double[nk];
        for(int k = 0; k < nk; k++){
        	Sy[k] = -(PX[k_plus[k]] - PX[k]);
        }

        S = new double[nk];
        nx = new double[nk];
        ny = new double[nk];
        for(int k = 0; k < nk; k++){
        	S[k] = Math.sqrt(Sx[k]*Sx[k]+Sy[k]*Sy[k]);
        	nx[k] = Sx[k]/S[k];
        	ny[k] = Sy[k]/S[k];
        }
		
        xs = 0;
        ys = 0;
        for(int k = 0; k < nk; k++){
        	xs = xs + PX[k]/nk;
        	ys = ys + PY[k]/nk;
        }
        
        O = 0;
        for(int k = 0; k < nk; k++){
        	int kp = k_plus[k];
            O = O + Math.abs((PX[k]-xs)*(PY[kp]-ys)-(PX[kp]-xs)*(PY[k]-ys))/2;
        }
        
        tez = 0;
        for(int k = 0; k < nk; k++){
        	tez = tez + S[k];
        }
        tez = nk*O/tez;
        
        // kvalita elementu
        double pic = 3.141592653589793;
        double Sc = 0;
        for(int k = 0; k < nk; k++){
        	Sc = Sc + S[k];
        }
        k_el = Sc/(2*Math.sqrt(pic*O));
    }

    //________________________________
    void bazove_funkce(int rad){
        if(nk == 3){ // trojuhelnik
            ne = 1; // pocet bazovych funkci
            if(rad > 1)
                for(int i = 2; i <= rad; i++){
                    ne = ne + i;
                }
            baze_xi = new int[ne];
            baze_eta = new int[ne];
            int s = 0;
            for(int i = 0; i < rad; i++){
                for(int j = 0; j <= i; j++){
                    baze_xi[s] = i-j;
                    baze_eta[s] = j;
                    s = s + 1;
                }
            }

            if(rad == 1){
                A_baze = new double[1][1];
                A_baze[0][0] = 1;
            }
            else{
                double[][] b = barycentr(rad,ne);
                double[][] V = new double[ne][ne];
                V = new double[ne][ne];
                for(int i = 0; i < ne; i++){
                    double xi = b[i][1];
                    double eta = b[i][2];
                    for(int j = 0; j < ne; j++){
                        V[i][j] = Math.pow(xi,baze_xi[j])*Math.pow(eta,baze_eta[j]);
                    }
                }
                A_baze = inverze(V);
            }
        }
        else{ // ctverec
            ne = rad*rad; // pocet bazovych funkci
            
            baze_xi = new int[ne];
            baze_eta = new int[ne];
            int s = 0;
            for(int i = 0; i < rad; i++){
                for(int j = 0; j <= i; j++){
                    baze_xi[s] = i-j;
                    baze_eta[s] = j;
                    s = s + 1;
                }
            }
            for(int i = rad-1; i > 0 ; i--){
                for(int j = 0; j < i; j++){
                    baze_xi[s] = rad-1 - j;
                    baze_eta[s] = rad - i + j;
                    s = s + 1;
                }
            }
            if(rad == 1){
                A_baze = new double[1][1];
                A_baze[0][0] = 1;
            }
            else{
                double[] b = linspace(1, rad);
                double[][] V = new double[ne][ne];
                V = new double[ne][ne];
                s = 0;
                for(int i = 0; i < rad; i++){
                    for(int j = 0; j < rad; j++){
                        double xi = b[i];
                        double eta = b[j];
                        for(int k = 0; k < ne; k++){
                            V[s][k] = Math.pow(xi,baze_xi[k])*Math.pow(eta,baze_eta[k]);
                        }
                        s = s + 1;
                    }
                }
                A_baze = inverze(V);
            }
        }
    }
    
    //__________________________________________________________________________
	void generuj_rad_elementu() throws IOException{
		
        // nacteni 1D integracnich vzorcu
        nacti_gauss(rad);
        
        // nacteni 2D integracnich vzorcu
        int rad_mod;
        switch(nk){
            case 3:
                rad_mod = 1;
                if(rad > 1){
                     rad_mod = 2*(rad-1)+1;
                }
                nacti_gauss_tri(rad_mod);
            break;
            
            case 4:
                rad_mod = 1;
                if(rad > 1){
                     rad_mod = 4*(rad-1)+1;
                }
                nacti_gauss_quad(rad_mod); 
            break;
        }
        Jac = new double[ni_int];
        for(int i = 0; i < ni_int; i++){
            Jac[i] = Tr.Jacobian(b_int[i][0],b_int[i][1]);
        }
	}
    
	//__________________________________________________________________________
	void generuj_Matice(){ // funkce pro generovani matic
        // pomocne matice
        Vb = new double[ni_int][ne];
        dbdx = new double[ni_int][ne];
        dbdy = new double[ni_int][ne];
        for(int p = 0; p < ni_int; p++){
            for(int m = 0; m < ne; m++){
                double xi = b_int[p][0];
                double eta = b_int[p][1];
                Vb[p][m] = baze(m,xi,eta);
                dbdx[p][m] = D_baze(m,xi,eta,'x');
                dbdy[p][m] = D_baze(m,xi,eta,'y');
            }
        }
        
		// matice hmotnosti a tuhosti v x a y
		M = new double[ne][ne];
        Is = new double[ne];
		for(int i = 0; i < ne; i++){
			for(int j = 0; j < ne; j++){
				for(int p = 0; p < ni_int; p++){
                    double xi = b_int[p][0];
                    double eta = b_int[p][1];
					M[i][j]  = M[i][j] + Jac[p]*w_int[p]*baze(i,xi,eta)*baze(j,xi,eta);
				}
			}
            for(int p = 0; p < ni_int; p++){
                double x = b_int[p][0];
                double y = b_int[p][1];
                Is[i] = Is[i] + Jac[p]*w_int[p]*baze(i,x,y)/O;
            }
		}
		Minv = inverze(M);
	}
    
    // trida transformace ============================================
    class Transformace{
        double[] A, B;
        double pic = 3.141592653589793;
        double[] xi_ref;
        double[] eta_ref;
        double xi_s, eta_s;
        
        Transformace(int nk, double[] PX, double[] PY){
            double[][] V = null;
            switch(nk){
                case 3:
                    xi_ref = new double[3];
                    eta_ref = new double[3];
                    xi_ref[0] = 0; xi_ref[1] =  1; xi_ref[2] =  0;
                    eta_ref[0] = 0; eta_ref[1] = 0; eta_ref[2] = 1;
                    xi_s = 1.0/3;
                    eta_s = 1.0/3;
                    V = new double[3][3];
                    for(int i = 0; i < 3; i++){
                        V[i][0] = xi_ref[i];
                        V[i][1] = eta_ref[i];
                        V[i][2] = 1;
                    }
                break;
                
                case 4:
                    xi_ref = new double[4];
                    eta_ref = new double[4];
                    xi_ref[0] = 0; xi_ref[1] =  1; xi_ref[2] =  1; xi_ref[3] =  0;
                    eta_ref[0] = 0; eta_ref[1] = 0; eta_ref[2] = 1; eta_ref[3] = 1;
                    xi_s = 0.5;
                    eta_s = 0.5;
                    V = new double[4][4];
                    for(int i = 0; i < 4; i++){
                        V[i][0] = xi_ref[i]*eta_ref[i];
                        V[i][1] = xi_ref[i];
                        V[i][2] = eta_ref[i];
                        V[i][3] = 1;
                    }
                break;
            }
            A = nasob_MV(inverze(V),PX);        
            B = nasob_MV(inverze(V),PY);
        }
        
        double Jacobian(double xi, double eta){
            if(nk == 3)
                return Math.abs(A[0]*B[1]-A[1]*B[0])/2;
            else
                return Math.abs((A[0]*eta+A[1])*(B[0]*xi+B[2])-(A[0]*xi+A[2])*(B[0]*eta+B[1]));
        }
        
        double[][] Jacobian_inverze(double xi, double eta){
            double [][] J = new double[2][2];
            switch(nk){
                case 3:
                    J[0][0] = A[0];
                    J[0][1] = A[1];
                    J[1][0] = B[0];
                    J[1][1] = B[1];
                break;

                case 4:
                    J[0][0] = A[0]*eta+A[1];
                    J[0][1] = A[0]*xi+A[2];
                    J[1][0] = B[0]*eta+B[1];
                    J[1][1] = B[0]*xi+B[2];
                break;
            }
            return inverze(J);
        }
        
        double itdx(double xi, double eta, double dxi, double deta){
            double[][] iJ = Jacobian_inverze(xi,eta);
            return (iJ[0][0]*dxi + iJ[1][0]*deta);
        }
        
        double itdy(double xi, double eta, double dxi, double deta){
            double[][] iJ = Jacobian_inverze(xi,eta);
            return (iJ[0][1]*dxi + iJ[1][1]*deta);
        }
    }
    // trida transformace ============================================
    
    
    // trida stena =================================================
    class Stena{
        int TT;
        double[] xp, yp;
        double[][] VpL, VpR, dxVpL, dxVpR, dyVpL, dyVpR;
        
        Stena(int TT, int i1, int i2, int j1, int j2){
            this.TT = TT;
            VpL = new double[nip][ne];
            dxVpL = new double[nip][ne];
            dyVpL = new double[nip][ne];
            for(int p = 0; p < nip; p++){
                double xi = (1+xip[p])/2*Tr.xi_ref[i1]+(1-xip[p])/2*Tr.xi_ref[i2];
                double eta = (1+xip[p])/2*Tr.eta_ref[i1]+(1-xip[p])/2*Tr.eta_ref[i2];
                for(int m = 0; m < ne; m++){
                    VpL[p][m] = baze(m,xi,eta);
                    dxVpL[p][m] = D_baze(m,xi,eta,'x');
                    dyVpL[p][m] = D_baze(m,xi,eta,'y');
                }
            }
            if(TT > -1){
                VpR = new double[nip][Sit[TT].ne];
                dxVpR = new double[nip][Sit[TT].ne];
                dyVpR = new double[nip][Sit[TT].ne];
                for(int p = 0; p < nip; p++){
                    double xi = (1+xip[p])/2*Sit[TT].Tr.xi_ref[j1]+(1-xip[p])/2*Sit[TT].Tr.xi_ref[j2];
                    double eta = (1+xip[p])/2*Sit[TT].Tr.eta_ref[j1]+(1-xip[p])/2*Sit[TT].Tr.eta_ref[j2];
                    for(int m = 0; m < Sit[TT].ne; m++){
                        VpR[p][m] = Sit[TT].baze(m,xi,eta);
                        dxVpR[p][m] = Sit[TT].D_baze(m,xi,eta,'x');
                        dyVpR[p][m] = Sit[TT].D_baze(m,xi,eta,'y');
                    }
                }
            }
        }
    }
    
    double[][] barycentr(int rad, int ne){
        double[][] b;
        if(rad == 1){
            b = new double[1][3];
            b[0][0] = 1.0/3;
            b[0][1] = 1.0/3;
            b[0][2] = 1.0/3;
        }
        else{
            double[] x = linspace(1,rad);
            b = new double[ne][3];
            int s = 0;
            for(int i = 0; i < rad; i++){
                double[] y = linspace(x[i],i+1);
                for(int j = 0; j < i+1; j++){
                    b[s][0] = 1-x[i];
                    b[s][1] = y[j];
                    b[s][2] = x[i] - y[j];
                    s = s + 1;
                }
            }
        }
        
        return b;
    }
    
    double[] linspace(double a, int n){
        double[] x = new double[n];
        if(n == 1)
            x[0] = a/2;
        else{
            for(int i = 0; i < n; i++){
                x[i] = a*i/(n-1);
            }
        }
        return x;
    }
    
	//__________________________________________________________________________
    double sign(double a){ // signum
    	if(a > 1)
    		return 1;
    	else
    		return -1;
    }
	
	//__________________________________________________________________________
    double nasob_skalar(double[] A, double[] B){ // skalarni nasobeni
    	int n = A.length;
    	double C = 0;
    	for(int i = 0; i < n; i++){
    		C = C + A[i]*B[i];
    	}
    	return C;
    }

	//_____________________________________________________________________________
    double bil_form(double[][] A, double[] x){ // bilinearni forma y = x*A*x'
    	int ma = A.length;
    	double[] B = new double[ma];
    	for(int i = 0; i < ma; i++){
    		B[i] = nasob_skalar(A[i],x);
    	}
    	return nasob_skalar(x,B);
    }

	//____________________________________________________________________________
    double[] spocti_Ws(){ // funkce vraci stredni hodnotu W na kontrolnim objemu K
    	double[] Ws = new double[nr];
		for(int j = 0; j < nr; j++){
            for(int m = 0; m < ne; m++){
                Ws[j] = Ws[j] + W[m][j]*Is[m];
            }
		}
		return Ws;
    }

    //__________________________________________________________________________
    double[][] nasob(double[][] A, double[][] B){ // nasobeni matic
    	int ma = A.length;
    	int na = A[0].length;
    	int nb = B[0].length;
    	double[][] C = new double[ma][nb];
    	for(int i = 0; i < ma; i++){
    		for(int j = 0; j < nb; j++){
    			for(int k = 0; k < na; k++){
    				C[i][j] = C[i][j] + A[i][k]*B[k][j];
    			}
    		}
    	}
    	return C;
    }
	
    //__________________________________________________________________________
    double[][] nasob_Mc(double[][] A, double b){ // nasobeni matic
    	int ma = A.length;
    	int na = A[0].length;
    	double[][] C = new double[ma][na];
    	for(int i = 0; i < ma; i++){
    		for(int j = 0; j < na; j++){
                C[i][j] = C[i][j] + A[i][j]*b;
    		}
    	}
    	return C;
    }
    
	//__________________________________________________________________________
    double[] nasob_MV(double[][] A, double[] b){ // nasobeni matice s vektorem
    	int ma = A.length;
    	int na = A[0].length;
    	double[] c = new double[ma];
    	for(int i = 0; i < ma; i++){
    		for(int j = 0; j < na; j++){
    			c[i] = c[i] + A[i][j]*b[j];
    		}
    	}
    	return c;
    }
	
	//__________________________________________________________________________
    double[][] secti(double[][] A, double[][] B){ // scitani matic
    	int ma = A.length;
    	int na = A[0].length;
    	double[][] C = new double[ma][na];
    	for(int i = 0; i < ma; i++){
    		for(int j = 0; j < na; j++){
    			C[i][j] = A[i][j] + B[i][j];
    		}
    	}
    	return C;
    }
    
    //__________________________________________________________________________
    double[][] prumer(double[][] A, double[][] B){ // vypocet aritmetickeho prumeru
    	int ma = A.length;
    	int na = A[0].length;
    	double[][] C = new double[ma][na];
    	for(int i = 0; i < ma; i++){
    		for(int j = 0; j < na; j++){
    			C[i][j] = (A[i][j] + B[i][j])/2;
    		}
    	}
    	return C;
    }
	
	//__________________________________________________________________________
    double[][] odecti(double[][] A, double[][] B){ // scitani matic
    	int ma = A.length;
    	int na = A[0].length;
    	double[][] C = new double[ma][na];
    	for(int i = 0; i < ma; i++){
    		for(int j = 0; j < na; j++){
    			C[i][j] = A[i][j] - B[i][j];
    		}
    	}
    	return C;
    }
	
	//__________________________________________________________________________
    double[] odmocni(double[] A){ // nasobeni obsahem a odmocneni
    	int m = A.length;
    	double[] C = new double[m];
    	for(int i = 0; i < m; i++){
    		C[i] = Math.sqrt(A[i]);
    	}
    	return C;
    }
	
    //__________________________________________________________________________
    double[][] inverze(double ain[][]){ // inverze matice
	    int n = ain.length;
	    double a[][] = new double[n][n];
	    for(int i = 0; i < n; i++){
	    	for(int j = 0; j < n; j++){
	    		a[i][j] = ain[i][j];
	    	}
	    }

	    double x[][] = new double[n][n];
	    double b[][] = new double[n][n];
	    int index[] = new int[n];
	    for (int i=0; i<n; ++i) b[i][i] = 1;

	 // Transform the matrix into an upper triangle
	    gaussian(a, index);

	 // Update the matrix b[i][j] with the ratios stored
	    for (int i=0; i<n-1; ++i)
	      for (int j=i+1; j<n; ++j)
	        for (int k=0; k<n; ++k)
	          b[index[j]][k]
	            -= a[index[j]][i]*b[index[i]][k];

	 // Perform backward substitutions
	    for (int i=0; i<n; ++i) {
	      x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
	      for (int j=n-2; j>=0; --j) {
	        x[j][i] = b[index[j]][i];
	        for (int k=j+1; k<n; ++k) {
	          x[j][i] -= a[index[j]][k]*x[k][i];
	        }
	        x[j][i] /= a[index[j]][j];
	      }
	    }
	  return x;
	}

	// Method to carry out the partial-pivoting Gaussian
	// elimination.  Here index[] stores pivoting order.

	void gaussian(double a[][],int index[]){
	    int n = index.length;
	    double c[] = new double[n];

	 // Initialize the index
	    for (int i=0; i<n; ++i) index[i] = i;

	 // Find the rescaling factors, one from each row
	    for (int i=0; i<n; ++i) {
	      double c1 = 0;
	      for (int j=0; j<n; ++j) {
	        double c0 = Math.abs(a[i][j]);
	        if (c0 > c1) c1 = c0;
	      }
	      c[i] = c1;
	    }

	 // Search the pivoting element from each column
	    int k = 0;
	    for (int j=0; j<n-1; ++j) {
	      double pi1 = 0;
	      for (int i=j; i<n; ++i) {
	        double pi0 = Math.abs(a[index[i]][j]);
	        pi0 /= c[index[i]];
	        if (pi0 > pi1) {
	          pi1 = pi0;
	          k = i;
	        }
	      }

	   // Interchange rows according to the pivoting order
	      int itmp = index[j];
	      index[j] = index[k];
	      index[k] = itmp;
	      for (int i=j+1; i<n; ++i) {
	        double pj = a[index[i]][j]/a[index[j]][j];

	     // Record pivoting ratios below the diagonal
	        a[index[i]][j] = pj;

	     // Modify other elements accordingly
	        for (int l=j+1; l<n; ++l)
	          a[index[i]][l] -= pj*a[index[j]][l];
	      }
	    }
  	}

    void nacti_gauss(int r) throws IOException{
        String str = "gauss/" + r + ".txt";
        File soubor = new File(str);
        FileReader fr = new FileReader(str);
        BufferedReader in = new BufferedReader(fr);
        String c;
        String radka;
        
        radka = in.readLine();
        nip = Integer.valueOf(radka).intValue();
        xip = new double[nip];
        w = new double[nip];
        for(int i = 0; i < nip; i++){
            radka = in.readLine();
            xip[i] = Double.valueOf(radka).doubleValue();
        }
        for(int i = 0; i < nip; i++){
            radka = in.readLine();
            w[i] = Double.valueOf(radka).doubleValue();
        }
        fr.close();
    }
    
    void nacti_gauss_tri(int r) throws IOException{
        String str = "gauss_tri/" + r + ".txt";
        File soubor = new File(str);
        FileReader fr = new FileReader(str);
        BufferedReader in = new BufferedReader(fr);
        String c;
        String radka;
        
        radka = in.readLine();
        ni_int = Integer.valueOf(radka).intValue();
       
        bp = new double[ni_int][3];
        w_int = new double[ni_int];
        for(int i = 0; i < ni_int; i++){
            radka = in.readLine();
            StringTokenizer st = new StringTokenizer(radka);
            for(int j = 0; j < 3; j++){
                c = st.nextToken();
                bp[i][j] = Double.valueOf(c).doubleValue();
            }
        }
        for(int i = 0; i < ni_int; i++){
            radka = in.readLine();
            w_int[i] = Double.valueOf(radka).doubleValue();
        }
        
        b_int = new double[ni_int][2];
        for(int i = 0; i < ni_int; i++){
            b_int[i][0] = bp[i][0]*Tr.xi_ref[0] + bp[i][1]*Tr.xi_ref[1] + bp[i][2]*Tr.xi_ref[2];
            b_int[i][1] = bp[i][0]*Tr.eta_ref[0] + bp[i][1]*Tr.eta_ref[1] + bp[i][2]*Tr.eta_ref[2];
        }
        fr.close();
    }
    
    void nacti_gauss_quad(int r) throws IOException{
        String str = "gauss/" + r + ".txt";
        File soubor = new File(str);
        FileReader fr = new FileReader(str);
        BufferedReader in = new BufferedReader(fr);
        String c;
        String radka;
        
        radka = in.readLine();
        int n1 = Integer.valueOf(radka).intValue();
       
        double[] x1 = new double[n1];
        double[] w1 = new double[n1];
        for(int i = 0; i < n1; i++){
            radka = in.readLine();
            x1[i] = Double.valueOf(radka).doubleValue();
        }
        for(int i = 0; i < n1; i++){
            radka = in.readLine();
            w1[i] = Double.valueOf(radka).doubleValue();
        }
        fr.close();
        
        ni_int = n1*n1;
        b_int = new double[ni_int][2];
        w_int = new double[ni_int];
        int s = 0;
        for(int i = 0; i < n1; i++){
            for(int j = 0; j < n1; j++){
                b_int[s][0] = (1+x1[i])/2;
                b_int[s][1] = (1+x1[j])/2;
                w_int[s] = w1[i]*w1[j]/4;
                s = s + 1;
            }
        }
    }
    
  	//__________________________________________________________________________
  	void tiskV(double[] V){ // tiskne vektor
  		if(V == null)
  			System.out.println(V);
  		else{
	  		int m = V.length;
	  		for(int i = 0; i < m; i++){
	  			System.out.println(V[i]);
	  		}
	  		System.out.println();
  		}
  	}
    
    //__________________________________________________________________________
  	void tiskV_int(int[] V){ // tiskne vektor
  		if(V == null)
  			System.out.println(V);
  		else{
	  		int m = V.length;
	  		for(int i = 0; i < m; i++){
	  			System.out.println(V[i]);
	  		}
	  		System.out.println();
  		}
  	}

  	//__________________________________________________________________________
  	void tiskM(double[][] A){ // tiskne matici
  		int m = A.length;
  		int n = A[0].length;
  		for(int i = 0; i < m; i++){
  			for(int j = 0; j < n; j++){
  				System.out.print(A[i][j]+" ");
  			}
  			System.out.println();
  		}
  		System.out.println();
  	}
}