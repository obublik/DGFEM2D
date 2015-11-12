import java.io.*;
import java.util.StringTokenizer;

class Rovnice{
    
    int nr = 4; // pocet rovnic
    int difuze;
    int zdroj;
    
    // okrajove podminky
    double pi; // tlak na vstupu
    double ri; // hustota na vstupu
    double alfa; // uhel nabehu
    double po; // vystupni tlak
    double kapa; // Poissonova konstanta
    double Re; // Reynoldsovo cislo
    double Pr; // Prandtlovo cislo
    double[] W_sup;
    
    // tolerance pro hodnotu hustoty
	double R_eps = 1e-2;
	double p_eps = 1e-2;
    
    Rovnice()  throws IOException{
        double[] hodnoty = nactiHodnoty(); // zde se nacitaji okrajove a ridici hodnoty ze souboru
	    pi = hodnoty[0]; 	// tlak na vstupu
	    ri = hodnoty[1];   // hustota na vstupu
	    alfa = hodnoty[2];	// uhel nabehu
	    po = hodnoty[3];		// vystupni tlak
        kapa = hodnoty[4];	// Poissonova konstanta
	    Re = hodnoty[5];	// Reynoldsovo cislo
        Pr = hodnoty[6];	// Prandtlovo cislo
	    int vstup_typ = (int)hodnoty[7]; // typ vstupni podminky
	    double[] W_sup2 = {hodnoty[8], hodnoty[9], hodnoty[10], hodnoty[11]};
        this.W_sup = W_sup2;
        
        difuze = 1;
        if(Re == -1)
            difuze = 0;
        
        int zdroj = 0;
        
	    if(vstup_typ == 0)
	    	W_sup[0] = -1;
    }
    
     //  nevazky tok stenou _____________________________________________________
	double[] Fn(double WL[], double WR[], double nx, double ny, int TT){

		double[] fn = new double[4];
        double[] fcL, fcR;
        double p, pL, pR;

		switch(TT){
        	case(-1): // stena
                p = tlak(WL);
        		fn[0] = 0;
        		fn[1] = p*nx;
        		fn[2] = p*ny;
        		fn[3] = 0;
        	break;

        	case(-2): // vstup
                fn = fp(WR,nx,ny);
        	break;
            
        	case(-3): // vystup
                fn = fp(WR,nx,ny);
        	break;
			
			case(-4): // nevazka stena
                p = tlak(WL);
        		fn[0] = 0;
        		fn[1] = p*nx;
        		fn[2] = p*ny;
        		fn[3] = 0;
        	break;
			
        	default: // vnitrni stena
                    pL = tlak(WL);
                    pR = tlak(WR);
                    double aL = Math.sqrt(kapa*pL/WL[0]);
                    double aR = Math.sqrt(kapa*pR/WR[0]);
                    double VnL = WL[1]/WL[0]*nx + WL[2]/WL[0]*ny;
                    double VnR = WR[1]/WR[0]*nx + WR[2]/WR[0]*ny;
                    fcL = fp(WL,nx,ny);
                    fcR = fp(WR,nx,ny);
                    double S = Math.max(Math.abs(VnL)+aL,Math.abs(VnR)+aR);
                    for(int j = 0; j < 4; j++){
                        fn[j] = (fcL[j] + fcR[j] - S*(WR[j]-WL[j]))/2;
                    }

        	break;
        }
	    return fn;
	}
	
	double[] fp(double[] W, double nx, double ny){
		double[] f = new double[4];
	    double V = W[1]/W[0]*nx + W[2]/W[0]*ny;
        double p = tlak(W);
	    f[0] = W[0]*V;
	    f[1] = W[1]*V + p*nx;
	    f[2] = W[2]*V + p*ny;
	    f[3] = (W[3]+p)*V;

	    return f;
	}
    
    // vazky tok stenou ________________________________________________________
	double[] Fvn(double WL[], double WR[], double dWL[], double dWR[], double nx, double ny, int TT){
		double[] fvnL = new double[4];
		double[] fvnR = new double[4];
		double[] fvn = new double[4];
		
		fvnL = fpv(WL,dWL,nx,ny);
		fvnR = fpv(WR,dWR,nx,ny);

        for(int j = 0; j < 4; j++){
	    	fvn[j] = (fvnL[j] + fvnR[j])/2;
	    }
        
        if(TT < 0){
            if(TT == -1)
                fvn[3] = 0;
            else
                for(int j = 0; j < 4; j++){
                    fvn[j] = 0;
                }
        }
        
	    return fvn;
	}
    
    double[] fpv(double[] W, double[] dW, double nx, double ny){
        double[] fvn = new double[4];
        double lam = -2./3; // Stokesuv vztah
        double r = W[0];
        double u = W[1]/r;
        double v = W[2]/r;
        double rx = dW[0];
        double ry = dW[4];
        double ux = 1/r*(dW[1] - rx*u);
        double uy = 1/r*(dW[5] - ry*u);
        double vx = 1/r*(dW[2] - rx*v);
        double vy = 1/r*(dW[6] - ry*v);
        double Ex = dW[3];
        double Ey = dW[7];
        double p = tlak(W);
        double px = (kapa-1)*(Ex-0.5*rx*(u*u+v*v)-r*(u*ux+v*vx));
        double py = (kapa-1)*(Ey-0.5*ry*(u*u+v*v)-r*(u*uy+v*vy));
        double prx = (r*px - p*rx)/(r*r);
        double pry = (r*py - p*ry)/(r*r);
        
        double txx = 2*ux + lam*(ux+vy);
        double txy = vx+uy;
        double tyy = 2*vy + lam*(ux+vy);
        
        fvn[0] = 0;
        fvn[1] = 1/Re*txx*nx + 1/Re*txy*ny;
        fvn[2] = 1/Re*txy*nx + 1/Re*tyy*ny;
        fvn[3] = (1/Re*(u*txx+ v*txy+ kapa/(kapa-1)/Pr*prx))*nx + (1/Re*(u*txy + v*tyy + kapa/(kapa-1)/Pr*pry))*ny;
        
        return fvn;
    }
	
	double tlak(double[] W){
		double p = (kapa-1)*(W[3]-(W[1]*W[1]+W[2]*W[2])/(2*W[0]));
		if(p < p_eps)
			p = p_eps;
		return p;
	}
    
    double[] hodnota_WR(double[] WL, double nx, double ny, int TT){
        double[] WR = new double[nr];
        double p = tlak(WL);
        switch(TT){
            case(-1): // stena
                if(Re != -1){
                    WR[0] = WL[0];
                    WR[1] = 0;
                    WR[2] = 0;
                    WR[3] = p/(kapa-1);
                }
                else{
                    WR[0] = WL[0];
                    WR[1] = WL[1];
                    WR[2] = WL[2];
                    WR[3] = WL[3];
                }
            break;

            case(-2):
                // osetreni NaN
//         		if(pi/p <= 1){
//         			p = pi;
//         		}
                if(W_sup[0] == -1){ // subsonicky vstup
                    double Minl;
                    if(p > pi){
                        Minl = 0;
                    }
                    else{
                        Minl = Math.sqrt((2/(kapa-1))*(-1 + Math.pow(pi/p,(kapa-1)/kapa)));
                    }
                    double Rinl = ri*Math.pow((1+((kapa-1)/2)*Minl*Minl),1/(1-kapa));
                    double Vinl = Minl*Math.sqrt((kapa*p)/Rinl);
                    double uinl = -Vinl*nx;
                    double vinl = -Vinl*ny;
                    double Einl = p/(kapa-1)+0.5*Rinl*Vinl*Vinl;
                
                	WR[0] = Rinl;
                	WR[1] = Rinl*uinl;
                	WR[2] = Rinl*vinl;
                	WR[3] = Einl;
                }
                else{ // supersonicky vstup
                	WR[0] = W_sup[0];
                    WR[1] = W_sup[1];
                    WR[2] = W_sup[2];
                    WR[3] = W_sup[3];
                }
            break;

            case(-3):
                double ro = WL[0];
                double uo = WL[1]/WL[0];
                double vo = WL[2]/WL[0];
                double ao = Math.sqrt(kapa*p/ro);
        		double Mo = Math.sqrt(uo*uo + vo*vo)/ao;
                double Eo, pout;
                if(Mo < 1){ // subsonicky vystup
            		pout = po;
                }
        		else{ // supersonicky vystup
        			pout = p;
				}
				Eo = pout/(kapa-1)+ro*(uo*uo+vo*vo)/2;
				WR[0] = ro;
                WR[1] = ro*uo;
                WR[2] = ro*vo;
                WR[3] = Eo;
            break;
            
            case(-4): // nevazka stena
                WR[0] = WL[0];
                WR[1] = WL[1];
                WR[2] = WL[2];
                WR[3] = WL[3];
            break;
        }
        return WR;
    }
    
    double[] P(double[] W, double[] dW){ // zdrojovy clen
        return null;
    }
    
    
    // Pomocne metody_______________________________________________________________
//______________________________________________________________________________
	static double[] nactiHodnoty() throws IOException{
	    int m = 12;
	    double [] u = new double[m];
	    File soubor = new File("parameters_equations.txt");
	    if(soubor.exists()){
	        FileReader fr = new FileReader("parameters_equations.txt");
	        BufferedReader in = new BufferedReader(fr);
	        for(int i = 0; i < m; i++){
	            String cislo;
	            cislo = in.readLine();
	            u[i] = Double.valueOf(cislo).doubleValue();
	        }
	        fr.close();
	    }
	    return(u);
	}
}
