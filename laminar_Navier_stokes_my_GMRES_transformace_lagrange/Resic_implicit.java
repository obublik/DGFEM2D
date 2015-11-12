import java.io.*;
import java.util.StringTokenizer;
import java.util.Arrays;
import java.util.Comparator;

/* tento program resi Navierovy-Stokesovy rovnice na nestrukturovane
 * siti pomoci nespojite Galerkinovy metody
 */
class Resic_implicit{
    static Casovac c;
    
	public static void main(String[] args) throws IOException{
		double[] hodnoty = nactiHodnoty(); // zde se nacitaji okrajove a ridici hodnoty ze souboru
	    int np      = (int)hodnoty[0];	// pocet bunek
	    int nt      = (int)hodnoty[1];	// pocet bodu
	    double CFL  = hodnoty[2];       // max CFL cislo
	    int rad     = (int)hodnoty[3];  // rad metody
	    double iter  = (int)hodnoty[4]; // pocet iteraci
	    int n_vl    = (int)hodnoty[5];  // pocet vlaken
        int int_iter    = (int)hodnoty[6];  // pocet vnitrnich iteraci
        double tol  = hodnoty[7]; // tolerance pro pridavnou viskozitu
	    double c_IP_max    = hodnoty[8];
        double eps_max    = hodnoty[9];  // maximalni pridavna viskozita
        double eps_0    = hodnoty[10];  // konstantni pridavna viskozita
        double rezid_old = 0;
        Rovnice Eq = new Rovnice();
        int nr = Eq.nr;
        
        int[]  typ = nactiTyp(nt);
        int[][] TP = nactiTP(nt,typ);
	    int[][] TT = nactiTT(nt,typ);
        double[][] PXY = nacti_Body(np);
        double[][] Wpoc = nactiW(nt,nr);

	   	// vytvoreni vypocetnich elementu
	   	Element[] Sit = new Element[nt];
	   	for(int i = 0; i < nt; i++){
	   		int nk = typ[i]; // typ elementu
	   		double[] PX = new double[nk];
	   		double[] PY = new double[nk];
	   		int[] TTe = new int[nk];
	   		int[] TPe = new int[nk];
            int rad_p = rad;
	   		for(int j = 0; j < nk; j++){
	   			PX[j] = PXY[TP[i][j]][0];
	   			PY[j] = PXY[TP[i][j]][1];
	   			TTe[j] = TT[i][j];
	   			TPe[j] = TP[i][j];
//                 if(TTe[j] == -1 && Eq.Re > 0)
// 	   				rad_p = 3;
	   		}
            Sit[i] = new Element(Sit, Eq, PX, PY, TPe, TTe, nk, rad_p, Wpoc[i], tol, c_IP_max, eps_max, eps_0);
	   	}
        
        // nastaveni CFL_max
        double CFL_max = CFL*nastav_CFL(Sit,nt);
        System.out.println("Maximalni CFL cislo je: " + CFL_max);
        
        // nastaveni globalnich indexu
        int sU = 0;
        int sM = 0;
        for(int i = 0; i < nt; i++){ 
            Sit[i].inicializuj();   // dopocet vztahu ktere nebylo mozne delat v konstruktoru
            sU = Sit[i].nastav_globalni_index_U(sU);
            Sit[i].alokuj_sousedy();
            Sit[i].index = i;
        }
        
        // vlastni vypocet
        Gmres G = new Gmres(Sit, sU, 35, 50, 1e-1, n_vl);
        
        c = new Casovac();
        CFL = 0.1;
        double t = 0;
        for(int op = 0; op < iter; op++){
            for(int i = 0; i < nt; i++){
                Sit[i].uloz_Wo();
            }
            
            // nastaveni dt
            double dt = 1000;
            for(int i = 0; i < nt; i++){
                double dt_loc = Sit[i].delta_t(CFL);
                if(dt_loc < dt)
                    dt = dt_loc;
            }
            for(int i = 0; i < nt; i++){
                 Sit[i].dt = dt;
            }
            if(op == 0){
                for(int i = 0; i < nt; i++){
                     Sit[i].dto = dt;
                }
            }
            
            boolean converged = false;
            double[] rezid_loc = new double[int_iter];
            for(int s = 0; s < int_iter; s++){  // vnitrni iterace (Newton)
                if(s == 0){ // limiter
                    for(int i = 0; i < nt; i++){
                        Sit[i].limiter();
                    }
                }
                
                // vytvoreni vlaken, paralelni sestaveni lokalnich matic a plneni globalni matice
                Vlakno[] parallel = new Vlakno[n_vl];
                
                // vlastni vypocet, parallelni beh
                for(int v = 0; v < n_vl; v++){
                    parallel[v] = new Vlakno(v,nt,n_vl,Sit);
                    parallel[v].start();
                }
                try{
                    for(int v = 0; v < n_vl; v++){
                        parallel[v].join();
                }
                }catch(java.lang.InterruptedException e){
                    System.out.println(e);
                }

                // reseni soustavy rovnic
                double[] x = new double[sU];
                converged = G.solve(x);
                
                double rezid = 0;
                if(converged){
                    for(int i = 0; i < x.length; i++){
                        rezid =  rezid + Math.abs(x[i])/nt;
                    }
                    rezid_loc[s] = rezid;
                    System.out.println("   " + s + "-ta vnitrni iterace, reziduum = " + rezid);
                }
                
                if(!converged){
                    for(int i = 0; i < nt; i++){
                        Sit[i].Wo_to_W();
                    }
                    CFL = CFL/2;
                    op = op - 1;
                    System.out.println("GMRES nekonverguje, CFL nastaveno na " + CFL);
                    break;
                }
                
                // ulozeni novych hodnot
                for(int i = 0; i < nt; i++){
                    Sit[i].uloz_reseni(x);
                    Sit[i].lim_zap_hodnot();
                }
                
                double iner_tol = 1e-4;
                if(rezid < iner_tol){
                    break;
                }
            }
            double rezid = 0;
            for(int i = 0; i < nt; i++){
                rezid = rezid + Sit[i].rezid()/nt;
            }
            
            if(converged){
                CFL = CFL + CFL_max/20;
                if(CFL > CFL_max)
                    CFL = CFL_max;
            }
            
            // zastavovaci podminka
            if(CFL < 0.1){
                System.out.println("Metoda nekonverguje.");
                break;
            }
            
            System.out.println(op + "-ta iterace, reziduum = " + rezid + ", CFL = " + CFL);
        }
		ulozWs(Sit,nt,nr);
		ulozWe(Sit,nt,nr);
		uloz_pomocny(Sit,nt);
	}
	
    // vlakno pro paralelni beh
	static class Vlakno extends Thread{
		Element[] Sit;
        int n_start;
        int n_vl;
        int nt;
        Vlakno(int n_start, int nt, int n_vl, Element[] Sit){
        	this.Sit = Sit;
            this.n_start = n_start;
            this.n_vl = n_vl;
            this.nt = nt;
        }

        public void run(){
            // paralelni sestavovani a plneni matic
            for(int i = n_start; i < nt; i = i+n_vl){
                Sit[i].sestav_lokalni_matice();
                Sit[i].vypocti_inverzi_predpodminovace();
            }
      	}
	}
    
    static class Gmres{
        Element[] Sit;
        int n;
        int m;
        int max_it;
        int n_vl;
        double tol;
        double[][] V,H;
        double[] cs, sn, e1, w, r;
        
        Gmres(Element[] Sit, int n, int m, int max_it, double tol, int n_vl){
            this.Sit = Sit;
            this.n = n;
            this.m = m;
            this.max_it = max_it;
            this.tol = tol;
            this.n_vl = n_vl;
            
            // initialize workspace
            V = new double[m+1][n];
            H = new double[m+1][m];
            cs = new double[m];
            sn = new double[m];
            e1 = new double[m+1];
            w = new double[n];
            r = new double[n];
        }
        
        boolean solve(double[] x){
            double norm_r, temp;
            double[] s, y;

            double bnrm2 = 0;
            for(int i = 0; i < Sit.length; i++){
                bnrm2 = bnrm2 + Sit[i].sqr();
            }
            bnrm2 = Math.sqrt(bnrm2);

            if(bnrm2 == 0)
                bnrm2 = 1;
            
            GMRES_reziduum(Sit,x,r,1,n_vl);
            double error = norm(r)/bnrm2;
            if(error < tol)
                return true;
            e1[0] = 1;

            for(int iter = 0; iter < max_it; iter++){          // begin iteration
                GMRES_reziduum(Sit,x,r,1,n_vl);
                norm_r = norm(r);
                for(int j = 0; j < n; j++){
                    V[0][j] = r[j]/norm_r;
                }
                s = nasob_vektor_cislem(e1,norm_r);
                for(int i = 0; i < m; i++){                        // construct orthonormal
                    GMRES_reziduum(Sit,V[i],w,0,n_vl);     // basis using Gram-Schmidt
                    for(int k = 0; k <= i; k++){
                        H[k][i] = nasob_skalar(w,V[k]);
                        for(int j = 0; j < n; j++){
                            w[j] = w[j] - V[k][j]*H[k][i];
                        }
                    }
                    H[i+1][i] = norm(w);
                    for(int j = 0; j < n; j++){
                        V[i+1][j] = w[j]/H[i+1][i];
                    }
                    for(int k = 0; k <= i-1; k++){                              // apply Givens rotation
                        temp  =  cs[k]*H[k][i] + sn[k]*H[k+1][i];
                        H[k+1][i] = -sn[k]*H[k][i] + cs[k]*H[k+1][i];
                        H[k][i]   = temp;
                    }
                    double[] rot = rotmat(H[i][i], H[i+1][i]); // form i-th rotation matrix
                    cs[i] = rot[0];
                    sn[i] = rot[1];
                    temp   = cs[i]*s[i];                            // approximate residual norm
                    s[i+1] = -sn[i]*s[i];
                    s[i]  = temp;
                    H[i][i] = cs[i]*H[i][i] + sn[i]*H[i+1][i];
                    H[i+1][i] = 0;
                    error = Math.abs(s[i+1])/bnrm2;
                    if (error <= tol){                        // update approximation
                       y = lsolve(H,s,i+1);                 // and exit
                       update_solution(V,y,x,i,n_vl);
                       break;
                    }
                 }
                 if(error <= tol)
                      break;

                 y = lsolve(H,s,m);
                 update_solution(V,y,x,m-1,n_vl);
                 GMRES_reziduum(Sit,x,r,1,n_vl);                      // compute residual
                 s[m] = norm(r);
                 error = s[m]/bnrm2;                     // check convergence
                 if(error <= tol)
                     break;
            }
            if(error <= tol)
                return true;
            else
                return false;
        }

        void GMRES_reziduum(Element[] Sit, double[] x, double[] r, int par, int n_vl){
            // vlastni vypocet, parallelni beh
            GMRES_vlakno[] parallel = new GMRES_vlakno[n_vl];
            for(int v = 0; v < n_vl; v++){
                parallel[v] = new GMRES_vlakno(v,n_vl,Sit,x,r,par);
                parallel[v].start();
            }
            try{
                for(int v = 0; v < n_vl; v++){
                    parallel[v].join();
            }
            }catch(java.lang.InterruptedException e){
                System.out.println(e);
            }
        }
        
        void update_solution(double[][] V, double[] y, double[] x, int i, int n_vl){
            // vlastni vypocet, parallelni beh
            US_vlakno[] parallel = new US_vlakno[n_vl];
            for(int v = 0; v < n_vl; v++){
                parallel[v] = new US_vlakno(v,n_vl,V,y,x,i);
                parallel[v].start();
            }
            try{
                for(int v = 0; v < n_vl; v++){
                    parallel[v].join();
            }
            }catch(java.lang.InterruptedException e){
                System.out.println(e);
            }
        }

        double[] rotmat(double a, double b){
        // Compute the Givens rotation matrix parameters for a and b.
            double c, s, temp;
            if(b == 0){
                c = 1;
                s = 0;
            }
            else if(Math.abs(b) > Math.abs(a)){
              temp = a/b;
              s = 1/Math.sqrt(1 + temp*temp);
              c = temp*s;
            }
            else{
              temp = b/a;
              c = 1/Math.sqrt(1 + temp*temp);
              s = temp*c;
            }
            return new double[] {c,s};
        }
        
        double nasob_skalar(double[] a, double[] b){
            int n = a.length;
            double s = 0;
            for(int i = 0; i < n; i++){
                s = s + a[i]*b[i];
            }
            return s;
        }

        double[] nasob_vektor_cislem(double[] a, double b){
            int n = a.length;
            double[] c = new double[n];
            for(int i = 0; i < n; i++){
                c[i]  = a[i]*b;
            }
            return c;
        }

        double norm(double[] a){
            double n = 0;
            int na = a.length;
            for(int i = 0; i < na; i++){
                n = n + a[i]*a[i];
            }
            n = Math.sqrt(n);

            return n;
        }
        
        // Gaussian elimination with partial pivoting
        double[] lsolve(double[][] A, double[] b, int N) {
            // int N  = b.length;
            for (int p = 0; p < N; p++) {
                // find pivot row and swap
                int max = p;
                for (int i = p + 1; i < N; i++) {
                    if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                        max = i;
                    }
                }
                double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
                double t = b[p]; b[p] = b[max]; b[max] = t;

                // pivot within A and b
                for (int i = p + 1; i < N; i++) {
                    double alpha = A[i][p] / A[p][p];
                    b[i] -= alpha * b[p];
                    for (int j = p; j < N; j++) {
                        A[i][j] -= alpha * A[p][j];
                    }
                }
            }

            // back substitution
            double[] x = new double[N];
            for (int i = N - 1; i >= 0; i--) {
                double sum = 0.0;
                for (int j = i + 1; j < N; j++) {
                    sum += A[i][j] * x[j];
                }
                x[i] = (b[i] - sum) / A[i][i];
            }
            return x;
        }
    }
    
     // vlakno pro paralelni beh
	static class GMRES_vlakno extends Thread{
		Element[] Sit;
        int n_start, n_vl, par, nt;
        double[] x, r;
        
        GMRES_vlakno(int n_start, int n_vl, Element[] Sit, double[] x, double[] r, int par){
        	this.Sit = Sit;
            this.n_start = n_start;
            this.n_vl = n_vl;
            this.x = x;
            this.r = r;
            this.par = par;
            nt = Sit.length;
        }

        public void run(){
            // paralelni sestavovani a plneni matic
            for(int i = n_start; i < nt; i = i+n_vl){
                Sit[i].GM_resid(x,r,par);
            }
      	}
	}
    
    // vlakno pro paralelni beh
	static class US_vlakno extends Thread{
        int n_start, n_vl, i, n;
        double[][] V;
        double[] y, x;
        
        US_vlakno(int n_start, int n_vl, double[][] V, double[] y, double[] x, int i){
            this.n_start = n_start;
            this.n_vl = n_vl;
            this.V = V;
            this.y = y;
            this.x = x;
            this.i = i;
            n = x.length;
        }

        public void run(){
            // paralelni sestavovani a plneni matic
            for(int j = n_start; j < n; j = j+n_vl){
               for(int k = 0; k <= i; k++){
                   x[j] = x[j] + V[k][j]*y[k];
               }
            }
      	}
	}
    
     static double nastav_CFL(Element[] Sit, int nt){
         double[] A = new double[nt];
         double A_min = 1e5;
         for(int i = 0; i < nt; i++){
             A[i] = Sit[i].O;
             if(A[i] < A_min)
                 A_min = A[i];
         }
         Arrays.sort(A);
         double A_med = A[(int)(0.5*nt)];
//          return Math.sqrt(A_med/A_min);
         return 1.0;
     }
    
// Pomocne metody_______________________________________________________________
//______________________________________________________________________________
	static double[] nactiHodnoty() throws IOException{
	    int m = 11;
	    double [] u = new double[m];
	    File soubor = new File("parameters_solver.txt");
	    if(soubor.exists()){
	        FileReader fr = new FileReader("parameters_solver.txt");
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

	static int[][] nactiTT(int nt, int[] typ) throws IOException{
        int[][] TT = new int[nt][];
        File soubor = new File("TT.txt");
        FileReader fr = new FileReader("TT.txt");
        BufferedReader in = new BufferedReader(fr);
        String c;
        String radka;
        for(int i = 0; i < nt; i++){
            TT[i] = new int[typ[i]];
            radka = in.readLine();
            StringTokenizer st = new StringTokenizer(radka);
            for(int k = 0; k < typ[i]; k++){
                c = st.nextToken();
                TT[i][k] = Integer.valueOf(c).intValue();
            }
        }
        fr.close();
        return(TT);
    }

	static int[][] nactiTP(int nt, int[] typ) throws IOException{
        int[][] TP = new int[nt][];
        File soubor = new File("TP.txt");
        FileReader fr = new FileReader("TP.txt");
        BufferedReader in = new BufferedReader(fr);
        String c;
        String radka;
        for(int i = 0; i < nt; i++){
            TP[i] = new int[typ[i]];
            radka = in.readLine();
            StringTokenizer st = new StringTokenizer(radka);
            for(int k = 0; k < typ[i]; k++){
                c = st.nextToken();
                TP[i][k] = Integer.valueOf(c).intValue();
            }
        }
        fr.close();
        return(TP);
    }
    
    static int[] nactiTyp(int nt) throws IOException{
        int[] typ = new int[nt];
        File soubor = new File("typ.txt");
        FileReader fr = new FileReader("typ.txt");
        BufferedReader in = new BufferedReader(fr);
        String c;
        String radka;
        for(int i = 0; i < nt; i++){
            radka = in.readLine();
            StringTokenizer st = new StringTokenizer(radka);
            c = st.nextToken();
            typ[i] = Integer.valueOf(c).intValue();
        }
        fr.close();
        return(typ);
    }

	static double[][] nactiW(int nt, int nr) throws IOException{
        double[][] W = new double[nt][nr];
        File soubor = new File("W.txt");
        FileReader fr = new FileReader("W.txt");
        BufferedReader in = new BufferedReader(fr);
        String c;
        String radka;
        for(int i = 0; i < nt; i++){
            radka = in.readLine();
            StringTokenizer st = new StringTokenizer(radka);
            for(int j = 0; j < nr; j++){
                c = st.nextToken();
                W[i][j] = Double.valueOf(c).doubleValue();
            }
        }
        fr.close();
        return(W);
    }

    static double[][] nacti_Body(int np) throws IOException{
        double[][] PXY = new double[np][2];
        File soubor = new File("PXY.txt");
        FileReader fr = new FileReader("PXY.txt");
        BufferedReader in = new BufferedReader(fr);
        String c1;
        String c2;
        String radka;
        for(int i = 0; i < np; i++){
            radka = in.readLine();
            StringTokenizer st = new StringTokenizer(radka);
            c1 = st.nextToken();
            c2 = st.nextToken();
            PXY[i][0] = Double.valueOf(c1).doubleValue();
            PXY[i][1] = Double.valueOf(c2).doubleValue();
        }
        fr.close();
        return(PXY);
    }

    static void ulozWs(Element[] Sit, int nt, int nr) throws IOException{
        File soubor = new File("W.txt");
        soubor.createNewFile();
        FileWriter fw = new FileWriter("W.txt");
        BufferedWriter out = new BufferedWriter(fw);
        String radka;
        for(int i = 0; i < nt; i++){
        	double[] Ws = Sit[i].spocti_Ws();
            radka =  String.valueOf(Ws[0]);
            for(int j = 1; j < nr; j++){
                radka =  radka+" "+ String.valueOf(Ws[j]);
            }
            out.write(radka);
            out.newLine();
        }
        out.close();
    }

    static void ulozWe(Element[] Sit, int nt, int nr) throws IOException{
        File soubor = new File("We.txt");
        soubor.createNewFile();
        FileWriter fw = new FileWriter("We.txt");
        BufferedWriter out = new BufferedWriter(fw);
        String radka;
        for(int i = 0; i < nt; i++){
        	int ne = Sit[i].ne;
        	for(int m = 0; m < ne; m++){
				double[][] W = Sit[i].W;
                radka =  String.valueOf(W[m][0]);
                for(int j = 1; j < nr; j++){
                    radka =  radka+" "+ String.valueOf(W[m][j]);
                }
            	out.write(radka);
            	out.newLine();
            }
        }
        out.close();

        soubor = new File("rad_elementu.txt");
        soubor.createNewFile();
        fw = new FileWriter("rad_elementu.txt");
        out = new BufferedWriter(fw);
        for(int i = 0; i < nt; i++){
        	int ne = Sit[i].rad;
        	radka =  String.valueOf(ne);
        	out.write(radka);
        	out.newLine();
        }
        out.close();
    }
	
	static void uloz_pomocny(Element[] Sit, int nt) throws IOException{
        File soubor = new File("pom.txt");
        soubor.createNewFile();
        FileWriter fw = new FileWriter("pom.txt");
        BufferedWriter out = new BufferedWriter(fw);
        String radka;
        for(int i = 0; i < nt; i++){
            radka = String.valueOf(Sit[i].pom);
            out.write(radka);
            out.newLine();
        }
        out.close();
    }
    
    //__________________________________________________________________________
    
    
  	static void tiskV(double[] V){ // tiskne vektor
  		if(V == null)
  			System.out.println(V);
  		else{
	  		int m = V.length;
	  		for(int i = 0; i < m; i++){
	  			System.out.print(V[i]+" ");
	  		}
	  		System.out.println();
  		}
  	}

  	//__________________________________________________________________________
  	static void tiskM(double[][] A){ // tiskne matici
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

class Casovac{
    long startTime;
    Casovac(){
    }
    
    void s(){
        startTime = System.currentTimeMillis();
    }
    
    void k(){
        System.out.println( System.currentTimeMillis() - startTime);
    }
}









