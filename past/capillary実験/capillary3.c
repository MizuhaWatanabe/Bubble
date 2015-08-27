# include<stdio.h>
# include<math.h>

void RungeKutta(double, double, double, double, double, double, double, double *, double *, double, double);
double differential(double, double, double, double, double, double, double, double, double);
double Fsigma(double, double, double);
double Fmu(double, double, double);
double Fg(double, double, double, double);

double r = 5.00e-04;
double rerror = 5.0e-05;
double g = 9.8;
double sigma = 0.071;
double rho = 1000;
double mu = 0.001;
double thetae = 0 * M_PI / 180;
double epsilon = 1;

double rplus, rminus,
       t, dt,
       l, lplus, lminus,
       dldt,
       linf, linfplus, linfminus,
       lerror,
       L, Lplus, Lminus,
       dLdt, dLdtplus, dLdtminus,
       fsigma, fmu, fg,
       Re, Fr, We;


int main(){
    FILE *fp;
    char *fname = "capillary3.csv";
    char *l1 = "times[s]";
    char *l2 = "height[m]";
    char *l3 = "height_up[m]";
    char *l4 = "height_down[m]";
    char *l5 = "velocity[m/s]";
    char *l6 = "capillary_force[N]";
    char *l7 = "viscous_force[N]";
    char *l8 = "gravity_force[N]";
    char *l9 = "Reynolds_Number";
    char *l10 = "Fluid_Number";
    char *l11 = "Weber_Number";

    fp = fopen(fname, "w");
    if(fp == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }
    fprintf(fp, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11);


    rplus = r + rerror;
    rminus = r - rerror;
    t = 0;
    dt = 0.001;
    l = lplus = lminus = 0;
    linf = 2 * sigma * cos(thetae) / (rho * g * r);
    linfplus = 2 * sigma * cos(thetae) / (rho * g * rplus);
    linfminus = 2 * sigma * cos(thetae) / (rho * g * rminus);
    L = 0.5 * ((epsilon - 1) * linf) * ((epsilon - 1) * linf);
    Lplus = 0.5 * ((epsilon - 1) * linfplus) * ((epsilon - 1) * linfplus);
    Lminus = 0.5 * ((epsilon - 1) * linfminus) * ((epsilon - 1) * linfminus);

    int i = 0;
    while(1){
        i += 1;
        RungeKutta(dt, sigma, rho, r, thetae, mu, g, &L, &dLdt, epsilon, linf);
        RungeKutta(dt, sigma, rho, rplus, thetae, mu, g, &Lplus, &dLdtplus, epsilon, linfplus);
        RungeKutta(dt, sigma, rho, rminus, thetae, mu, g, &Lminus, &dLdtminus, epsilon, linfminus);

        fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", t, l, lplus, lminus, dldt, fsigma, fmu, fg, Re, Fr, We);

        t = t + dt;
        l = pow(2 * L, 0.5) - (epsilon - 1) * linf;
        lplus = pow(2 * Lplus, 0.5) - (epsilon - 1) * linfplus;
        lminus = pow(2 * Lminus, 0.5) - (epsilon - 1) * linfminus;
        dldt = dLdt / (l + (epsilon - 1) * linf);
        fsigma = Fsigma(r, sigma, thetae);
        fmu = Fmu(mu, dldt, l);
        fg = Fg(rho, r, l, g);
        Re = 2 * sigma * r * dldt / mu;
        Fr = dldt / pow(r * g, 0.5);
        We = rho * r * dldt * dldt / sigma;
        lerror = fabs(linf - l) / l;

        if(lerror < 0.001 || i > 9999){
            fprintf(fp, "\n");
            break;
        }
    }

    fclose(fp);
    printf("%i回の計算で収束しました.\n", i);
    printf("最終上昇水位[m]:%e\n", l);
    printf("%sファイルへの書き込みが完了しました.\n", fname);
    return 0;
}


void RungeKutta(double dt, double sigma, double rho, double r, double thetae, double mu, double g, double *L, double *dLdt, double epsilon, double linf){
    double k, k1, k2, k3, k4;
    k1 = dt * differential(sigma, rho, r, thetae, mu, g, *L, epsilon, linf);
    k2 = dt * differential(sigma, rho, r, thetae, mu, g, *L + k1 * 0.5, epsilon, linf);
    k3 = dt * differential(sigma, rho, r, thetae, mu, g, *L + k2 * 0.5, epsilon, linf);
    k4 = dt * differential(sigma, rho, r, thetae, mu, g, *L + k3, epsilon, linf);
    k = (k1 + k2 * 2 + k3 * 2 + k4) / 6;
    *L = *L + k;
    *dLdt = k / dt;
}


double differential(double sigma, double rho, double r, double thetae, double mu, double g, double L, double epsilon, double linf){
    double tension;
    double gravity;
    double viscous;
    tension = r * sigma * cos(thetae) / (4 * mu);
    gravity = r * r * rho * g * (pow(2 * L, 0.5) - (epsilon - 1) * linf) / (8 * mu);
    viscous = tension - gravity;
    return viscous;
}


double Fsigma(double r, double sigma, double thetae){
    double Fsigma;
    Fsigma = 2 * M_PI *r * sigma * cos(thetae);
    return Fsigma;
}


double Fmu(double mu, double dldt, double l){
    double Fmu;
    Fmu = -8 * M_PI * mu * dldt * l;
    return Fmu;
}


double Fg(double rho, double r, double l, double g){
    double Fg;
    Fg = -1 * M_PI * rho * r * r * l * g;
    return Fg;
}
