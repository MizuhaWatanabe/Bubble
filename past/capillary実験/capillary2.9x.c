# include<stdio.h>
# include<math.h>

void RungeKutta(double, double *, double, double, double, double, double, double, double *, double *);
double differential(double, double, double, double, double, double, double, double);
double Fsigma(double, double, double);
double Fi(double, double, double, double, double);
double Fmu(double, double, double);
double Fg(double, double, double, double);

double r = 2.4876e-04;
double rerror = 5e-05;
double g = 9.8;
double sigma = 0.071;
double rho = 1000;
double mu = 0.001;
double thetae = 0 * M_PI / 180;     /* thetae has rad dimention */

double rplus, rminus,
t, dt,
l, lplus, lminus,
dldt, dldtplus, dldtminus,
d2ldt2,
linf, lerror,
L, Lplus, Lminus,
V, Vplus, Vminus,
dVdt, dVdtplus, dVdtminus,
fsigma, fi, fmu, fg,
Re;


int main(){
    FILE *fp;
    char *fname = "capillary2.9x.csv";
    char *l1 = "times[s]";
    char *l2 = "height[m]";
    char *l3 = "height_up[m]";
    char *l4 = "height_down[m]";
    char *l5 = "velocity[m/s]";
    char *l6 = "capillary_force[N]";
    char *l7 = "inertia_force[N]";
    char *l8 = "viscous_force[N]";
    char *l9 = "gravity_force[N]";
    char *l10 = "Reynolds_Number";

    fp = fopen(fname, "w");
    if(fp == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }
    fprintf(fp, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", l1, l2, l3, l4, l5, l6, l7, l8, l9, l10);


    rplus = r + rerror;
    rminus = r - rerror;
    t = 0;
    dt = 0.001;
    l = lplus = lminus = 0;
    dldt = pow(2 * sigma * cos(thetae) / (rho * r), 0.5);
    dldtplus = pow(2 * sigma * cos(thetae) / (rho * rplus), 0.5);
    dldtminus = pow(2 * sigma * cos(thetae) / (rho * rminus), 0.5);
    L = Lplus = Lminus = 0;
    V = Vplus = Vminus = 0;
    linf = 2 * sigma * cos(thetae) / (rho * g * r);

    int i = 0;
    while(1){
        i += 1;
        RungeKutta(dt, &V, sigma, rho, r, thetae, mu, g, &L, &dVdt);
        RungeKutta(dt, &Vplus, sigma, rho, rplus, thetae, mu, g, &Lplus, &dVdtplus);
        RungeKutta(dt, &Vminus, sigma, rho, rminus, thetae, mu, g, &Lminus, &dVdtminus);

        fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", t, l, lplus, lminus, dldt, fsigma, fi, fmu, fg, Re);

        t = t + dt;
        l = pow(2 * L, 0.5);
        lplus = pow(2 * Lplus, 0.5);
        lminus = pow(2 * Lminus, 0.5);
        dldt = V / l;
        d2ldt2 = (dVdt - (dldt * dldt)) / l;
        fsigma = Fsigma(r, sigma, thetae);
        fi = Fi(rho, r, l, dldt, d2ldt2);
        fmu = Fmu(mu, dldt, l);
        fg = Fg(rho, r, l, g);
        Re = 2 * sigma * r * dldt / mu;
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

void RungeKutta(double dt, double *V, double sigma, double rho, double r, double thetae, double mu, double g, double *L, double *dVdt){
    double k, k1, k2, k3, k4,
    m, m1, m2, m3, m4;
    k1 = dt * *V;
    m1 = dt * differential(sigma, rho, r, thetae, mu, *V, g, *L);
    k2 = dt * (*V + m1 * 0.5);
    m2 = dt * differential(sigma, rho, r, thetae, mu, *V + m1 * 0.5, g, *L + k1 * 0.5);
    k3 = dt * (*V + m2 * 0.5);
    m3 = dt * differential(sigma, rho, r, thetae, mu, *V + m2 * 0.5, g, *L + k2 * 0.5);
    k4 = dt * (*V + m3);
    m4 = dt * differential(sigma, rho, r, thetae, mu, *V + m3, g, *L + k3);
    k = (k1 + k2 * 2 + k3 * 2 + k4) / 6;
    m = (m1 + m2 * 2 + m3 * 2 + m4) / 6;
    *L = *L + k;
    *V = *V + m;
    *dVdt = m / dt;
}

double differential(double sigma, double rho, double r, double thetae, double mu, double V, double g, double L){
    double tension;
    double viscous;
    double gravity;
    double inertia;
    tension = 2 * sigma * cos(thetae) / (rho * r);
    viscous = 8 * mu * V / (rho * r * r);
    gravity = g * pow(2 * L, 0.5);
    inertia = tension - viscous - gravity;
    return inertia;
}

double Fsigma(double r, double sigma, double thetae){
    double Fsigma;
    Fsigma = 2 * M_PI *r * sigma * cos(thetae);
    return Fsigma;
}

double Fi(double rho, double r, double l, double dldt, double d2ldt2){
    double Fi;
    Fi = -1 * M_PI * rho * r * r * (dldt * dldt + l * d2ldt2);
    return Fi;
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
