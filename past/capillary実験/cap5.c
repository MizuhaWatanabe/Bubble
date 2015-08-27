# include<stdio.h>
# include<math.h>

void RungeKutta(double, double, double, double, double, double *, double, double, double, double *);
double differential(double, double, double, double, double, double, double, double);
double Fsigma(double, double, double);
double Fmu(double, double, double);
double Fg(double, double, double, double);
double Ff(double, double, double, double);

double r = 2.4876e-04;
double rerror = 5.0e-05;
double g = 9.8;
double sigma = 0.070;
double rho = 1000;
double mu = 0.001;
double thetae = 0 * M_PI / 180;
double xi = 8.18;

double rplus, rminus,
       t, dt,
       l, lplus, lminus,
       dldt, dldtplus, dldtminus,
       linf, lerror,
       fsigma, fmu, fg, ff,
       Re, Fr, We;


int main(){
    FILE *fp;
    char *fname = "cap5.csv";
    char *l1 = "times[s]";
    char *l2 = "height[m]";
    char *l3 = "height_up[m]";
    char *l4 = "height_down[m]";
    char *l5 = "velocity[m/s]";
    char *l6 = "capillary_force[N]";
    char *l7 = "viscous_force[N]";
    char *l8 = "gravity_force[N]";
    char *l9 = "friction_force[N]";
    char *l10 = "Reynolds_Number";
    char *l11 = "Fluid_Number";
    char *l12 = "Weber_Number";

    fp = fopen(fname, "w");
    if(fp == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }
    fprintf(fp, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12);


    rplus = r + rerror;
    rminus = r - rerror;
    t = 0;
    dt = 0.001;
    l = lplus = lminus = 0;
    dldt = differential(xi, mu, r, rho, l, g, sigma, thetae);
    linf = 2 * sigma * cos(thetae) / (rho * g * r);
    fsigma = Fsigma(r, sigma, thetae);
    fmu = Fmu(mu, dldt, l);
    fg = Fg(rho, r, l, g);
    ff = Ff(xi, r, rho, dldt);
    Re = 2 * rho * r * dldt / mu;
    Fr = dldt / pow(r * g, 0.5);
    We = rho * r * dldt * dldt / sigma;

    int i = 0;
    while(1){
        i += 1;

        fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", t, l, lplus, lminus, dldt, fsigma, fmu, fg, ff, Re, Fr, We);

        RungeKutta(dt, xi, mu, r, rho, &l, g, sigma, thetae, &dldt);
        RungeKutta(dt, xi, mu, rplus, rho, &lplus, g, sigma, thetae, &dldtplus);
        RungeKutta(dt, xi, mu, rminus, rho, &lminus, g, sigma, thetae, &dldtminus);

        fsigma = Fsigma(r, sigma, thetae);
        fmu = Fmu(mu, dldt, l);
        fg = Fg(rho, r, l, g);
        ff = Ff(xi, r, rho, dldt);
        Re = 2 * rho * r * dldt / mu;
        Fr = dldt / pow(r * g, 0.5);
        We = rho * r * dldt * dldt / sigma;
        t = t + dt;

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


void RungeKutta(double dt, double xi, double mu, double r, double rho, double *l, double g, double sigma, double thetae, double *dldt){
    double k, k1, k2, k3, k4;
    k1 = dt * differential(xi, mu, r, rho, *l, g, sigma, thetae);
    k2 = dt * differential(xi, mu, r, rho, *l + k1 * 0.5, g, sigma, thetae);
    k3 = dt * differential(xi, mu, r, rho, *l + k2 * 0.5, g, sigma, thetae);
    k4 = dt * differential(xi, mu, r, rho, *l + k3, g, sigma, thetae);
    k = (k1 + k2 * 2 + k3 * 2 + k4) / 6;
    *l = *l + k;
    *dldt = k / dt;
}


double differential(double xi, double mu, double r, double rho, double l, double g, double sigma, double thetae){
    double viscous;
    double gravity;
    double tension;
    double friction;
    viscous = 8 * mu * l / (r * r * rho);
    gravity = g * l;
    tension = 2 * sigma * cos(thetae) / (r * rho);
    friction = (1 / xi) * (-1 * viscous + pow(viscous * viscous - 2 * xi * (gravity - tension), 0.5));
    return friction;
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


double Ff(double xi, double r, double rho, double dldt){
    double Ff;
    Ff = -1 * xi * M_PI * r * r * rho * dldt * dldt / 2;
    return Ff;
}
