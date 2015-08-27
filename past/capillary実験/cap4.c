# include<stdio.h>
# include<math.h>

void RungeKutta(double, double, double, double, double *, double, double, double *, double *, double, double);
double Ldifferential(double, double, double, double, double, double, double);
double thetadifferential(double, double, double, double, double, double, double, double);
double Fsigma(double, double, double);
double Fmu(double, double, double);
double Fg(double, double, double, double);

double r = 2.5e-04;
double rerror = 5.0e-05;
double g = 9.8;
double sigma = 0.071;
double rho = 1000;
double mu = 0.001;
double thetae = 10 * M_PI / 180;
double epsilon = 15;

double rplus, rminus,
       t, dt,
       l, lplus, lminus,
       dldt,
       linf, lerror,
       L, Lplus, Lminus,
       dLdt, dLdtplus, dLdtminus,
       theta, thetaplus, thetaminus,
       fsigma, fmu, fg,
       Re;


int main(){
    FILE *fp;
    char *fname = "capillary4.csv";
    char *l1 = "times[s]";
    char *l2 = "height[m]";
    char *l3 = "height_up[m]";
    char *l4 = "height_down[m]";
    char *l5 = "velocity[m/s]";
    char *l6 = "theta[rad]";
    char *l7 = "capillary_force[N]";
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
    L = Lplus = Lminus = 0;
    theta = thetaplus = thetaminus = 89 * M_PI / 180;
    linf = 2 * sigma / (rho * g * r);

    int i = 0;
    while(1){
        i += 1;
        RungeKutta(dt, sigma, rho, r, &theta, mu, g, &L, &dLdt, epsilon, thetae);
        RungeKutta(dt, sigma, rho, rplus, &thetaplus, mu, g, &Lplus, &dLdtplus, epsilon, thetae);
        RungeKutta(dt, sigma, rho, rminus, &thetaminus, mu, g, &Lminus, &dLdtminus, epsilon, thetae);

        fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", t, l, lplus, lminus, dldt, theta, fsigma, fmu, fg, Re);

        t = t + dt;
        l = pow(2 * L, 0.5);
        lplus = pow(2 * Lplus, 0.5);
        lminus = pow(2 * Lminus, 0.5);
        dldt = dLdt / l;
        fsigma = Fsigma(r, sigma, theta);
        fmu = Fmu(mu, dldt, l);
        fg = Fg(rho, r, l, g);
        Re = 2 * sigma * r * dldt / mu;

        if(i > 9999){
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


void RungeKutta(double dt, double sigma, double rho, double r, double *theta, double mu, double g, double *L, double *dLdt, double epsilon, double thetae){
    double k, k1, k2, k3, k4,
           m, m1, m2, m3, m4;
    k1 = dt * Ldifferential(sigma, rho, r, *theta, mu, g, *L);
    m1 = dt * thetadifferential(r, epsilon, *theta, thetae, mu, rho, g, sigma);
    k2 = dt * Ldifferential(sigma, rho, r, *theta + m1 * 0.5, mu, g, *L + k1 * 0.5);
    m2 = dt * thetadifferential(r, epsilon, *theta + m1 * 0.5, thetae, mu, rho, g, sigma);
    k3 = dt * Ldifferential(sigma, rho, r, *theta + m2 * 0.5, mu, g, *L + k2 * 0.5);
    m3 = dt * thetadifferential(r, epsilon, *theta, thetae + m2 * 0.5, mu, rho, g, sigma);
    k4 = dt * Ldifferential(sigma, rho, r, *theta + m3, mu, g, *L + k3);
    m1 = dt * thetadifferential(r, epsilon, *theta + m3, thetae, mu, rho, g, sigma);
    k = (k1 + k2 * 2 + k3 * 2 + k4) / 6;
    m = (m1 + m2 * 2 + m3 * 2 + m4) / 6;
    *L = *L + k;
    *dLdt = k / dt;
    *theta = *theta + m;
}


double Ldifferential(double sigma, double rho, double r, double theta, double mu, double g, double L){
    double tension;
    double gravity;
    double viscous;
    tension = r * sigma * cos(theta) / (4 * mu);
    gravity = r * r * rho * g * pow(2 * L, 0.5) / (8 * mu);
    viscous = tension - gravity;
    return viscous;
}

double thetadifferential(double r, double epsilon, double theta, double thetae, double mu, double rho, double g, double sigma){
    double header, theta1, theta2, theta3, thetadif;
    header = -4 * r / (3 * epsilon);
    theta1 = sin(theta) * (cos(thetae) - cos(theta));
    theta2 = (r * r * rho * g / (8 * mu)) + sigma * tan(theta) * (cos(thetae) - cos(theta)) / (3 * epsilon * mu);
    theta3 = sigma * (pow(cos(theta), -2) + tan(theta) * sin(theta)) / (3 * epsilon * mu);
    thetadif = header * theta1 * theta2 /(theta3 + tan(theta) * theta2);
    return thetadif;
}


double Fsigma(double r, double sigma, double theta){
    double Fsigma;
    Fsigma = 2 * M_PI *r * sigma * cos(theta);
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
