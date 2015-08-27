#include<stdio.h>
#include<math.h>

int main(void){
    double sigma = 0.071;
    double r = 0.00025;
    double g = 9.8;
    double rho = 1000;
    double mu = 0.001;
    double thetae = 0 * M_PI / 180;              /* theta has rad dimention */
    double A = 15;                   /* A is 15~20 according to the paper*/

    double linf = (2 * sigma) / (r * g * rho);
    double k = (pow(r, 2) * g * rho) / (8 * mu);
    double eps = sigma / (6 * A * mu);

    double l, dl, t, dt, L, dL, dL0, dL1, dL2, dL3, y1, y2, y3, theta, dtheta, dtheta0, dtheta1, dtheta2, dtheta3, z1, z2, z3;

    double dLdt(double k, double linf, double theta, double L);
    double dthetadt(double k, double linf, double eps, double theta, double thetae);

    FILE *fp;
    char *fname = "capillary2.csv";
    char *l1 = "time[s]";
    char *l2 = "height[m]";
    char *l3 = "dl[m/s]";
    char *l4 = "angle[rad]";
    char *l5 = "dthetadt[rad/s]";

    fp = fopen(fname, "w");
    if(fp == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }

    L = 0;
    l = 0;
    theta = 63.5 * M_PI / 180;
    t = 0;
    dt = 0.002;

    fprintf(fp, "%s,%s,%s,%s,%s\n", l1, l2, l3, l4, l5);

    int i = 0;
    while( i < 10000){
        i += 1;
        dL0 = dt * dLdt(k, linf, theta, L);
        dtheta0 = dt * dthetadt(k, linf, eps, theta, thetae);
        y1 = L + (dL0 * 0.5);
        z1 = theta + (dtheta0 * 0.5);
        dL1 = dt * dLdt(k, linf, z1, y1);
        dtheta1 = dt * dthetadt(k, linf, eps, z1, thetae);
        y2 = L + (dL1 * 0.5);
        z2 = theta + (dtheta1 * 0.5);
        dL2 = dt * dLdt(k, linf, z2, y2);
        dtheta2 = dt * dthetadt(k, linf, eps, z2, thetae);
        y3 = L + dL2;
        z3 = theta + z2;
        dL3 = dt * dLdt(k, linf, z3, y3);
        dtheta3 = dt * dthetadt(k, linf, eps, z3, thetae);

        dL = (dL0 + 2 * dL1 + 2 * dL2 + dL3) / 6;
        dl = pow(2 * dL, 0.5);
        dtheta = (dtheta0 + 2 * dtheta1 + 2 * dtheta2 + dtheta3) / 6;
        L = L + dL;

        fprintf(fp, "%e,%e,%e,%e,%e\n", t, l, dl, theta, dtheta);
        t = t + dt;
        l = pow(2 * L, 0.5);
        theta = theta + dtheta;
    }

    fclose(fp);
    printf("%i回の計算で収束しました.\n", i);
    printf("平衡上昇水位[m]:%e\n", linf);
    printf("最終上昇水位[m]:%e\n", l);
    printf("%sファイルへの書き込みが完了しました.\n", fname);
    return 0;
}

double dLdt(double k, double linf, double theta, double L){
    double dLdt;
    dLdt = k * (linf * cos(theta) - pow(2 * L, 0.5));
    return dLdt;
}

double dthetadt(double k, double linf, double eps, double theta, double thetae){
    double thetaparams;
    double thetaparams2;
    double dthetadt;
    thetaparams = eps * theta *(pow(theta, 2) - pow(thetae, 2));
    thetaparams2 = eps * (3 * pow(theta, 2) - pow(thetae, 2));
    dthetadt = (-1 /(k * linf)) * pow(k + thetaparams, 2) * thetaparams / (sin(theta) * (k + thetaparams) + thetaparams2);
    return dthetadt;
}
