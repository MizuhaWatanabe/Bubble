#include<stdio.h>
#include<math.h>

int main(void){
    double gamma = 0.071;
    double r = 0.00024876;
    double g = 9.8;
    double rho = 1000;
    double eta = 0.001;
    double theta = 9.63;              /* theta has degree dimention */
    double pie = 3.1415926535;
    double alpha = 0.168;           /* alpha = 0.168 according to the refference */
    double beta = 4.66;            /* beta = 4.66 according to the refference */

    double linf = (2 * gamma * cos(theta * pie / 180)) / (r * g * rho);
    double k = (pow(r, 2) * g * rho) / (8 * eta);

    double l, L, dL, dldt, s, t, dt, dL0, dL1, dL2, dL3, y1, y2, y3, F, Feta, Fg;

    double differential(double k, double linf, double alpha, double beta, double t, double L);

    FILE *fp;
    char *fname = "capillary.csv";
    char *l1 = "time[s]";
    char *l2 = "height[m]";
    char *l3 = "velocity[m/s]";
    char *l4 = "capillary_force[N]";
    char *l5 = "viscous_force[N]";
    char *l6 = "gravity_force[N]";

    fp = fopen(fname, "w");
    if(fp == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }

    L = 0;
    l = 0;
    t = 0;
    dt = 0.001;

    fprintf(fp, "%s,%s,%s,%s,%s,%s\n", l1, l2, l3, l4, l5, l6);

    int i = 0;
    while( i < 10000){
        i += 1;
        dL0 = dt * differential(k, linf, alpha, beta, t, L);
        s = t + (dt * 0.5);
        y1 = L + (dL0 * 0.5);
        dL1 = dt * differential(k, linf, alpha, beta, s, y1);
        y2 = L + (dL1 * 0.5);
        dL2 = dt * differential(k, linf, alpha, beta, s, y2);
        s = t + (dt * 0.5);
        y3 = L + dL2;
        dL3 = dt * differential(k, linf, alpha, beta, t, y3);
        dL = (dL0 + 2 * dL1 + 2 * dL2 + dL3) / 6;
        L = L + dL;

        dldt = dL / (l * dt);
        F = 2 * M_PI * r * gamma * cos(theta * M_PI / 180);
        Feta = 8 * M_PI * eta * pow(2 * L, 0.5) * dldt;
        Fg = M_PI * rho * pow(r, 2) * g * pow(2 * L, 0.5);

        fprintf(fp, "%e,%e,%e,%e,%e,%e\n", t, l, dldt, F, Feta, Fg);
        t = t + dt;
        l = pow(2 * L, 0.5);
    }

    fclose(fp);
    printf("%i回の計算で収束しました.\n", i);
    printf("平衡上昇水位[m]:%e\n", linf);
    printf("最終上昇水位[m]:%e\n", l);
    printf("%sファイル書き込みが完了しました.\n", fname);


    return 0;
}

double differential(double k, double linf, double alpha, double beta, double t, double L){
    double dL;
    dL = k * (linf * (1 - alpha * exp(-1 * beta * t)) - pow(2 * L, 0.5));
    return dL;
}
