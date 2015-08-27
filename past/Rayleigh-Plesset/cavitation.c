#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define rho 1.00e3
#define nu 1.00e-6
#define sigma 7.10e-2
#define r0 10e-6
#define v0 0.00
#define dt 1.00e-10
#define total_step 999999

double t, r, v, p, r_ratio;
int step;

double Pressure(double);
double RayleighPlesset(double, double, double);
void RungeKutta(double *, double *, double);

int main(){
    FILE *fp;
    char *fname = "cavitation.csv";
    char *label1 = "times[s]";
    char *label2 = "r_ratio";
    char *label3 = "pressure[Pa]";
    char *label4 = "velocity[m/s]";
    char *label5 = "radius[m]";

    fp = fopen(fname, "w");
    if(fp == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }
    fprintf(fp, "%s,%s,%s,%s,%s\n", label1, label2, label3, label4, label5);

    t = 0, r = r0, v = v0;

    while(1){
        step += 1;
        p = Pressure(t);
        r_ratio = r / r0;
        RungeKutta(&r, &v, t);
        fprintf(fp, "%e,%e,%e,%e,%e\n", t, r_ratio, p, v, r);
        t = t + dt;
        if(isnan(r) || step > total_step){
            break;
        }
    }

    fclose(fp);
    printf("%i回の計算で収束しました.\n", step);
    printf("%sファイルへの書き込みが完了しました.\n", fname);
    return 0;
}

double Pressure(double t){
    double Pressure;
    Pressure = sin(M_PI * 2 * 50000 * t) * 0.5e6;
    /*if( t > 2e-6 && t < 4e-6){
        Pressure = -30e3;
    }
    else{
        Pressure = 0;
    }*/
    /*Pressure = (rand() / (double)RAND_MAX * 2 - 1.0) * 1e5;*/
    return Pressure;
}

double RayleighPlesset(double r, double v, double t){
    double firstTerm = -3 * v * v / (2 * r);
    double secondTerm = -4 * nu * v / (r * r);
    double thirdTerm = -2 * sigma / (r * r * rho);
    double fourthTerm = Pressure(t) / (r * rho);
    double dvdt = firstTerm + secondTerm + thirdTerm + fourthTerm;
    printf("%e,%e,%e,%e\n", firstTerm, secondTerm, thirdTerm, fourthTerm);

    return dvdt;
}

void RungeKutta(double *r, double *v, double t){
    double dr, dr0, dr1, dr2, dr3, dv, dv0, dv1, dv2, dv3;
    dr0 = dt * *v;
    dv0 = dt * RayleighPlesset(*r, *v, t);
    dr1 = dt * (*v + dv0 * 0.5);
    dv1 = dt * RayleighPlesset(*r + dr0 * 0.5, *v + dv0 * 0.5, t + dt * 0.5);
    dr2 = dt * (*v + dv1 * 0.5);
    dv2 = dt * RayleighPlesset(*r + dr1 * 0.5, *v + dv1 * 0.5, t + dt * 0.5);
    dr3 = dt * (*v + dv2);
    dv3 = dt * RayleighPlesset(*r + dr2, *v + dv2, t + dt);
    dr = (dr0 + dr1 * 2 + dr2 * 2 + dr3) / 6;
    dv = (dv0 + dv1 * 2 + dv2 * 2 + dv3) / 6;
    *r = *r + dr;
    *v = *v + dv;
}
