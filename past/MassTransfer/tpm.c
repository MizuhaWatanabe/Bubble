#include<stdio.h>
#include<math.h>

int main(void){
    double P0 = 1;
    double T = 293;
    double R = 0.00008207;
    double pie = 3.1415926535;
    double lmolmass = 0.032;            /* l means liquid phase */
    double H = 73100;
    double D = 0.0000000026;
    double lsur = 0.071;
    double lvis = 0.001002;
    double ldens = 1000;
    double g = 9.8;
    double e = 0.00000000070832;        /* Permittivity of water */
    double vol = 5;                     /* Surface Voltage of bubble */
    double He0 = 1;                     /* Initial Henry Pressure */

    double r, t, dt, Py, Pe, P, V, A, n, c, d, U, kl, L, VL, cL, PL, N, dr;         /* L means Boundary Film */

    double YoungLaplace(double P0, double lsur, double r);
    double SphereVolume(double pie, double r);
    double SphereSurface(double pie, double r);
    double IdealGasLaw(double P, double V, double R, double T);
    double Concentration(double n, double V);
    double Density(double n, double lmolmass, double V);
    double Stokes(double ldens, double d, double g, double r, double lvis);
    double TransferCoefficient(double D, double r, double U);
    double Trasnfermol(double kl, double A, double P, double PL, double H);
    double Differential(double N, double R, double T, double pie, double r, double P0, double lsur, double e, double vol);
    double FilmThick(double D, double kl);
    double Henry(double H, double cL);
    double ElectroCapillary(double e, double vol, double r);                /* ElectroCapillary Pressure */

    FILE *fp;
    char *fname = "bubble.csv";
    char *l1 = "time[s]";
    char *l2 = "dr/dt[m/s]";
    char *l3 = "radius[m]";
    char *l4 = "laplace[atm]";
    char *l5 = "electrocap.[atm]";
    char *l6 = "pressure[atm]";
    char *l7 = "volume[m3]";
    char *l8 = "surface[m2]";
    char *l9 = "n[mol]";
    char *l10 = "conc.[mol/m3]";
    char *l11 = "dens.[kg/m3]";
    char *l12 = "velocity[m/s]";
    char *l13 = "kl[m/s]";
    char *l14 = "trans.[mol/s]";
    char *l15 = "film thickness[m]";
    char *l16 = "film volume[m3]";
    char *l17 = "film conc.[mol/m3]";
    char *l18 = "Henry[atm]";

    fp = fopen(fname, "w");
    if(fp == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }

    r = 0.001;
    t = 0;
    fprintf(fp, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18);

    int i=0;
    while(r > 1.0e-10 && i < 100000 ){
        i += 1;
        Py = YoungLaplace(P0, lsur, r);
        Pe = ElectroCapillary(e, vol, r);
        P = Py + Pe;
        V = SphereVolume(pie, r);
        A = SphereSurface(pie, r);
        n = IdealGasLaw(P, V, R, T);
        c = Concentration(n, V);
        d = Density(n, lmolmass, V);

        U = Stokes(ldens, d, g, r, lvis);
        kl = TransferCoefficient(D, r, U);
        N = Trasnfermol(kl, A, P, PL, H);
        dr = Differential(N, R, T, pie, r, P0, lsur, e, vol);

        L = FilmThick(D, kl);
        VL = SphereVolume(pie, (r + L)) - V;
        dt = pow(L, 2) / (6 * D);
        cL = N * dt / VL;
        PL = He0 + Henry(H, cL);

        fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", t, dr, r, Py, Pe, P, V, A, n, c, d, U, kl, N, L, VL, cL, PL);
        t = t + dt;
        r = r + dr * dt;
    }

    fclose(fp);
    printf("%i回の計算で収束しました.\n", i);
    printf("最終気泡径[m]:%e\n", r);
    printf("%sファイル書き込みが完了しました.\n", fname);
    return 0;
}


double YoungLaplace(double P0, double lsur, double r){
    double Py;
    Py = P0 + (2 * lsur / (r * 101300));
    return Py;
}

double SphereVolume(double pie, double r){
    double V;
    V = 4 * pie * pow(r, 3) / 3;
    return V;
}

double SphereSurface(double pie, double r){
    double A;
    A = 4 * pie * pow(r, 2);
    return A;
}

double IdealGasLaw(double P, double V, double R, double T){
    double n;
    n = (P * V) / (R * T);
    return n;
}

double Concentration(double n, double V){
    double c;
    c = n / V;
    return c;
}

double Density(double n, double lmolmass, double V){
    double d;
    d = (n * lmolmass) / V;
    return d;
}

double Stokes(double ldens, double d, double g, double r, double lvis){
    double U;
    U = ( 2 * (ldens - d) * g * pow(r, 2)) / (9 * lvis);
    return U;
}

double TransferCoefficient(double D, double r, double U){
    double kl;
    kl = (D / (2 * r)) * (1 + pow(1 + (2 * r * U / D), 0.3333333));
    return kl;
}

double Trasnfermol(double kl, double A, double P, double PL, double H){
    double N;
    N =kl * A * (P - PL) * 101300 / H;
    return N;
}

double Differential(double N, double R, double T, double pie, double r, double P0, double lsur, double e, double vol){
    double dr;
    dr = (-3 * N * R * T) / (4 * pie * (3 * P0 * pow(r, 2) + (4 * lsur * r / 101300) - (e * pow(vol, 2) / (2 * 101300))));
    return dr;
}

double FilmThick(double D, double kl){
    double L;
    L = D / kl;
    return L;
}

double Henry(double H, double cL){
    double PL;
    PL = H * cL / 101300;
    return PL;
}

double ElectroCapillary(double e, double vol, double r){
    double Pe;
    Pe = -1 * e * pow((vol / r), 2) / (2 * 101300);
    return Pe;
}
