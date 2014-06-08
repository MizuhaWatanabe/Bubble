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
    double e = 0.0000000000088542;
    double vol = 0;                     /* Surface Voltage of bubble */

    double r;
    double t;
    double dt;
    double P;
    double V;
    double A;
    double n;
    double c;
    double d;
    double U;
    double kl;
    double L;       /* L means Boundary Film */
    double VL;
    double cL;
    double PL;
    double N;
    double dr;

    double YoungLaplace(double P0, double lsur, double r);
    double SphereVolume(double pie, double r);
    double SphereSurface(double pie, double r);
    double IdealGasLaw(double P, double V, double R, double T);
    double Concentration(double n, double V);
    double Density(double n, double lmolmass, double V);
    double Stokes(double ldens, double d, double g, double r, double lvis);
    double TransferCoefficient(double D, double r, double U);
    double Trasnfermol(double kl, double A, double P, double Pl, double H);
    double Differential(double N, double R, double T, double pie, double r, double P0, double lsur);
    double FilmThick(double D, double kl);
    double Henry(double H, double cL);
    double ElectroCapillary(double e, double vol, double r);                                                /* ElectroCapillary Pressure */



    printf("time[s]        dr/dt[m/s]      radius[m]      pressure[atm]  volume[m3]     surface[m2]    n[mol]         conc.[mol/m3]  dens.[kg/m3]   velocity[m/s]  kl[m/s]        trans.[mol/s]  Henry[atm]\n");


    r = 0.00001;
    t = 0;

    int i=0;
    while(i < 10000){
        i += 1;
        P = YoungLaplace(P0, lsur, r) + ElectroCapillary(e, vol, r);
        V = SphereVolume(pie, r);
        A = SphereSurface(pie, r);
        n = IdealGasLaw(P, V, R, T);
        c = Concentration(n, V);
        d = Density(n, lmolmass, V);

        U = Stokes(ldens, d, g, r, lvis);
        kl = TransferCoefficient(D, r, U);
        N = Trasnfermol(kl, A, P, PL, H);
        dr = Differential(N, R, T, pie, r, P0, lsur);

        L = FilmThick(D, kl);
        VL = SphereVolume(pie, (r + L)) - V;
        dt = pow(L, 2) / (6 * D);
        cL = N * dt / VL;
        PL = 0.2 + Henry(H, cL);


        printf("%e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e\n", t, dr, r, P, V, A, n, c, d, U, kl, N, PL);
        t = t + dt;
        r = r + dr * dt;
    }

    printf("\ntime[s]        dr/dt[m/s]      radius[m]      pressure[atm]  volume[m3]     surface[m2]    n[mol]         conc.[mol/m3]  dens.[kg/m3]   velocity[m/s]  kl[m/s]        trans.[mol/s]  Henry[atm]\n");

    printf("\n%iTimes Calculation Complete.\n", i);
}

double YoungLaplace(double P0, double lsur, double r){
    double P;
    P = P0 + (2 * lsur / (r * 101300));
    return P;
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

double Trasnfermol(double kl, double A, double P, double Pl, double H){
    double N;
    N =kl * A * (P - Pl) * 101300 / H;
    return N;
}

double Differential(double N, double R, double T, double pie, double r, double P0, double lsur){
    double dr;
    dr = (-3 * N * R * T) / (4 * pie * r * (3 * P0 * r + (4 * lsur / 101300)));
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
