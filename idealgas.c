#include<stdio.h>
#include<math.h>

int main(){
    double P0 = 1;
    double T = 293;
    double R = 0.00008207;
    double pie = 3.1415926535;
    double lmolmass = 0.032;
    double H = 73100;
    double D = 0.0000000026;
    double lsur = 0.071;
    double lvis = 0.001002;
    double ldens = 1000;
    double g = 9.8;

    double r = 0.00001;
    double t = 0;
    double dt = 0.001;
    double pressure = P0 + (2 * lsur / (r * 101300));
    double volume = 4 * pie * pow(r, 3) / 3;
    double surface = 4 * pie * pow(r, 2);
    double n = (pressure * volume) / (R * T);
    double conc = n / volume;
    double dens = n * lmolmass / volume;

    double stokeseq = ( 2 * (ldens - dens) * g * pow(r, 2)) / (9 * lvis);
    double kl = (D / (2 * r)) * (1 + pow(1 + (2 * r * stokeseq / D), 0.3333333));
    double transfereq = kl * surface * (pressure - P0) * 101300 / H;
    double dr = (-3 * transfereq * R * T) / (4 * pie * r * (3 * P0 * r + (4 * lsur / 101300)));

    printf("time[s]        dr/dt[m/s]      radius[m]      pressure[atm]  volume[m3]     surface[m2]   n[mol]         conc.[mol/m3]  dens.[kg/m3]   velocity[m/s]  kl[m/s]        transfer[mol/s]\n");
    int i;
    for(i=0;i<3000;i++){
        printf("%e   %e   %e   %e   %e   %e  %e   %e   %e   %e   %e   %e\n", t, dr, r, pressure, volume, surface, n, conc, dens, stokeseq, kl, transfereq);
        t = t + dt;
        r = r + dr * dt;
        pressure = P0 + (2 * lsur / (r * 101300));
        volume = 4 * pie * pow(r, 3) / 3;
        surface = 4 * pie * pow(r, 2);
        n = (pressure * volume) / (R * T);
        conc = n / volume;
        dens = n * lmolmass / volume;

        stokeseq = ( 2 * (ldens - dens) * g * pow(r, 2)) / (9 * lvis);
        kl = (D / (2 * r)) * (1 + pow(1 + (2 * r * stokeseq / D), 0.3333333));
        transfereq = kl * surface * (pressure - P0) * 101300 / H;
        dr = (-3 * transfereq * R * T) / (4 * pie * r * (3 * P0 * r + (4 * lsur / 101300)));
    }

}
