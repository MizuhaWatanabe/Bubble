#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define mass       1.0
#define gamma      1.0
#define c_rand     1.0
#define dt         0.1
#define total_step 100000

double fr();

int main(){
    double x, y;
    double u, v;
    double f_mu, f_rand;
    double ran;
    int step;

    FILE *fout;
    char *fname = "BrownMotion.csv";
    char *title1 = "x";
    char *title2 = "y";
    char *title3 = "u";
    char *title4 = "v";
    char *title5 = "f_mu";
    char *title6 = "F_rand";

    fout = fopen(fname, "w");
    if(fout == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }


    x = 0.0;
    y = 0.0;
    u = 0.0;
    v = 0.0;

    for(step=0; step < total_step; step++){
        fprintf(fout, "%e,%e,%e,%e,%e,%e\n", x, y, u, v, f_mu, f_rand);

        f_mu = - 1* gamma * u;
        f_rand = fr();
        x += dt * u;
        y += dt * v;
        u += dt * (-1 * gamma * u + fr()) / mass;
        v += dt * (-1 * gamma * v + fr()) / mass;
        ran = rand() / (double)RAND_MAX * 2;
        printf("%e\n", ran);
    }

    fclose(fout);

    return 0;
}

double fr(){
    double factor;
    factor = sqrt((3 * c_rand) / (2 * dt));
    return (rand() / (double)RAND_MAX * 2 - 1.0) * factor;
}
