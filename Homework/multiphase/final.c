#include<stdio.h>
#include<math.h>

#define rhoBubble 1.293
#define rhoWater 1.0e3
#define rBubble 5.0e-8
#define nuWater 1.0e-3
#define g 9.81
#define w0 0.0
#define t0 0.0
#define dt 1.0e-20
#define loop 10000

double w, dw, dw0, dw1, dw2, dw3;
double t;
double A, B, C;
int i;
FILE *fp;
char *fname = "output.csv";

int main(){
    fp = fopen(fname, "w");
    if(fp == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }

    A = (4.0 / 3.0) * M_PI * rBubble * rBubble * rBubble * (rhoBubble + 0.5 * rhoWater);
    B = 6.0 * M_PI * nuWater * rBubble;
    C = (4.0 / 3.0) * M_PI * rBubble * rBubble * rBubble * (rhoWater - rhoBubble) * g;
    printf("%e, %e, %e, %e\n", A, B, C, C/B);

    w = w0;
    t = t0;

    i = 0;
    while( i < loop ){
        fprintf(fp, "%e,%e\n", t, w);
        dw0 = (dt / A) * ( -1 * B * w + C);
        dw1 = (dt / A) * ( -1 * B * (w + 0.5 * dw0) + C);
        dw2 = (dt / A) * ( -1 * B * (w + 0.5 * dw1) + C);
        dw3 = (dt / A) * ( -1 * B * (w + 0.5 * dw2) + C);
        dw = (dw0 + 2.0 * dw1 + 2.0 * dw2 + dw3) / 6.0;
        w += dw;
        t += dt;
        i++;
    }

    fclose(fp);
    return 0;
}
