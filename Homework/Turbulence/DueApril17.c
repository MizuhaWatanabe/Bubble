#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main(void){
    double const u0 = 0.1;
    double const a = 1;
    int const dt = 1;
    int const T = 1000;
    double const du = 5.0e-2;

    double u;
    int t;
    double range;

    FILE *fp;
    char *fname = "DueApril17.csv";
    char *l1 = "time[s]";
    char *l2 = "u[m/s]";

    fp = fopen(fname, "w");
    if(fp == NULL){
        printf("%sファイルが開けません\n", fname);
        return -1;
    }

    fprintf(fp, "%s,%s\n", l1, l2);

    u = u0;
    t = 0;
    range = u / du;
    range = (int)range;
    fprintf(fp, "%d,%e,%e\n", t, u, range);
    while(t < T){
	u = a * u * (1 - u);
	t = t + dt;
        range = u / du;
	range = (int)range;
        fprintf(fp, "%d,%e,%e\n", t, u, range);
        printf("%d,%e,%e\n", t, u, range);
    }

    fclose(fp);
    printf("%sファイル書き込みが完了しました\n", fname);
    return 0;
}
