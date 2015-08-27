#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N     20000        /* 時間ステップ数 */
#define X     100          /* 空間グリッド数 */
#define time  2.0          /* 時間の計算範囲 */
#define space 1.0          /* 空間の計算範囲 */
#define D     0.1          /* 拡散係数 */
#define write 200          /* 何ステップに1回書き込むか */

double u[N][X] = {};
double du[N][X] = {};
double x[X] = {};
double S[N][X] = {};
double dt;
double dx;
double a[4] = { 0, 0.5, 0.5, 1 };
double b[4] = { (double)1/6, (double)1/3, (double)1/3, (double)1/6};
double varid;
int    i, n, k, m;
FILE   *fp;
char   fname[32];

int main(){

    /* 数値安定性の計算。不安定ならエラー */
    dt = time / N;
    dx = space / X;
    varid = D * dt / (dx * dx);
    if(varid > 0.5){
        printf("ERROR: Calculation is unstable.(%e)\n", varid);
        return -1;
    }
    else{
        printf("Calculation is stable.(%e)\n", varid);
    }

    /* outputフォルダの初期化 */
    system("rm -rf output");
    system("mkdir output");

    /* 初期条件 */
    for(i = 0; i < X; i++){
        u[0][i] = 1;
    }

    /* Runge-Kutta法の計算 */
    for(n = 0; n < N; n++){

        /* 空間の初期化 */
        for(i = 0; i < X; i++){
            x[i] = i * dx;
        }

        /* outputファイルの書き出し */
        if(n % write == 0){
            sprintf(fname, "output/output_%d.csv", m);
            fp = fopen(fname, "w");
            if(fp == NULL){
                printf("cannot open%s\n", fname);
                return -1;
            }
            fprintf(fp, "x,u\n");
            for(i = 0; i < X; i++){
                fprintf(fp, "%e,%e\n", x[i], u[n][i]);
            }
            fclose(fp);
            m += 1;
        }


        /* 第一種境界条件 */
        u[n][X - 1] = exp(-0.1 * (double)n / (double)N);

        for(k = 0; k < 4; k++){
            for(i = 1; i < X - 1; i++){
                u[n][i] += a[k] * du[n][i];
                x[i] += a[k] * dx;
                du[n][i] = dt * D * ((u[n][i + 1] + u[n][i - 1] - 2 * u[n][i])  / (dx * dx) + (u[n][i + 1] - u[n][i - 1]) / (x[i] * dx));
                S[n][i] += b[k] * du[n][i];
            }
        }

        /* 第二種境界条件 */
        u[n][0] = 2 * u[n][1] - u[n][2];

        /* 時間発展 */
        for(i = 0; i < X ; i++){
            u[n + 1][i] = u[n][i] + S[n][i];
        }
    }

    return 0;
}
