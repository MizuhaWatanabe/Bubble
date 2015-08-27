#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nx 128
#define ny 48
#define nz 128
#define Re_b 5600

float x[nx], y[ny], z[nz];
float u[nx][ny][nz], v[nx][ny][nz], w[nx][ny][nz];
float p[nx][ny][nz];
float T[nx][ny][nz];
int   i, j, k;
float U[ny], V[ny], W[ny], P[ny];
float u_dash[nx][ny][nz], v_dash[nx][ny][nz], w_dash[nx][ny][nz], p_dash[nx][ny][nz];
float uu[ny], vv[ny], ww[ny], uv[ny];
float u_rms[ny], v_rms[ny], w_rms[ny];
float left_term[ny], right_term[ny], sum[ny];
float u_dash_bar[ny], v_dash_bar[ny], w_dash_bar[ny];
float uy_dash_bar[ny], uz_dash_bar[ny];
float vx_dash_bar[ny], vz_dash_bar[ny];
float wx_dash_bar[ny], wy_dash_bar[ny];
float epsilon[ny], Pk[ny];
float laplaceX[nx - 2][nz - 2];
float laplaceZ[nx - 2][nz - 2];
float Q[nx - 2][nz - 2];
float etak[ny], tauk[ny], uk[ny];

int main(){

    FILE *fp;
    fp = fopen("xyz.txt", "r");
    if(fp == NULL){
        printf("There are no input files called xyz.txt\n");
        return -1;
    }

    for(i = 0; i < nx; i++){
        fscanf(fp, "%e", &x[i]);
    }
    for(j = 0; j < ny; j++){
        fscanf(fp, "%e", &y[j]);
    }
    for(k = 0; k < nz; k++){
        fscanf(fp, "%e", &z[k]);
    }

    fp = fopen("uvwpt.txt", "r");
    if(fp == NULL){
        printf("There are no input files called uvwpt.txt\n");
        return -1;
    }

    for(k = 0; k < nz; k++){
        for(j = 0; j < ny; j++){
            for(i = 0; i < nx; i++){
                fscanf(fp, "%e", &u[i][j][k]);
                fscanf(fp, "%e", &v[i][j][k]);
                fscanf(fp, "%e", &w[i][j][k]);
                fscanf(fp, "%e", &p[i][j][k]);
                fscanf(fp, "%e", &T[i][j][k]);
            }
        }
    }

    char  *fname = "output.csv";
    fp = fopen(fname, "w");

    for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
            for(k = 0; k < nz; k++){
                U[j] += u[i][j][k] / (nx * nz);
                V[j] += v[i][j][k] / (nx * nz);
                W[j] += w[i][j][k] / (nx * nz);
                P[j] += p[i][j][k] / (nx * nz);
                uu[j] += u[i][j][k] * u[i][j][k] / (nx * nz);
                vv[j] += v[i][j][k] * v[i][j][k] / (nx * nz);
                ww[j] += w[i][j][k] * w[i][j][k] / (nx * nz);
                u_rms[j] = pow(uu[j] - U[j] * U[j], 0.5);
                v_rms[j] = pow(vv[j] - V[j] * V[j], 0.5);
                w_rms[j] = pow(ww[j] - W[j] * W[j], 0.5);
            }
        }
    }
    for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
            for(k = 0; k < nz; k++){
                uv[j] += (u[i][j][k] - U[j]) * (v[i][j][k] - V[j]) / (nx * nz);
                u_dash[i][j][k] = u[i][j][k] - U[j];
                v_dash[i][j][k] = v[i][j][k] - V[j];
                w_dash[i][j][k] = w[i][j][k] - W[j];
                p_dash[i][j][k] = p[i][j][k] - P[j];
            }
        }
    }

    fprintf(fp, "etak,tauk,uk'\n");
    for(j = 0; j < ny - 1; j ++){
        for(i = 0; i < nx - 1; i++){
            for(k = 0; k < nz - 1; k++){
                u_dash_bar[j] += (u_dash[i + 1][j][k] - u_dash[i][j][k]) * (u_dash[i + 1][j][k] - u_dash[i][j][k]) / ((x[i + 1] - x[i]) * (x[i + 1] - x[i]) * ((nx - 1) * (nz - 1)));
                v_dash_bar[j] += (v_dash[i][j + 1][k] - v_dash[i][j][k]) * (v_dash[i][j + 1][k] - v_dash[i][j][k]) / ((y[j + 1] - y[j]) * (y[j + 1] - y[j]) * ((nx - 1) * (nz - 1)));
                w_dash_bar[j] += (w_dash[i][j][k + 1] - w_dash[i][j][k]) * (w_dash[i][j][k + 1] - w_dash[i][j][k]) / ((z[k + 1] - z[k]) * (z[k + 1] - z[k]) * ((nx - 1) * (nz - 1)));
                uy_dash_bar[j] += (u_dash[i][j + 1][k] - u_dash[i][j][k]) * (u_dash[i][j + 1][k] - u_dash[i][j][k]) / ((y[j + 1] - y[j]) * (y[j + 1] - y[j]) * ((nx - 1) * (nz - 1)));
                uz_dash_bar[j] += (u_dash[i][j][k + 1] - u_dash[i][j][k]) * (u_dash[i][j][k + 1] - u_dash[i][j][k]) / ((z[k + 1] - z[k]) * (z[k + 1] - z[k]) * ((nx - 1) * (nz - 1)));
                vx_dash_bar[j] += (v_dash[i + 1][j][k] - v_dash[i][j][k]) * (v_dash[i + 1][j][k] - v_dash[i][j][k]) / ((x[i + 1] - x[i]) * (x[i + 1] - x[i]) * ((nx - 1) * (nz - 1)));
                vz_dash_bar[j] += (v_dash[i][j][k + 1] - v_dash[i][j][k]) * (v_dash[i][j][k + 1] - v_dash[i][j][k]) / ((z[k + 1] - z[k]) * (z[k + 1] - z[k]) * ((nx - 1) * (nz - 1)));
                wx_dash_bar[j] += (w_dash[i + 1][j][k] - w_dash[i][j][k]) * (w_dash[i + 1][j][k] - w_dash[i][j][k]) / ((x[i + 1] - x[i]) * (x[i + 1] - x[i]) * ((nx - 1) * (nz - 1)));
                wy_dash_bar[j] += (w_dash[i][j + 1][k] - w_dash[i][j][k]) * (w_dash[i][j + 1][k] - w_dash[i][j][k]) / ((y[j + 1] - y[j]) * (y[j + 1] - y[j]) * ((nx - 1) * (nz - 1)));
            }
        }
        Pk[j] = -1 * uv[j] * (U[j + 1] - U[j]) / (y[j + 1] - y[j]);
        epsilon[j] = u_dash_bar[j] + v_dash_bar[j] + w_dash_bar[j] + uy_dash_bar[j] + uz_dash_bar[j] + vx_dash_bar[j] + vz_dash_bar[j] + wx_dash_bar[j] + wy_dash_bar[j];

        etak[j] = pow(epsilon[j], -0.25);
        tauk[j] = pow(epsilon[j], -0.5);
        uk[j] = pow(epsilon[j], 0.25);
        fprintf(fp, "%e,%e,%e,%e\n", y[j], etak[j], tauk[j], uk[j]);
    }

    for(k = 1; k < nz - 1; k++){
        for(i = 1; i < nx - 1; i++){
            laplaceX[i][k] = (p[i + 1][ny][k] + p[i - 1][ny][k] - 2 * p[i][ny][k]) / ((x[i + 1] - x[i]) * (x[i] - x[i - 1]));
            laplaceZ[i][k] = (p[i][ny][k + 1] + p[i][ny][k - 1] - 2 * p[i][ny][k]) / ((z[k + 1] - z[k]) * (z[k] - z[k - 1]));
            Q[i][k] = -1 * (laplaceX[i][k] + laplaceZ[i][k]);
//            fprintf(fp, "%e,%e,%e,%e\n", x[i], z[k], Q[i][k], p[i][ny][k]);
        }
    }

    fclose(fp);
    system("open output.csv");

    return 0;
}
