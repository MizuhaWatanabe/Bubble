# include<stdio.h>
# include<math.h>

# define LOOP 15
# define N 100
# define D 50

int main(void){

    const double d = 2.5;
    const double t = 1;
    const double delta = d / D;
    const double dt = t / N;
    const double Conv = 1.0e-6;
    const double uinf = 1;

    const double rho = 40;
    const double mu = 1.0e-1;

    double u[N][D][D] = {};
    double v[N][D][D] = {};
    double p[N][D][D] = {};
    int    barrier[N][D][D] = {};
    double time[N] = {};

    double MaxPress;
    double MaxErr;
    double CurErr;
    double Prev_Press;
    double Re, CL;

    int n, i, j;
    int reset;
    double timestump;
    FILE *f;


    /* Start Initial Conditions */

    for(i = 0; i < D; i++){
        for(j = 0; j < 2; j++){
            u[0][i][j] = uinf;
            v[0][i][j] = 0;
        }
    }

    for(i = 0; i < D; i++){
        for(j = 2; j < D; j++){
            u[0][i][j] = uinf;
            v[0][i][j] = 0;
        }
    }

    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            p[0][i][j] = 0;
        }
    }

    /* End Initial Conditions */

    /* Start Setting Barrier */
    for(n = 0; n < N; n++){
        for(i = 23; i < 28; i++){
            for(j = 10; j < 15; j++){
                barrier[n][i][j] = 1;
                u[n][i][j] = v[n][i][j] = 0;
            }
        }
    }

    for(n = 0; n < N; n++){
        for(i = 23; i < 28; i++){
            u[n][i][10] = u[n][i][14] = 0;
            v[n][i][10] = v[n][i][14] = 0;
            p[n][i][10] = p[n][i][9];
            p[n][i][14] = p[n][i][15];
        }
    }

    for(n = 0; n < N; n++){
        for(j = 11; i < 14; i++){
            u[n][23][j] = u[n][27][j] = 0;
            v[n][23][j] = v[n][27][j] = 0;
            p[n][23][j] = p[n][22][j];
            p[n][27][j] = p[n][28][j];
        }
    }

    for(n = 0;n < N; n++){
        for(i = 24; i < 27; i++){
            u[n][i][11] = -1 * u[n][i][9];
            u[n][i][13] = -1 * u[n][i][15];
            v[n][i][11] = -1 * v[n][i][9];
            v[n][i][13] = -1 * v[n][i][15];
        }
    }

    for(n = 0; n < N; n++){
        u[n][24][12] = -1 * u[n][22][12];
        u[n][26][12] = -1 * u[n][28][12];
        v[n][24][12] = -1 * v[n][22][12];
        v[n][26][12] = -1 * v[n][28][12];
    }


    reset = 0;
    time[0] = 0;

    for(reset = 0; reset < LOOP; reset++){
        for(n = 0; n < N - 1; n++){

            /* Start Calculating velocity */
            for(i = 2; i < D - 2; i++){
                for(j = 2; j < D - 2; j++){
                    if(barrier[n][i][j] == 0){
                        u[n + 1][i][j] = u[n][i][j] + (dt / (12 * rho * delta * delta)) * (-1 * rho * delta * (u[n][i][j] * (-1 * u[n][i + 2][j] + 8 * (u[n][i + 1][j] - u[n][i - 1][j]) + u[n][i - 2][j]) + 3 * fabs(u[n][i][j]) * (u[n][i + 2][j] - 4 * u[n][i + 1][j] + 6 * u[n][i][j] - 4 * u[n][i - 1][j] + u[n][i - 2][j]) + v[n][i][j] * (-1 * u[n][i][j + 2] + 8 * (u[n][i][j + 1] - u[n][i][j - 1]) + u[n][i][j - 2]) + 3 * fabs(v[n][i][j]) * (u[n][i][j + 2] - 4 * u[n][i][j + 1] + 6 * u[n][i][j] - 4 * u[n][i][j - 1] + u[n][i][j - 2])) - 6 * delta * (p[n][i + 1][j] - p[n][i - 1][j]) + 12 * mu * (u[n][i + 1][j] + u[n][i - 1][j] + u[n][i][j + 1] + u[n][i][j - 1] - 4 * u[n][i][j]));

                        v[n + 1][i][j] = v[n][i][j] + (dt / (12 * rho * delta * delta)) * (-1 * rho * delta * (u[n][i][j] * (-1 * v[n][i + 2][j] + 8 * (v[n][i + 1][j] - v[n][i - 1][j]) + v[n][i - 2][j]) + 3 * fabs(u[n][i][j]) * (v[n][i + 2][j] - 4 * v[n][i + 1][j] + 6 * v[n][i][j] - 4 * v[n][i - 1][j] + v[n][i - 2][j]) + v[n][i][j] * (-1 * v[n][i][j + 2] + 8 * (v[n][i][j + 1] - v[n][i][j - 1]) + v[n][i][j - 2]) + 3 * fabs(v[n][i][j]) * (v[n][i][j + 2] - 4 * v[n][i][j + 1] + 6 * v[n][i][j] - 4 * v[n][i][j - 1] + v[n][i][j - 2])) - 6 * delta * (p[n][i][j + 1] - p[n][i][j - 1]) + 12 * mu * (v[n][i + 1][j] + v[n][i - 1][j] + v[n][i][j + 1] + v[n][i][j - 1] - 4 * v[n][i][j]));
                    }
                }
            }
            /* End Calculating velocity */

            /* Start velocity Boundary Conditions */
            for(i = 1; i < D - 1; i++){
                if(barrier[n][i][j] == 0){
                    u[n + 1][i][1] = uinf;
                    u[n + 1][i][D - 2] = u[n + 1][i][D - 3];
                    v[n + 1][i][1] = 0;
                    v[n + 1][i][D - 2] = 0;

                    u[n + 1][i][0] = uinf;
                    u[n + 1][i][D - 1] = -1 * u[n + 1][i][D - 3] + 2 * u[n + 1][i][D - 2];
                    v[n + 1][i][0] = 0;
                    v[n + 1][i][D - 1] = 0;
                }
            }

            for(j = 2; j < D - 2; j++){
                if(barrier[n][i][j] == 0){
                    u[n + 1][1][j] = u[n + 1][2][j];
                    u[n + 1][D - 2][j] = u[n + 1][D - 3][j];
                    v[n + 1][1][j] = 0;
                    v[n + 1][D - 2][j] = 0;

                    u[n + 1][0][j] = -1 * u[n + 1][2][j] + 2 * u[n + 1][1][j];
                    u[n + 1][D - 1][j] = -1 * u[n + 1][D - 3][j] + 2 * u[n + 1][D - 2][j];
                    v[n + 1][0][j] = 0;
                    v[n + 1][D - 1][j] = 0;
                }
            }
            /* End velocity Boundary Conditions */

            /* Start Solving Poisson Equation */
            MaxPress = 1.0e-10;
            do{
                MaxErr = CurErr = 0.0;
                for(i = 2; i < D - 2; i++){
                    for(j = 2; j < D - 2; j++){
                        if(barrier[n + 1][i][j] == 0){
                            Prev_Press = p[n + 1][i][j];
                            p[n + 1][i][j] = 0.25 * (p[n + 1][i + 1][j] + p[n + 1][i - 1][j] + p[n + 1][i][j + 1] + p[n + 1][i][j - 1]) + (0.0625 * rho / dt) * (dt * ((u[n + 1][i + 1][j] - u[n + 1][i - 1][j]) * (u[n + 1][i + 1][j] - u[n + 1][i - 1][j]) + (u[n + 1][i][j + 1] - u[n + 1][i][j - 1]) * (v[n + 1][i + 1][j] - v[n + 1][i - 1][j]) + (v[n + 1][i][j + 1] - v[n + 1][i][j - 1]) * (v[n + 1][i][j + 1] - v[n + 1][i][j - 1])) - 2 * delta * (u[n + 1][i + 1][j] - u[n + 1][i - 1][j] + v[n + 1][i][j + 1] - v[n + 1][i][j - 1]));
                            if(MaxPress < fabs(p[n + 1][i][j])) MaxPress = p[n + 1][i][j];
                            CurErr = (fabs(p[n + 1][i][j] - Prev_Press)) / MaxPress;
                            if(MaxErr < CurErr) MaxErr = CurErr;
                        }
                    }
                }
            }while(MaxErr > Conv);
            /* End Solving Poisson Equation */

            /* Start Pressure Boundary Conditions */
            for(i = 1; i < D - 1; i++){
                if(barrier[n + 1][i][j] == 0){
                    p[n + 1][i][1] = 0;
                    p[n + 1][i][D - 2] = -1 * p[n + 1][i][D - 4] + 2 * p[n + 1][i][D - 3];
                }
            }

            for(j = 1; j < D - 1; j++){
                if(barrier[n + 1][i][j] == 0){
                    p[n + 1][1][j] = -1 * p[n + 1][3][j] + 2 * p[n + 1][2][j];
                    p[n + 1][D - 2][j] = -1 * p[n + 1][D - 4][j] + 2 * p[n + 1][D - 3][j];
                }
            }
            /* End Pressure Boundary Conditions */

            time[n + 1] = time[n] + dt;
        }

        /* Start Resetting */
        for(i = 0; i < D; i++){
            for(j = 0; j < D; j++){
                u[0][i][j] = u[N - 1][i][j];
                v[0][i][j] = v[N - 1][i][j];
                p[0][i][j] = p[N - 1][i][j];
            }
        }
        time[0] = time[N - 1];
        /* End Resetting */

        printf("%e\n", time[N - 1]);

    }

    f = fopen("MAC2.csv", "wt");
    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            fprintf(f, "%e,%e\n%e,%e\n\n", j * delta, i * delta, j * delta +  0.5 * u[(int)(N - 1) / 2][i][j], i * delta + 0.5 * v[(int)(N - 1) / 2][i][j]);
        }
    }
    fclose(f);


    f = fopen("MAC.csv", "wt");
    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            fprintf(f, "%e,", u[(int) (N - 1) / 8][i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            fprintf(f, "%e,", u[(int) 2 * (N - 1) / 8][i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            fprintf(f, "%e,", u[(int) 3 * (N - 1) / 8][i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            fprintf(f, "%e,", u[(int) 4 * (N - 1) / 8][i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            fprintf(f, "%e,", u[(int) 5 * (N - 1) / 8][i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            fprintf(f, "%e,", u[(int) 6 * (N - 1) / 8][i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            fprintf(f, "%e,", u[(int) 7 * (N - 1) / 8][i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
    for(i = 1; i < D - 1; i++){
        for(j = 1; j < D - 1; j++){
            fprintf(f, "%e,", u[N - 1][i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
    fclose(f);

    Re = rho * d * uinf / mu;
    CL = uinf *dt /delta;
    timestump = time[N - 1];

    printf("%es後まで計算しました。\n", timestump);
    if(CL < 1){
        printf("クーラン条件を満たします\nCL = %e\n", CL);
    }else{
        printf("クーラン条件を満たしません\nCL = %e\n", CL);
    }
    printf("Re = %e\n", Re);


    return 0;

}
