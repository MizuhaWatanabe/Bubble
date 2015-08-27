# include<stdio.h>
# include<math.h>

# define LOOP 50
# define N 10
# define X 100
# define Y 50

int main(void){

    const double d = 2.5;
    const double t = 0.1;
    const double delta = d / Y;
    const double dt = t / N;
    const double Conv = 1.0e-6;
    const double uinf = 1;

    const double rho = 1;
    const double mu = 1.0e-2;

    double u[N][Y][X] = {};
    double v[N][Y][X] = {};
    double p[N][Y][X] = {};
    int    barrier[N][Y][X] = {};

    double MaxPress;
    double MaxErr;
    double CurErr;
    double Prev_Press;
    double Re, CL;

    int n, i, j, r, w;
    double time[N] = {};
    int reset;
    double timestump;


    double input[3 * Y * X];
    FILE *finput;
    FILE *f;


    /* Start Initial Conditions

    for(i = 0; i < Y; i++){
        for(j = 0; j < 2; j++){
            u[0][i][j] = uinf;
            v[0][i][j] = 0;
        }
    }

    for(i = 0; i < Y; i++){
        for(j = 2; j < X; j++){
            u[0][i][j] = uinf;
            v[0][i][j] = 0;
        }
    }

    for(i = 1; i < Y - 1; i++){
        for(j = 1; j < X - 1; j++){
            p[0][i][j] = 0;
        }
    }

    End Initial Conditions */


    /* Start Setting Barrier */
    for(n = 0; n < N; n++){
        for(i = 23; i < 28; i++){
            for(j = 10; j < 15; j++){
                barrier[n][i][j] = 1;
    /*            u[n][i][j] = v[n][i][j] = 0;  */
            }
        }
    }


/*
    for(i = 23; i < 28; i++){
        u[0][i][10] = u[0][i][14] = 0;
        v[0][i][10] = v[0][i][14] = 0;
        p[0][i][10] = p[0][i][9];
        p[0][i][14] = p[0][i][15];
    }

    for(j = 11; i < 14; i++){
        u[0][23][j] = u[0][27][j] = 0;
        v[0][23][j] = v[0][27][j] = 0;
        p[0][23][j] = p[0][22][j];
        p[0][27][j] = p[0][28][j];
    }

    for(i = 24; i < 27; i++){
        u[0][i][11] = -1 * u[0][i][9];
        u[0][i][13] = -1 * u[0][i][15];
        v[0][i][11] = -1 * v[0][i][9];
        v[0][i][13] = -1 * v[0][i][15];
    }

    u[0][24][12] = -1 * u[0][22][12];
    u[0][26][12] = -1 * u[0][28][12];
    v[0][24][12] = -1 * v[0][22][12];
    v[0][26][12] = -1 * v[0][28][12];
    End Setting Barrier */


    /* Start Reading Input File */
    finput = fopen("input.csv", "r");
    if(finput == NULL){
        printf("cannot open\n");
        return -1;
    }

    for(i = 0; i < 3 * Y * X; i++){
        fscanf(finput, "%le", &(input[i]));
    }

    fclose(finput);

    for(i = 0; i < Y; i++){
        for(j = 0; j < X; j++){
            u[0][i][j] = input[3 * Y * i + 3 * j];
            v[0][i][j] = input[3 * Y * i + 3 * j + 1];
            p[0][i][j] = input[3 * Y * i + 3 * j + 2];
        }
    }
    /* End Reading Input File */


    reset = 0;
    time[0] = 0;

    for(reset = 0; reset < LOOP; reset++){
        for(n = 0; n < N - 1; n++){

            /* Start Calculating velocity */
            for(i = 2; i < Y - 2; i++){
                for(j = 2; j < X - 2; j++){
                    if(barrier[n][i][j] == 0){
                        u[n + 1][i][j] = u[n][i][j] + (dt / (12 * rho * delta * delta)) * (-1 * rho * delta * (u[n][i][j] * (-1 * u[n][i + 2][j] + 8 * (u[n][i + 1][j] - u[n][i - 1][j]) + u[n][i - 2][j]) + 3 * fabs(u[n][i][j]) * (u[n][i + 2][j] - 4 * u[n][i + 1][j] + 6 * u[n][i][j] - 4 * u[n][i - 1][j] + u[n][i - 2][j]) + v[n][i][j] * (-1 * u[n][i][j + 2] + 8 * (u[n][i][j + 1] - u[n][i][j - 1]) + u[n][i][j - 2]) + 3 * fabs(v[n][i][j]) * (u[n][i][j + 2] - 4 * u[n][i][j + 1] + 6 * u[n][i][j] - 4 * u[n][i][j - 1] + u[n][i][j - 2])) - 6 * delta * (p[n][i + 1][j] - p[n][i - 1][j]) + 12 * mu * (u[n][i + 1][j] + u[n][i - 1][j] + u[n][i][j + 1] + u[n][i][j - 1] - 4 * u[n][i][j]));

                        v[n + 1][i][j] = v[n][i][j] + (dt / (12 * rho * delta * delta)) * (-1 * rho * delta * (u[n][i][j] * (-1 * v[n][i + 2][j] + 8 * (v[n][i + 1][j] - v[n][i - 1][j]) + v[n][i - 2][j]) + 3 * fabs(u[n][i][j]) * (v[n][i + 2][j] - 4 * v[n][i + 1][j] + 6 * v[n][i][j] - 4 * v[n][i - 1][j] + v[n][i - 2][j]) + v[n][i][j] * (-1 * v[n][i][j + 2] + 8 * (v[n][i][j + 1] - v[n][i][j - 1]) + v[n][i][j - 2]) + 3 * fabs(v[n][i][j]) * (v[n][i][j + 2] - 4 * v[n][i][j + 1] + 6 * v[n][i][j] - 4 * v[n][i][j - 1] + v[n][i][j - 2])) - 6 * delta * (p[n][i][j + 1] - p[n][i][j - 1]) + 12 * mu * (v[n][i + 1][j] + v[n][i - 1][j] + v[n][i][j + 1] + v[n][i][j - 1] - 4 * v[n][i][j]));
                    }
                }
            }
            /* End Calculating velocity */

            /* Start velocity Boundary Conditions */
            for(i = 1; i < Y - 1; i++){
                if(barrier[n][i][j] == 0){
                    u[n + 1][i][1] = uinf;
                    u[n + 1][i][X - 2] = u[n + 1][i][X - 3];
                    v[n + 1][i][1] = 0;
                    v[n + 1][i][X - 2] = 0;

                    u[n + 1][i][0] = uinf;
                    u[n + 1][i][X - 1] = -1 * u[n + 1][i][X - 3] + 2 * u[n + 1][i][X - 2];
                    v[n + 1][i][0] = 0;
                    v[n + 1][i][X - 1] = 0;
                }
            }

            for(j = 2; j < X - 2; j++){
                if(barrier[n][i][j] == 0){
                    u[n + 1][1][j] = u[n + 1][2][j];
                    u[n + 1][Y - 2][j] = u[n + 1][Y - 3][j];
                    v[n + 1][1][j] = 0;
                    v[n + 1][Y - 2][j] = 0;

                    u[n + 1][0][j] = -1 * u[n + 1][2][j] + 2 * u[n + 1][1][j];
                    u[n + 1][Y - 1][j] = -1 * u[n + 1][Y - 3][j] + 2 * u[n + 1][Y - 2][j];
                    v[n + 1][0][j] = 0;
                    v[n + 1][Y - 1][j] = 0;
                }
            }
            /* End velocity Boundary Conditions */

            /* Start velocity Barrier Conditions */
            for(i = 23; i < 28; i++){
                u[n + 1][i][10] = u[n + 1][i][14] = 0;
                v[n + 1][i][10] = v[n + 1][i][14] = 0;
            }

            for(j = 11; i < 14; i++){
                u[n + 1][23][j] = u[n + 1][27][j] = 0;
                v[n + 1][23][j] = v[n + 1][27][j] = 0;
            }

            for(i = 24; i < 27; i++){
                u[n + 1][i][11] = -1 * u[n + 1][i][9];
                u[n + 1][i][13] = -1 * u[n + 1][i][15];
                v[n + 1][i][11] = -1 * v[n + 1][i][9];
                v[n + 1][i][13] = -1 * v[n + 1][i][15];
            }

            u[n + 1][24][12] = -1 * u[n + 1][22][12];
            u[n + 1][26][12] = -1 * u[n + 1][28][12];
            v[n + 1][24][12] = -1 * v[n + 1][22][12];
            v[n + 1][26][12] = -1 * v[n + 1][28][12];
            /* End velocity Barrier Conditions */

            /* Start Solving Poisson Equation */
            MaxPress = 1.0e-10;
            do{
                MaxErr = CurErr = 0.0;
                for(i = 2; i < Y - 2; i++){
                    for(j = 2; j < X - 2; j++){
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
            for(i = 1; i < Y - 1; i++){
                if(barrier[n + 1][i][j] == 0){
                    p[n + 1][i][1] = 0;
                    p[n + 1][i][X - 2] = -1 * p[n + 1][i][X - 4] + 2 * p[n + 1][i][X - 3];
                }
            }

            for(j = 1; j < X - 1; j++){
                if(barrier[n + 1][i][j] == 0){
                    p[n + 1][1][j] = -1 * p[n + 1][3][j] + 2 * p[n + 1][2][j];
                    p[n + 1][Y - 2][j] = -1 * p[n + 1][Y - 4][j] + 2 * p[n + 1][Y - 3][j];
                }
            }
            /* End Pressure Boundary Conditions */

            /* Start Pressure Barrier Conditions */
            for(i = 23; i < 28; i++){
                p[n + 1][i][10] = p[n + 1][i][9];
                p[n + 1][i][14] = p[n + 1][i][15];
            }

            for(j = 11; i < 14; i++){
                p[n + 1][23][j] = p[n + 1][22][j];
                p[n + 1][27][j] = p[n + 1][28][j];
            }
            /* End Pressure Barrier Conditions */

            time[n + 1] = time[n] + dt;
        }

        /* Start Resetting */
        for(i = 0; i < Y; i++){
            for(j = 0; j < X; j++){
                u[0][i][j] = u[N - 1][i][j];
                v[0][i][j] = v[N - 1][i][j];
                p[0][i][j] = p[N - 1][i][j];
            }
        }
        time[0] = time[N - 1];
        /* End Resetting */

        printf("%.2e\n", time[N - 1]);

    }

    f = fopen("MAC2.csv", "wt");
    for(i = 1; i < Y - 1; i++){
        for(j = 1; j < X - 1; j++){
            fprintf(f, "%e,%e\n%e,%e\n\n", j * delta, i * delta, j * delta +  0.3 * u[(int)(N - 1)][i][j], i * delta + 0.3 * v[(int)(N - 1)][i][j]);
        }
    }
    fclose(f);


    f = fopen("MAC.csv", "wt");
    for(w = 1; w < 9; w++){
        for(i = 1; i < Y - 1; i++){
            for(j = 1; j < X - 1; j++){
                fprintf(f, "%e,", u[(int) (N - 1) * w / 8][i][j]);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }
    fclose(f);

    f = fopen("input.csv", "wt");
    for(i = 0; i < Y; i++){
        for(j = 0; j < X; j++){
            fprintf(f, "%e\n", u[N - 1][i][j]);
            fprintf(f, "%e\n", v[N - 1][i][j]);
            fprintf(f, "%e\n", p[N - 1][i][j]);
        }
    }
    fclose(f);

    Re = rho * d * uinf / mu;
    CL = uinf * dt /delta;
    timestump = time[N - 1];

    printf("前回から%.2es後まで計算しました。\n", timestump);
    if(CL < 1){
        printf("クーラン条件を満たします\nCL = %.2e\n", CL);
    }else{
        printf("クーラン条件を満たしません\nCL = %.2e\n", CL);
    }
    printf("Re = %.2e\n", Re);


    return 0;

}
