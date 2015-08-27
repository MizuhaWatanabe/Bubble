# include<stdio.h>
# include<math.h>


# define N 100

int main(void){

    const double X = 1.0;
    const int center = (int)(N / 2);
    const double delta = X / N;
    const double Conv = 1.0e-6;

    double phi[N][N];
    double MaxPhi;
    double MaxErr;
    double CurErr;
    double Prev_phi;
    double vx, vy;
    int i, j;
    int loop;
    FILE *f;


    /* Start Boundary Conditions */
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            phi[i][j] = 0;
        }
    }

    phi[(int)center / 2][(int)center / 2] = 1.0;
    phi[(int)3 * center / 2][(int)3 * center / 2] = 1.0;
    for(i = 0; i < N; i++){
        phi[i][i] = 0.1;
    }

    /* End Boundary Conditions */


    loop = 0;
    MaxPhi = 1.0e-10;

    do{
        if(!(loop%1000)) printf("%05d %e\n", loop, MaxPhi);
        MaxErr = CurErr = 0.0;
        for(i = 1; i < N -1; i++){
            for(j = 1; j < N -1; j++){
                Prev_phi = phi[i][j];
                phi[i][j] = 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                if(MaxPhi < fabs(phi[i][j])) MaxPhi = phi[i][j];
                CurErr = (fabs(phi[i][j] - Prev_phi)) / MaxPhi;
                if(MaxErr < CurErr) MaxErr = CurErr;
            }
        }
        loop++;
    }while(MaxErr > Conv);

    f = fopen("Phi.csv", "wt");
    fprintf(f, " ,");
    for(i = 0; i < N; i++){
            fprintf(f, "%e,",  delta * i);
        }
    fprintf(f, "\n");
    for(j = 0; j < N; j++){
        fprintf(f, "%e,", delta * j);
        for(i = 0; i < N; i++){
            fprintf(f, "%e,", phi[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);

    return 0;
}
