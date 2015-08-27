#include <stdio.h>
#include <math.h>
#include <time.h>

#define R 1e-3      /* pipe radius */
#define L 1e-2      /* pipe length */

int iData;          /* some data */
int numberOfData;   /* all data */

void leastSquare(double [], double [], double *, double *, double *);

int main(){

    /* read Q and deltaP in input file */
    FILE *fpInput;
    fpInput = fopen("input.csv", "r");
    if(fpInput == NULL){
        printf("There are no input files\n");
        return -1;
    }
    fscanf(fpInput, "%d", &numberOfData);

    double Q[numberOfData];
    double deltaP[numberOfData];

    for(iData = 0; iData < numberOfData; iData++){
        fscanf(fpInput, "%lf,%lf", &Q[iData], &deltaP[iData]);
    }

    /* calculate gamma, tauWall and their logarhyhms */
    double gamma[numberOfData];
    double tauWall[numberOfData];
    double lnGamma[numberOfData];
    double lnTauWall[numberOfData];

    for(iData = 0; iData < numberOfData; iData++){
        gamma[iData] = 4 * Q[iData] / (M_PI * R * R * R);
        tauWall[iData] = R * deltaP[iData] / (2 * L);
        lnGamma[iData] = log(gamma[iData]);
        lnTauWall[iData] = log(tauWall[iData]);
    }

    /* calculate alpha, beta by least square method */
    double alpha;      /* slope of line */
    double beta;       /* intercept of line */
    double rSquared;   /* r squared value of apploximation*/

    leastSquare(lnGamma, lnTauWall, &alpha, &beta, &rSquared);

    /* calculate gammaDotWall and nu0 */
    double gammaDotWall[numberOfData];
    double nu0;

    for(iData = 0; iData < numberOfData; iData++){
        gammaDotWall[iData] = (3 * alpha + 1) * pow(tauWall[iData] / beta, 1 / alpha) / (4 * alpha);
    }

    nu0 = beta * pow(4 * alpha / (3 * alpha + 1), alpha);

    /* write output file */
    FILE *fpOutput;
    char fname[32];
    time_t t;
    t = time(NULL);
    strftime(fname, sizeof(fname), "results/%Y%m%d%H%M%S.csv", localtime(&t));
    char *label1 = "nu0[Pa s]";
    char *label2 = "n[-]";
    char *label3 = "R^2[-]";
    char *label4 = "Volume Rate[m3/s]";
    char *label5 = "delta P[Pa]";
    char *label6 = "lnGamma[/s]";
    char *label7 = "lnTauWall[Pa]";
    char *label8 = "tauWall[Pa]";
    char *label9 = "gammaDotWall[/s]";

    fpOutput = fopen(fname, "w");
    if(fpOutput == NULL){
        printf("cannot open %s\n", fname);
        return -1;
    }

    fprintf(fpOutput, "%s,%e\n", label1, nu0);
    fprintf(fpOutput, "%s,%e\n", label2, alpha);
    fprintf(fpOutput, "%s,%e\n\n", label3, rSquared);

    fprintf(fpOutput, "%s,%s,%s,%s,%s,%s\n", label4, label5, label6, label7, label8, label9);
    for(iData = 0; iData < numberOfData; iData++){
        fprintf(fpOutput, "%e,%e,%e,%e,%e,%e\n",
                Q[iData], deltaP[iData],
                lnGamma[iData], lnTauWall[iData],
                gammaDotWall[iData], tauWall[iData]);
    }

    fclose(fpOutput);

    /* write on the shell */
    printf("%20s%20e\n", label1, nu0);
    printf("%20s%20e\n", label2, alpha);
    printf("%20s%20e\n\n", label3, rSquared);

    printf("%20s%20s%20s%20s%20s%20s\n", label4, label5, label6, label7, label8, label9);
    for(iData = 0; iData < numberOfData; iData++){
        printf("%20e%20e%20e%20e%20e%20e\n",
                Q[iData], deltaP[iData],
                lnGamma[iData], lnTauWall[iData],
                gammaDotWall[iData], tauWall[iData]);
    }

    return 0;
}


void leastSquare(double x[], double y[], double *alpha, double *beta, double *rSquared){

    /* calculate alpha and beta */
    double n, xi, yi, xixi, xiyi;
    n = xi = yi = xixi = xiyi = 0.0;

    for(iData = 0; iData < numberOfData; iData++){
        n    += 1.0;
        xi   += x[iData];
        yi   += y[iData];
        xixi += x[iData] * x[iData];
        xiyi += x[iData] * y[iData];
    }

    *alpha = (n * xiyi - xi * yi) / (n * xixi - xi * xi);
    *beta = (xixi * yi - xiyi * xi) / (n * xixi - xi * xi);

    /* calculate r squared value */
    double xMean, yMean;
    xMean = xi / n;
    yMean = yi / n;

    double dxdy, dxdx, dydy;

    for(iData = 0; iData < numberOfData; iData++){
        dxdy += (x[iData] - xMean) * (y[iData] - yMean);
        dxdx += (x[iData] - xMean) * (x[iData] - xMean);
        dydy += (y[iData] - yMean) * (y[iData] - yMean);
    }

    *rSquared = (dxdy * dxdy) / (dxdx * dydy);
}
