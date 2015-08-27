#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define maxData 10000
#define fps 1316
#define WindowFlag 3 // 0:rectangular, 1:hann, 2:hamming, 3:blackman

int main()
{
    int k, n, N;
    double f[maxData + 1];

    FILE *source;
    FILE *fourier;
    char *title1 = "DataNumber";
    char *title2 = "RawWave";
    char *title3 = "Window";
    char *title4 = "Input";
    char *title5 = "Real";
    char *title6 = "Imaginary";
    char *title7 = "Frequency";
    char *title8 = "PowerSpectrum";

    source = fopen("source.csv", "r");
    fourier = fopen("fourier.csv", "w");

    for(N = 0; N < maxData; N++){
        if(fscanf(source, "%lf", &f[N]) == EOF){
            break;
        }
    }
    if(N % 2 != 0){
        printf("ERROR:データの個数は偶数にして下さい.\n");
        return -1;
    }

    double Input[N], Freq[N], Window[N], ReF[N], ImF[N], Pow[N];

    fprintf(fourier, "%s,%s,%s,%s,%s,%s,%s,%s\n", title1, title2, title3, title4, title5, title6, title7, title8);
    for(n = 0; n < N; n++){
        ReF[n] = ImF[n] = 0.0;
        for(k = 0; k < N; k++){
            if(WindowFlag == 0){
                Window[k] = 1;
            }
            else if(WindowFlag == 1){
                Window[k] = 0.5 - 0.5 * cos(2 * M_PI * k / (N - 1));
            }
            else if(WindowFlag == 2){
                Window[k] = 0.54 - 0.46 * cos(2 * M_PI * k / (N - 1));
            }
            else if(WindowFlag == 3){
                Window[k] = 0.42 - 0.5 * cos(2 * M_PI * k /(N - 1)) + 0.08 * cos(4 * M_PI * k / (N - 1));
            }
            else{
                printf("ERROR:WindowFlagの値が不適切です.\n");
                return -1;
            }

            ReF[n] += Window[k] * f[k] * ( cos(2 * M_PI * k * n / N));
            ImF[n] += Window[k] * f[k] * (-sin(2 * M_PI * k * n / N));
        }
        Input[n] = Window[n] * f[n];
        Freq[n] = (double)n * fps / N;
        Pow[n] = ReF[n] * ReF[n] + ImF[n] * ImF[n];
        fprintf(fourier, "%d,%f,%f,%f,%f,%f,%f,%f\n", n, f[n], Window[n], Input[n], ReF[n], ImF[n], Freq[n], Pow[n]);
    }

    fclose(source);
    fclose(fourier);

    system("open fourier.csv");
    return 0;
}
