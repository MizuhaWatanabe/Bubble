#include <stdio.h>
#include <math.h>

int main()
{
    double r;
    double phi;

    for(r = 0.9; r < 4; r += 0.01){
        phi = 4.0*(pow(r, -12)-pow(r, -6));
        printf("%f %f\n", r, phi);
    }
}
