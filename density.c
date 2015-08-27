#include <stdio.h>
#include <math.h>

#define Max_particle 20000

#define particle_dis 0.01
#define Radius_for_num_density (2.1*particle_dis)
#define Ghost -1
#define Fluid 0
#define Wall 2
#define Dummy_Wall 3

int    NumberOFParticles;
double Time;
int    Type[Max_particle];
double X[Max_particle],
       Y[Max_particle],
       Z[Max_particle];
double Vx[Max_particle],
       Vy[Max_particle],
       Vz[Max_particle];
double Ax[Max_particle],
       Ay[Max_particle],
       Az[Max_particle];
double Pressure[Max_particle];
double NUM_DEN[Max_particle];
double dt;

/*weight is original letters. dis and re are concluded*/
double weight(double, double);
void   calculateNUM_DEN(void);

/*start line of the main*/
int main(int argc, char **argv)
{
    int  iParticle, iTimestep, nTimestep;
    FILE *fp;
    char filename[2048];
    int  iFile;

    nTimestep=100;
    dt=0.001;
    iFile=0;

    /*read*/
    fp=fopen("input.grid", "r");
    fscanf(fp, "%lf", &Time);
    fscanf(fp, "%d", &NumberOFParticles);
        for(iParticle=0; iParticle<NumberOFParticles; iParticle++){
            fscanf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
            &Type[iParticle],
            &X[iParticle],&Y[iParticle],&Z[iParticle],
            &Vx[iParticle],&Vy[iParticle],&Vz[iParticle],
            &Pressure[iParticle],
            &NUM_DEN[iParticle]
            );
        }
        fclose(fp);

        /*culculate the data at the first time*/
    calculateNUM_DEN();//refer the end culculation fomula

        for(iTimestep=0; iTimestep<=nTimestep; iTimestep++){
        /*write*/
            if((iTimestep%50) == 0){
                sprintf(filename, "output_%04d.prof", iFile);
                fp=fopen(filename,"w");
                fprintf(fp, "%lf\n", Time);
                fprintf(fp, "%d\n", NumberOFParticles);
                for(iParticle=0; iParticle<NumberOFParticles; iParticle++){
                    fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                            Type[iParticle],
                            X[iParticle],Y[iParticle],Z[iParticle],
                            Vx[iParticle],Vy[iParticle],Vz[iParticle],
                            Pressure[iParticle],
                            NUM_DEN[iParticle]
                            );
                }
                fclose(fp);
                iFile++;
            }

    /*calculate and substitute resule at each time*/
    /*initialize accelerarion*/
            for(iParticle=0; iParticle<NumberOFParticles; iParticle++){
                Ax[iParticle]=0.0;
                Ay[iParticle]=0.0;
                Az[iParticle]=0.0;
            }
    /*gravity*/
            for(iParticle=0; iParticle<NumberOFParticles; iParticle++){
                if(Type[iParticle]==Fluid){
                Ax[iParticle]=0.0;
                Ay[iParticle]=-9.8;
                Az[iParticle]=0.0;
                }
            }
    /*update velocity*/
            for(iParticle=0; iParticle<NumberOFParticles; iParticle++){
                Vx[iParticle]+=Ax[iParticle]*dt;
                Vy[iParticle]+=Ay[iParticle]*dt;
                Vz[iParticle]+=Az[iParticle]*dt;
            }
    /*update positon*/
            for(iParticle=0; iParticle<NumberOFParticles; iParticle++){
                X[iParticle]+=Vx[iParticle]*dt;
                Y[iParticle]+=Vy[iParticle]*dt;
                Z[iParticle]+=Vz[iParticle]*dt;
            }
            calculateNUM_DEN();//refer the end culculation fomura
            Time+=dt;
        }
 	return 0;
}

/*maim is finished*/
/*calculate about NUM_DEN*/
double weight( double dis, double re){
    double weightIJ;

    if(dis>=re){
        weightIJ=0.0;
    }else{
        weightIJ=(re/dis)-1.0;
    }
    return weightIJ;
    }

void calculateNUM_DEN(void){
    int i,j;
    double xij, yij,zij;
    double dis, dis2;
    double w;

    for(i=0; i<NumberOFParticles; i++){
        NUM_DEN[i]=0.0;
        for(j=0; j<NumberOFParticles; j++){
            if(j==i) continue;//if j=1, skip the bellow process.
            xij=X[j]-X[i];
            yij=Y[j]-Y[i];
            zij=Z[j]-Z[i];
            dis2=pow(xij,2)+pow(yij,2)+pow(zij,2);
            dis=sqrt(dis2);//distance between base and another
            w=weight(dis2, Radius_for_num_density);
            NUM_DEN[i]+=w;
        }
    }
}
