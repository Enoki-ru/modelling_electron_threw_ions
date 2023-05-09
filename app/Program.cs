#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 5
#define M 3 

#define PI 3.14159265
#define e0 8.85418717e-12
#define q -1.60217662e-19 //ELECTRON_CHARGE
#define me 9.10938356e-31 //MASS_ELECTRON

#define TIME_STEP 1e-12
#define NUM_TIME_STEPS 1000000


//float NUM_TIME_STEPS=TIME/TIME_STEP;

double coord_step(double (coord_ion)[5][3],double coord_0){
    double accel;
    for (int i =0; i<L; i++){
        for (int j =0; j<M; j++){
            accel+=1/(coord_ion[i][j]-coord_0);
        }
    }
    
    return accel;
}


int main() {  
    double a, z, vx0, vy0, vz0, x0, y0, z0;

    
    a=2e-9;
    z=3*1.60217662e-19;
    vx0=0.0;
    vy0=0.0;
    vz0=2e6;
    x0=1e-10;
    y0=1e-10;
    z0=0.0;
    
    double K=(q*z/(4*PI*e0*me));
    double grid_length = 2*a;

    double x_ion [L][M]; 
    double y_ion [L][M];
    double z_ion [L][M];
    for (int i =0; i<L; i++){          
        for (int j =0; j<M; j++){        
            if (i==0){ 
                x_ion[i][j]=0;
                y_ion[i][j]=0;
                z_ion[i][j]=j*a;
            }
            if (i==1){ 
                x_ion[i][j]=a;
                y_ion[i][j]=a;
                z_ion[i][j]=j*a;
            }
            if (i==2){ 
                x_ion[i][j]=a;
                y_ion[i][j]=-a;
                z_ion[i][j]=j*a;
            }
            if (i==3){ 
                x_ion[i][j]=-a;
                y_ion[i][j]=-a;
                z_ion[i][j]=j*a;
            }
            if (i==4){ 
                x_ion[i][j]=-a;
                y_ion[i][j]=a;
                z_ion[i][j]=j*a;
            }
        }
    }
    
    FILE *fx = fopen("x.txt", "w"); 
    FILE *fy = fopen("y.txt", "w"); 
    FILE *fz = fopen("z.txt", "w"); 
  
    

    for (double t = 0; t < NUM_TIME_STEPS; t+=TIME_STEP) {
        
        double ax = K*coord_step(x_ion,x0);
        double ay = K*coord_step(y_ion,y0);
        double az = K*coord_step(z_ion,z0);

        double xn = x0 + vx0 * TIME_STEP + 0.5 * ax * pow(TIME_STEP, 2);
        double yn = y0 + vy0 * TIME_STEP + 0.5 * ay * pow(TIME_STEP, 2);
        double zn = z0 + vz0 * TIME_STEP + 0.5 * az * pow(TIME_STEP, 2);

        double vxn = vx0 + ax * TIME_STEP;
        double vyn = vy0 + ay * TIME_STEP;
        double vzn = vz0 + az * TIME_STEP;
        
        fprintf(fx, "%e;", x0); 
        fprintf(fy, "%e;", y0);
        fprintf(fz, "%e;", z0);
        
        x0 = xn;
        y0 = yn;
        z0 = zn;
        
        vx0 = vxn;
        vy0 = vyn;
        vz0 = vzn;
        
        if (abs(x0) == a) break;
        if (abs(y0) == a) break;
        if (abs(z0) == grid_length) break;
    }

    fclose(fx);
    fclose(fy);
    fclose(fz);

    return 0;
}


