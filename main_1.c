#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265
#define e0 8.85418717e-12
#define q 1.60217662e-19 
#define me 9.10938356e-31 

#define L 5 //количество ионов в слое
#define M 3 //количество слоёв по OZ

#define TIME_STEP 1e-17
#define NUM_TIME_STEPS 1000


int main() {  
    double a, z, x0, y0, z0, vx0, vy0, vz0;
  
    printf("Введите расстояние между ионами: ");  
    scanf("%lf", &a);
    printf("Введите заряд иона: ");
    scanf("%lf", &z);
    
    printf("Введите нач.скорость электрона Ox: ");
    scanf("%lf", &vx0);
    printf("Введите нач.скорость электрона Oy: ");
    scanf("%lf", &vy0);
    printf("Введите нач.скорость электрона Oz: ");
    scanf("%lf", &vz0);
    
    printf("Введите нач. координату электрона X: ");
    scanf("%lf", &x0);
    printf("Введите нач. координату электрона Y: ");
    scanf("%lf", &y0);
    printf("Введите нач. координату электрона Z: ");
    scanf("%lf", &z0);
/* 
    a=2e-9;
    z=4*1.60217662e-19;
    
    vx0=1e6;
    vy0=1e6;
    vz0=1e6;

    x0=-3e-9;
    y0=-3e-9;
    z0=0.0;
*/

    z = z*q;
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
    
     
    
    FILE *fiy = fopen("ions.txt", "w"); 
    for (int j=0; j<M; j++){                                       
            for (int i=0; i<L; i++){
                fprintf(fiy, "%e\t%e\t%e\t\n",  x_ion[i][j], y_ion[i][j], z_ion[i][j]);
            }
    }
    
   fclose(fiy);
    
    FILE *fx = fopen("x.txt", "w"); 
    
    double ax, ay, az, r, xn, yn, zn, vxn, vyn, vzn;
    double K=(q*z/(4*PI*e0*me));
    double ax2, ay2, az2, ax3, ay3, az3;


    for (double t = 0; t < NUM_TIME_STEPS*TIME_STEP; t+=TIME_STEP) {
        ax = 0;
        ay = 0;
        az = 0;
        for (int i=0; i<M; i++){
            for (int j=0; j<L; j++){
            
                r = sqrt(pow(x_ion[i][j]-x0,2)+pow(y_ion[i][j]-y0,2)+pow(z_ion[i][j]-z0,2));
           
                ax+=K*(x_ion[i][j]-x0)/pow(r,3);
                ay+=K*(y_ion[i][j]-y0)/pow(r,3);
                az+=K*(z_ion[i][j]-z0)/pow(r,3);
            }
        }
        
        double vxn2 = vx0 + ax * TIME_STEP/2;
        double vyn2 = vy0 + ay * TIME_STEP/2;
        double vzn2 = vz0 + az * TIME_STEP/2;

        double x2 = x0 + vxn2 * TIME_STEP/2 + 0.5 * ax * pow(TIME_STEP/2, 2);
        double y2 = y0 + vyn2 * TIME_STEP/2 + 0.5 * ay * pow(TIME_STEP/2, 2);
        double z2 = z0 + vzn2 * TIME_STEP/2 + 0.5 * az * pow(TIME_STEP/2, 2);

        ax2 = 0;
        ay2 = 0;
        az2 = 0;
        
        for (int i=0; i<M; i++){
            for (int j=0; j<L; j++){
            
                r = sqrt(pow(x_ion[i][j]-x2,2)+pow(y_ion[i][j]-y2,2)+pow(z_ion[i][j]-z2,2));
           
                ax2+=K*(x_ion[i][j]-x2)/pow(r,3);
                ay2+=K*(y_ion[i][j]-y2)/pow(r,3);
                az2+=K*(z_ion[i][j]-z2)/pow(r,3);
            }
        }
        double vxn3 = vx0 - ax * TIME_STEP+2*ax2*TIME_STEP;
        double vyn3 = vy0 - ay * TIME_STEP+2*ay2*TIME_STEP;
        double vzn3 = vz0 - az * TIME_STEP+2*ay2*TIME_STEP;

        double xn = x2 + vxn3 * TIME_STEP + 0.5 * ax2 * pow(TIME_STEP, 2);
        double yn = y2 + vyn3 * TIME_STEP + 0.5 * ay2 * pow(TIME_STEP, 2);
        double zn = z2 + vzn3 * TIME_STEP + 0.5 * az2 * pow(TIME_STEP, 2);


        ax3 = 0;
        ay3 = 0;
        az3 = 0;

        for (int i=0; i<M; i++){
            for (int j=0; j<L; j++){
            
                r = sqrt(pow(x_ion[i][j]-xn,2)+pow(y_ion[i][j]-yn,2)+pow(z_ion[i][j]-zn,2));
           
                ax3+=K*(x_ion[i][j]-xn)/pow(r,3);
                ay3+=K*(y_ion[i][j]-yn)/pow(r,3);
                az3+=K*(z_ion[i][j]-zn)/pow(r,3);
            }
        }

        double axkonch=(ax+4*ax2+ax3)/6;
        double aykonch=(ay+4*ay2+ay3)/6;
        double azkonch=(az+4*az2+az3)/6;

        double vxn = vx0+TIME_STEP*axkonch;
        double vyn = vy0 +TIME_STEP*aykonch;
        double vzn = vz0 +TIME_STEP*azkonch;
        
        if(fabs(xn) >2*a) break;
        if(fabs(yn) >2*a) break;
        if(fabs(zn) >3*a) break;
        
        fprintf(fx, "%e\t%e\t%e\t\n", x0, y0, z0); 

        x0 = xn;
        y0 = yn;
        z0 = zn;
        
        vx0 = vxn;
        vy0 = vyn;
        vz0 = vzn;
    }

    fclose(fx);

    return 0;
}





