/* MAQUINA DE INDUCCION ALIMENTADA A LINEA CON CONTROLADOR DIFUSO PI SID_TS

Usando como salida de datos archivos .dat
Programado por Orlando Trejo - USB
Para la Materia: Electrónica de Potencia 2 - Prof. José Restrepo

Parametros de la Máquina:
Rs=0.435 Ls=0.0713 Lr=0.0713 Lm0=0.0693 Lls=0.002 Llr=0.002 Rr=0.816
Voltaje línea-lí­nea 220 Vrms

Librerías:
<LibConstantes.h>

Repositorio:
https://otrejot@bitbucket.org/otrejot/maquina-ind-con-control-difuso

*/

/* //DSP LIBRERIAS

#include<21369.h>
#include<def21369.h>

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include <LibConstantes.h> //EL HEADER DEBE ESTAR EN LA CARPETA DE PROYECTO
//si no esta la libreria no importa, aca estan las cttes para portabilidad
#define M_SQRT23        0.81649658092772603446
#define M_EE             2.71828182845904509079
#define M_LOG2E         1.4426950408889634074
#define M_LOG10E        0.43429448190325182765
#define M_LN2           0.69314718055994530942
#define M_LN10          2.30258509299404568402
#define M_2PI                   6.2831853071795864769
#define M_PI            3.14159265358979323846
#define M_PI_2          1.57079632679489661923
#define M_PI_4          0.78539816339744830962
#define M_1_PI          0.31830988618379067154
#define M_2_PI          0.63661977236758134308
#define M_2PI_3                 2.0943951023932
#define M_2_SQRTPI      1.12837916709551257390
#define M_SQRT2         1.41421356237309504880
#define M_SQRT2_2       0.70710678118654752440
#define M_1_SQRT6               0.40824829046386301637
#define PI              M_PI
#define PI2             M_PI_2
#define M_SQRT3         1.73205080756887719317
#define M_SQRT3_2       0.86602540378443859658
#define M_SQRT6_2       1.22474487139158894066
#define CTE3                    0.816497
#define CTE1                    0.408248
#define CTE2                    0.707107


// Parametros para adquisicion

#define LEN                     5000
#define STRIDE                  1       //Para el dsp

//Frecuencia de muestreo

#define Fs                      10000
#define Fs_1                    1e-4
#define Ts_i                    Fs_1

#define PWMMAX                  1    
#define Vdc                     400

/* // DSP MAPEO DE LAS VARIABLES EN MEMORIA

section ("seg_sdram") float datos1[LEN];
section ("seg_sdram") float datos2[LEN];
section ("seg_sdram") float datos3[LEN];
section ("seg_sdram") float datos4[LEN];
section ("seg_sdram") float datos5[LEN];
section ("seg_sdram") float datos6[LEN];
section ("seg_sdram") float datos7[LEN];
section ("seg_sdram") float datos8[LEN];
section ("seg_sdram") float datos9[LEN];
section ("seg_sdram") float datos10[LEN];

 */

// FUNCIONES PARA CONTROLADOR DIFUSO

float PI_MF(float x, float a, float c);
float Z_MF(float x, float l, float r);
float SID_TS_I(float error, float varerror);
float SID_TS_W(float error, float varerror);

// FUNCIONES DE INTEGRACION

void ode45(int NUM, float *yout,float t0,float tfinal,float *x0,
           void (*objetivo)(float *,float *));      // MIO
void ode  (float *yout,float t0,float tfinal,float *x0,float paso,
           void (*maquina)(float *,float *));         // DEL PROF

// FUNCIONES DEL MODELO DE LA MAQUINA

void maquina_model(float *X, float *Y);

// ECUACION DIFERENCIAL PARA ESTIMACION DE FLUJO

void flujo_est(float *M, float *N);

/* // DSP PLL
void InitPLL_SDRAM(void);
*/

// VARIABLES GLOBALES

void (*objetivo)(float, float) = NULL; //Dos apuntadores: Entrada y Salida
void (*maquina)(float, float) = NULL;  // Para el ODE del Prof
float kr,ks,k1,k2,k3,k4,k5,k6,k7;
float c1,c2,c3,c4,c5,c6,c7;
float Vsx,Vsy;
float Te,Tl=0.0,J=0.089,P=2;           //Parametros de la Maquina


int main(void)
{
// DECLARACIONES DE VARIABLES

    // Tiempo y Ciclos de Trabajo

    float x[]={0,0,0,0,0},y[]={0,0,0,0,0};
    float m[]={0.0,0.0},n[]={0.0,0.0};
    float Da=0,Db=0,Dc=0,Dmax,Dmin,Dx,DT;
    int   mflag;
    float t1,t2,t3,t4,t5,t6,t_a,SLOPE,A,B,C;

    // Parametros de la Maquina

    float Rs=0.435,Ls=0.0713,Lr=0.0713,Lm0=0.0693,Rr=0.816;
    float Tr,Lm;

    // Variables electricas

    float Psi_m=0,Psi_ang=0,Psi_x=0,Psi_y=0;
    float Psi_E;
    float Leq;
    int it;

    // Variables del Modelo de la Maquina

    float theta;
    float angf,dangf,VP;
    float fe,Vll;
    float Isa,Isb,Isc,wr;
    float wr_ref,wr_ant=0.0;
    float Isq,Isd,Isx,Isy;

    // Para el Control

    float xin,theta_e=0;
    float Isa_ref,Isb_ref,Isc_ref;
    float Err_a=0,Err_b=0,Err_c=0,Verr_a=0,Verr_b=0,Verr_c=0;
    float wr_err,wr_buff=0.0,Verr_wr=0.0;
    float Ebuff_a=0,Ebuff_b=0,Ebuff_c;

    // Salida de Datos por Archivos

    FILE *ferr,*fverr;             // Para Control
    FILE *fIsa,*fIsaRef,*fDa,*fWr,*fTe;
    FILE *fPsi,*fPsiE;

    // Otras

    float dummy;

/* // DSP INSTRUCCIONES
    InitPLL_SDRAM();
*/

// INICIALIZACION DE VARIABLES

    Tr=Lr/Rr;
    Lm=Lm0;

    Vsx=0;
    Vsy=0;
    Te=0;
    theta=0.0;

    for(it=0;it<5;it++)   x[it]=y[it]=0;

    Vll=220;    //Alimentacion linea linea RMS
    fe=60;
    wr_ref=100.0;

    // El voltaje linea-linea es 220 Vrms => Vap=220*sqrt(2)/sqrt(3)
    VP=Vll*2.0/(3.0*Vdc); // Vsx=sqrt(2/3)*(va+vb.a+vc.a^2)=sqrt(2/3)*(1.5*va)
    angf=0;
    dangf=M_2PI*Ts_i*fe;//dangf=M_PI*0.012;
    SLOPE=Ts_i/(2.0*PWMMAX);
    wr=0;
    xin=0;

    // Variables de Control

    mflag=0;

// MODELO DE LA MAQUINA (en coordenadas de estator)

    Lm=Lm0;
    ks=1.0/(Ls*Lr-Lm*Lm);
    kr=1.0/Lr;
    k1=Lr*ks;      // Leq=1/k1;
    Leq=1.0/k1;     //<---Leq
    Tr=Lr/Rr;
    k2=(Lr*Lr*Rs+Lm*Lm*Rr)*ks*kr;
    k4=Lm*ks;
    k6=Rr*kr;
    k5=k6*Lm;
    k3=k5*ks;
    k7=Lm/Lr;

// PARA LA ESTIMACION DE FLUJO DE ROTOR

    c1=Rr/Lr;
    c2=c1*Lm;

// ABRIR ARCHIVOS

    ferr=fopen("ErrorW.dat", "w");
    fIsaRef=fopen("Isa_ref.dat", "w");
    fverr=fopen("VErrorW.dat", "w");
    fIsa=fopen("Isa.dat", "w");
    fDa=fopen("Da.dat", "w");
    fWr=fopen("Wr.dat", "w");
    fTe=fopen("Te.dat", "w");
    fPsi=fopen("Psi_r_m.dat", "w");
    fPsiE=fopen("Psi_E.dat", "w");

/* // DIBUJAR LOS CONJUNTOS DE PERTENENCIA DEL ERROR Y SU VARIACION

    int id;
    float Xd,Yne,Yze,Ype;

    float Fs_dif;
    int LEN_dif;

    Fs_dif = 1.0e-2;
    LEN_dif = 20000;

    FILE *fne,*fze,*fpe;
    fne=fopen("CPderrorN.dat", "w");
    fze=fopen("CPderrorZ.dat", "w");
    fpe=fopen("CPderrorP.dat", "w");

    for (id=0;id<LEN_dif;id++)
    {
       Xd = -100.0+Fs_dif*id;

       Yne = Z_MF(Xd,-110.0,-0.1);
       Yze = PI_MF(Xd,55.0,0.0);
       Ype = 1-Z_MF(Xd,0.1,110.0);

       fprintf(fne, "%f\n\r",Yne); //Imprimo los resultados
       fprintf(fze, "%f\n\r",Yze); //Imprimo los resultados
       fprintf(fpe, "%f\n\r",Ype); //Imprimo los resultados

    }

    fclose(fne);
    fclose(fze);
    fclose(fpe);

*/

/* //DEBUGGER PARA PROBAR SID_TS TOMANDO UNA MUESTRA DEL ERROR Y SU VAR

    float Perror,Pvarerror;

    while(1)
    {
        printf("ERROR = ");
        scanf("%f",&Perror);
        printf("\n\r");
        printf("VARERROR = ");
        scanf("%f",&Pvarerror);
        printf("\n\r");

        if ((Perror==123)||(Pvarerror==123)) break; // sale del debugger

        printf("RESULTADO = %f\n\r",SID_TS_W(Perror,Pvarerror));

    }

*/

// LAZO PRINCIPAL DEL PROGRAMA

    for(it=0;it<LEN*STRIDE;it++)
      {

        wr_err = wr_ref-wr;

        if (it!=0)  Verr_wr = wr_err - wr_buff;
        wr_buff = wr_err;

        theta_e=atan2f(x[3],x[2]);

        Isa=M_SQRT23*x[0];
        Isb=M_SQRT2_2*x[1]-M_1_SQRT6*x[0];
        Isc=-(Isa+Isb);

        //if (it<5000) Isq=7.0;
        //else Isq=0;

        //Fijo el Flujo
        wr_ant = SID_TS_W(wr_err,Verr_wr);
        Isq = wr_ant;


         Isd = 9.0;

        if((Psi_m<0.6)&&(!mflag))
        {
           Isq=0.0;
           Isd=20.0;
        }
        else    mflag=1;

        m[1]=Isd;

        Isx=Isd*cos(theta_e)-Isq*sin(theta_e);
        Isy=Isd*sin(theta_e)+Isq*cos(theta_e);

        // Referencias para el control

        Isa_ref=M_SQRT23*Isx;
        Isb_ref=M_SQRT2_2*Isy-M_1_SQRT6*Isx;
        Isc_ref=-(Isa_ref+Isb_ref);

        Err_a = Isa_ref-Isa;
        Err_b = Isb_ref-Isb;
        Err_c = Isc_ref-Isc;

        fprintf(ferr, "%f\n\r",wr_err); //Imprimo los resultados

        if(it!=0)
        {
           Verr_a = Err_a - Ebuff_a;
           Verr_b = Err_b - Ebuff_b;
           Verr_c = Err_c - Ebuff_c;
        }

        Ebuff_a = Err_a;
        Ebuff_b = Err_b;
        Ebuff_c = Err_c;

        fprintf(fverr, "%f\n\r",Verr_wr); //Imprimo los resultados
        fprintf(fIsaRef, "%f\n\r",Isa_ref); //Imprimo los resultados

/*        // Controlador Tipo BangBang

        if(Isa<Isa_ref) Da=PWMMAX;
        else            Da=0;

        if(Isb<Isb_ref) Db=PWMMAX;
        else            Db=0;

        if(Isc<Isc_ref) Dc=PWMMAX;
        else            Dc=0;
*/

        // Integral del Controlador Difuso

        //Da_buff=SID_TS_I(Err_a,Verr_a);
        //Db_buff=SID_TS_I(Err_b,Verr_b);
        //Dc_buff=SID_TS_I(Err_c,Verr_c);

        Da=SID_TS_I(Err_a,Verr_a);
        Db=SID_TS_I(Err_b,Verr_b);
        Dc=SID_TS_I(Err_c,Verr_c);

        fprintf(fDa, "%f\n\r",Da); //Imprimo los resultados

        // Calculo de los tiempos para PWM

        DT=Da+Db+Dc;
        Dmax=fmax(fmax(Da,Db),Dc);
        Dmin=fmin(fmin(Da,Db),Dc);
        Dx=DT-Dmax-Dmin;

        t1=Dmin*SLOPE;
        t2=Dx*SLOPE;
        t3=Dmax*SLOPE;
        t4=Ts_i-t3;
        t5=Ts_i-t2;
        t6=Ts_i-t1;
        t_a=0.0;

        // Aplicacion de los ciclos de trabajo

        if(t1!=0)       //if(Dmin != 0)
        {
           A=B=C=1;
           Vsx=0;     //Vdc*(A-0.5*B-0.5*C);
           Vsy=0;     //Vdc*M_SQRT3_2*(B-C);

           ode(y,t_a,t1,x,t1,maquina_model);

           x[0]=y[0];
           x[1]=y[1];
           x[2]=y[2];
           x[3]=y[3];
           x[4]=y[4];
           t_a=t1;
        }


        if(t1!=t2)     // Analisis de t1->t2
        {
           A=B=C=1;

           if(Da<=Dmin) A=0;
           if(Db<=Dmin) B=0;
           if(Dc<=Dmin) C=0;

           Vsx=Vdc*(A-0.5*B-0.5*C);
           Vsy=Vdc*M_SQRT3_2*(B-C);

           ode(y,t_a,t2,x,(t2-t_a),maquina_model);

           x[0]=y[0];
           x[1]=y[1];
           x[2]=y[2];
           x[3]=y[3];
           x[4]=y[4];
           t_a=t2;
         }


        if(t2!=t3)     // Analisis de t2->t3
        {
           A=B=C=1;

           if(Da<=Dx) A=0;
           if(Db<=Dx) B=0;
           if(Dc<=Dx) C=0;

           Vsx=Vdc*(A-0.5*B-0.5*C);
           Vsy=Vdc*M_SQRT3_2*(B-C);

           ode(y,t_a,t3,x,(t3-t_a),maquina_model);

           x[0]=y[0];
           x[1]=y[1];
           x[2]=y[2];
           x[3]=y[3];
           x[4]=y[4];
           t_a=t3;
        }

        if(t3!=t4)      // Analisis de t3->t4
        {
           A=B=C=0;

           Vsx=0;       //Vdc*(A-0.5*B-0.5*C);
           Vsy=0;       //Vdc*M_SQRT3_2*(B-C);

           ode(y,t_a,t4,x,(t4-t_a),maquina_model);

           x[0]=y[0];
           x[1]=y[1];
           x[2]=y[2];
           x[3]=y[3];
           x[4]=y[4];
           t_a=t4;
        }

        if(t4!=t5)      // Analisis de t4->t5
        {
           A=B=C=1;

           if(Da<Dmax) A=0;
           if(Db<Dmax) B=0;
           if(Dc<Dmax) C=0;

           Vsx=Vdc*(A-0.5*B-0.5*C);
           Vsy=Vdc*M_SQRT3_2*(B-C);

           ode(y,t_a,t5,x,(t5-t_a),maquina_model);

           x[0]=y[0];
           x[1]=y[1];
           x[2]=y[2];
           x[3]=y[3];
           x[4]=y[4];
           t_a=t5;
        }

        if(t5!=t6)      // Analisis de t5->t6
        {
           A=B=C=1;

           if(Da<Dx) A=0;
           if(Db<Dx) B=0;
           if(Dc<Dx) C=0;

           Vsx=Vdc*(A-0.5*B-0.5*C);
           Vsy=Vdc*M_SQRT3_2*(B-C);

           ode(y,t_a,t6,x,(t6-t_a),maquina_model);

           x[0]=y[0];
           x[1]=y[1];
           x[2]=y[2];
           x[3]=y[3];
           x[4]=y[4];

           t_a=t6;
        }

        if(t6!=Ts_i)    // Analisis de t6->Ts_i
        {
           A=B=C=1;

           Vsx=0;     //Vdc*(A-0.5*B-0.5*C);
           Vsy=0;     //Vdc*M_SQRT3_2*(B-C);

           ode(y,t_a,Ts_i,x,Ts_i-t6,maquina_model);

           x[0]=y[0];
           x[1]=y[1];
           x[2]=y[2];
           x[3]=y[3];
           x[4]=y[4];
        }

        Te=P*k7*(y[1]*y[2]-y[0]*y[3]);

        if(it%STRIDE==0)            // Para Sub-Muestrear
        {
           Isa=x[0]*M_SQRT23;
           Isb=-x[0]*M_1_SQRT6+x[1]*M_SQRT2_2;
           wr=x[4];

           Psi_x=x[2];
           Psi_y=x[3];

           ode45(2,n,0,1e-4,m,flujo_est);

           m[0]=n[0];

           Psi_E = m[0];

           Psi_m=sqrtf(Psi_x*Psi_x+Psi_y*Psi_y);
           Psi_ang=atan2(Psi_y,Psi_x);

           /*    // DSP SALIDA DE DATOS

           datos1[it/STRIDE]=Isa;
           datos2[it/STRIDE]=Isa_ref;
           datos3[it/STRIDE]=wr;
           datos4[it/STRIDE]=Te;
           //datos5[it/STRIDE]=Km;
           //datos6[it/STRIDE]=ewr;
           datos5[it/STRIDE]=Psi_x=x[2];
           datos6[it/STRIDE]=Psi_y=x[3];
           datos7[it/STRIDE]=Psi_m=sqrtf(Psi_x*Psi_x+Psi_y*Psi_y);
           //datos9[it/STRIDE]=flim;
           //datos10[it/STRIDE]=Vsy;

           */
        }

        fprintf(fIsa, "%f\n\r",Isa);    //Imprimo los resultados
        fprintf(fWr, "%f\n\r",wr);      //Imprimo los resultados
        fprintf(fTe, "%f\n\r",Te);      //Imprimo los resultados
        fprintf(fPsi, "%f\n\r",Psi_m);      //Imprimo los resultados
        fprintf(fPsiE, "%f\n\r",Psi_E);      //Imprimo los resultados

     }  // FIN DEL LAZO PRINCIPAL

    // Cierro los archivos
    fclose(ferr);
    fclose(fverr);
    fclose(fIsa);
    fclose(fIsaRef);
    fclose(fWr);
    fclose(fTe);
    fclose(fPsi);
    fclose(fPsiE);

    dummy = M_PI+CTE3;

    system ("PAUSE");

 /*   // DSP TRAMPA DE ARENA

    for(;;)
      {
            asm("nop;");
            asm("nop;");
            asm("nop;");
            asm("nop;");
      }
*/

    return 0;
}

//DEFINICIONES DE LAS FUNCIONES (haria una biblioteca de difuso..)

float PI_MF(float x, float a, float c)
{
    //FUNCION DE PERTENENCIA PI
    //ENTRADAS: coordenada x, ancho a y centro c
    //SALIDAS: Valor de la MF en el punto x
    //LIMITANTES: a > 0

    if(a<0)
    {
        printf("Hubo un error en función PI_MF, ancho a<0\n\t");
        return 0.0;

    }else
        {

        if (x <= c-a)                       return 0.0;
        if ((x <= (2*c-a)/2))               return 2*pow((x-c+a)/a,2);
        if ((x <=c) && (x > (2*c-a)/2))     return 1-2*pow((c-x)/a,2);
        if ((x <= (2*c+a)/2) && (x > c))    return 1-2*pow((x-c)/a,2);
        if ((x <= c+a) && (x > (2*c+a)/2))  return 2*pow((c+a-x)/a,2);
        if (x > c+a)                        return 0.0;

        }

    printf("Hubo un error en función PI_MF, no entro en nigun caso !\n\r!");
    return 0.0;

}

float Z_MF(float x, float l, float r)
{
    //FUNCION DE PERTENENCIA SIGMA !NEGATIVA!
    //ENTRADAS: coordenada x, pto_inflex l y pto_inflex r
    //SALIDAS: Valor de la MF en el punto x
    //LIMITANTES: l < r
    //NOTA: Para obtener la positiva se hace 1 - Z_MF

    if(l>=r)
    {
        printf("Hubo un error en función Z_MF, ancho l<r\n\t");
        return 0.0;

    }else
        {

            if (x <= l)                     return 1.0;
            if ((x <= (l+r)/2) && (x > l))  return 1-2*pow((x-l)/(r-l),2);
            if ((x <= r) && (x >(l+r)/2))   return 2*pow((r-x)/(r-l),2);
            if (x > r)                      return 0.0;

        }

    printf("Hubo un error en función Z_MF, no entro en nigun caso!!\n\r");
    return 0.0;

}

float SID_TS_I(float error, float varerror)
{
    //SISTEMA DE INFERENCIA DIFUSO TAKAGI SUGENO
    //ENTRADAS: error, variacion_error
    //SALIDAS: valor crisp
    //LIMITANTES: clip inicial, hay que derivar
    //NOTA: Esto es para ciertas MF, se pueden modificar...

    float eMFn,eMFz,eMFp;
    float vMFn,vMFz,vMFp;

    float R[9]={ 0.25, 0.25, 0.5,      //REGLAS DIFUSAS
                 0.25,  0.5, 0.75,
                  0.5,  0.75, 0.75};

    float w[9]={0.0}; //Para la composición de CP
    float z=0.0,wsum=0.0;

    int ir;

    //FUZZIFICACION

    eMFn = Z_MF(error,-4.0,-0.05);
    eMFz = PI_MF(error,2.0,0.0);
    eMFp = 1-Z_MF(error,0.05,4.0);

    vMFn = Z_MF(varerror,-10.0,-0.05);
    vMFz = PI_MF(varerror,8.0,0.0);
    vMFp = 1-Z_MF(varerror,0.05,10.0);

    //COMPOSICION DE LOS VALORES DIFUSOS (se usa producto)

    w[0]=eMFn*vMFn;     //e N y ve N = NG
    w[1]=eMFn*vMFz;     //e N y ve Z = NP
    w[2]=eMFn*vMFp;     //e N y ve P = Z
    w[3]=eMFz*vMFn;     //e Z y ve N = NP
    w[4]=eMFz*vMFz;     //e Z y ve Z = Z
    w[5]=eMFz*vMFp;     //e Z y ve P = PP
    w[6]=eMFp*vMFn;     //e P y ve N = Z
    w[7]=eMFp*vMFz;     //e P y ve Z = PP
    w[8]=eMFp*vMFp;     //e P y ve P = PG

    //TAKAGI-SUGENO

    for(ir=0;ir<3*3;ir++)   wsum += w[ir];
    for(ir=0;ir<3*3;ir++)   z += w[ir]*R[ir]/wsum;

    return z;

}

float SID_TS_W(float werror, float wvarerror)
{
    //SISTEMA DE INFERENCIA DIFUSO TAKAGI SUGENO para la velocidad
    //ENTRADAS: error, variacion_error
    //SALIDAS: valor crisp
    //LIMITANTES: clip inicial, hay que derivar
    //NOTA: Esto es para ciertas MF, se pueden modificar...

    float eMFn,eMFz,eMFp;
    float vMFn,vMFz,vMFp;

    float R[9]={-50.0, -50.0, -80.0,      //REGLAS DIFUSAS
                25.0, 1.9, -2.0,
                55.0, 35.0, 35.0};

    float w[9]={0.0}; //Para la composición de CP
    float z=0.0,wsum=0.0;

    int ir;

    //NORMALIZACION

    wvarerror *=100.0;

    //FUZZIFICACION

    eMFn = Z_MF(werror,-30.0,-0.01);
    eMFz = PI_MF(werror,30,0.0);
    eMFp = 1-Z_MF(werror,0.01,30.0);

    vMFn = Z_MF(wvarerror,-2.0,-0.01);
    vMFz = PI_MF(wvarerror,1.0,0.0);
    vMFp = 1-Z_MF(wvarerror,0.01,2.0);

    //COMPOSICION DE LOS VALORES DIFUSOS (se usa producto)

    w[0]=eMFn*vMFn;     //e N y ve N = NG
    w[1]=eMFn*vMFz;     //e N y ve Z = NP
    w[2]=eMFn*vMFp;     //e N y ve P = Z
    w[3]=eMFz*vMFn;     //e Z y ve N = NP
    w[4]=eMFz*vMFz;     //e Z y ve Z = Z
    w[5]=eMFz*vMFp;     //e Z y ve P = PP
    w[6]=eMFp*vMFn;     //e P y ve N = Z
    w[7]=eMFp*vMFz;     //e P y ve Z = PP
    w[8]=eMFp*vMFp;     //e P y ve P = PG

    //TAKAGI-SUGENO

    for(ir=0;ir<3*3;ir++)
    {
        wsum += w[ir];
    }

    for(ir=0;ir<3*3;ir++)
    {
        z += w[ir]*R[ir]/wsum;
    }

    return z;

}

void maquina_model(float *X, float *Y)
{
 float Isx,Isy,Flrx,Flry,w;
 //float dIsx,dIsy,dFlrx,dFlry,dw;

 Isx=X[0];
 Isy=X[1];
 Flrx=X[2];
 Flry=X[3];
 w=P*X[4];

 Y[0]=k1*Vsx-k2*Isx+k3*Flrx+k4*w*Flry;
 Y[1]=k1*Vsy-k2*Isy+k3*Flry-k4*w*Flrx;
 Y[2]=k5*Isx-k6*Flrx-w*Flry;
 Y[3]=k5*Isy-k6*Flry+w*Flrx;
 Te=P*k7*(Isy*Flrx-Isx*Flry);
 Y[4]=(Te-Tl)/J;
}

void flujo_est(float *M, float *N)
{
   float Psi_est;

   Psi_est=M[0];

   N[0]=(c2*M[1] - c1*Psi_est);
}

void ode45(int NUM, float *yout,float t0,float tfinal,float *x0,
           void (*objetivo)(float *,float *))
{
    //INTEGRADOR PARA RESOLUCION DE ECUACIONES DIFERENCIALES
    //ENTRADAS: NUM de variables, *resultado, t_inicio, t_fin, *cond.inicial,
    //          paso de int, apuntador a la función objetivo
    //SALIDAS: ninguna, se pasan parámetros por referencia en *yout
    //LIMITANTES: Runge Kutta de 4to orden con error de 5to orden
    //            el *objetivo debe tener parametros consistentes in-out
    //NOTA: x0[0]=Isx, x0[1]=Isy, x0[2]=Flrx; x0[3]=Flry; x0[4]=wr

    int i;
    float t,h;
    float s1[NUM],s2[NUM],s3[NUM],s4[NUM],xt[NUM];

    t = t0;
    h = (tfinal - t);

    objetivo(x0,s1);
    for(i=0;i<NUM;i++) xt[i]=x0[i]+h*s1[i]/2.0;

    objetivo(xt,s2);
    for(i=0;i<NUM;i++) xt[i]=x0[i]+h*s2[i]/2.0;

    objetivo(xt,s3);
    for(i=0;i<NUM;i++) xt[i]=x0[i]+h*s3[i];

    objetivo(xt,s4);
    for(i=0;i<NUM;i++) yout[i]=x0[i]+h*(s1[i]+2.0*s2[i]+2.0*s3[i]+s4[i])/6.0;
}

void ode(float *yout,float t0,float tfinal,float *x0,float paso,
         void (*maquina)(float *,float *))
// x0[0]=Isx, x0[1]=Isy, x0[2]=Flrx; x0[3]=Flry; x0[4]=wr
{
 int i;
 float t,hmax,h;
 float s1[5],s2[5],s3[5],s4[5],xt[5];
 t = t0;
 hmax = (tfinal - t)/paso;
 h = paso;
   // Compute the slopes
   //   s1 = feval(ypfun, t, y); s1 = s1(:);
 maquina(x0,s1);
 for(i=0;i<5;i++)
    xt[i]=x0[i]+h*s1[i]/2.0;

   //  s2 = feval(ypfun, t+h, y+h*s1); s2 = s2(:);
 maquina(xt,s2);
   //  s3 = feval(ypfun, t+h/2, y+h*(s1+s2)/4); s3 = s3(:);
 for(i=0;i<5;i++)
    xt[i]=x0[i]+h*s2[i]/2.0;

 maquina(xt,s3);
 for(i=0;i<5;i++)
    xt[i]=x0[i]+h*s3[i];
 maquina(xt,s4);
 for(i=0;i<5;i++)
   yout[i]=x0[i]+h*(s1[i]+2.0*s2[i]+2.0*s3[i]+s4[i])/6.0;
}
