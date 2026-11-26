#include <stdio.h>
#include <math.h>
#include <string.h>
#define pi 3.1415926
#define G 9.80
enum
{
    X=0,Y,Z,
    Q0,Q1,Q2,Q3,
    X_D,Y_D,Z_D,
    P,Q,R
}; /*13 state variables*/
typedef struct
{
    double m;
    double d;
    double Ix,Iy,Iz;
    double x[17];/*13 state variables and 4 rotor speeds*/
} Rigid6DOF;

typedef struct
{
    double U1,U2,U3,U4;
} Wrench;

typedef struct
{
    double sigma[4];
    double Cr;
    double omegab;
    double Tm;
    double Ct,Cm;
    double omegass[4];
} Actuator;

void rb_init(Rigid6DOF *rb)
{
    memset(rb,0,sizeof(*rb));
    rb->m=1.0;
    rb->Ix=0.0025;
    rb->Iy=0.0025;
    rb->Iz=0.0045;
    rb->d=0.17;
    rb->x[Q0]=1.0;
}


void ac_init(Actuator *ac)
{
    memset(ac,0,sizeof(*ac));
    ac->Cr=1250;
    ac->omegab=0;
    ac->Tm=0.05;
    ac->Ct=5*1e-6;
    ac->Cm=2.5*1e-7;
}

void ac_operate(const Rigid6DOF *rb,Wrench *w,const Actuator *a)
{
    const double *x=rb->x;
    w->U1=a->Ct*(x[13]*x[13]+x[14]*x[14]+x[15]*x[15]+x[16]*x[16]);
    w->U2=sqrt(2)/2*rb->d*a->Ct*(-x[13]*x[13]-x[14]*x[14]+x[15]*x[15]+x[16]*x[16]);
    w->U3=sqrt(2)/2*rb->d*a->Ct*(x[13]*x[13]-x[14]*x[14]-x[15]*x[15]+x[16]*x[16]);
    w->U4=a->Cm*(x[13]*x[13]-x[14]*x[14] +x[15]*x[15]-x[16]*x[16]);
}

void rb_deriv(const Rigid6DOF *rb,Actuator *a,double dx[17])
{
    const double *x=rb->x;
    double omega[4];

    for(int i=0;i<4;i++)
    {
        omega[i]=x[13+i];
        a->omegass[i]=a->Cr*a->sigma[i]+a->omegab;
        dx[13+i]=(a->omegass[i]-omega[i])/a->Tm;
    }
    double q0=x[Q0],q1=x[Q1],q2=x[Q2],q3=x[Q3];
    double x_d=x[X_D],y_d=x[Y_D],z_d=x[Z_D];
    double p=x[P],q=x[Q],r=x[R];
    double mod=sqrt(q0*q0+q1*q1+q2*q2+q3*q3);

    q0=q0/mod;
    q1=q1/mod;
    q2=q2/mod;
    q3=q3/mod;
    dx[0]=x_d;
    dx[1]=y_d;
    dx[2]=z_d;
    Wrench w; /*declare wrench*/
    ac_operate(rb,&w,a);
    dx[X_D]=2*w.U1*(q1*q3+q0*q2)/rb->m;
    dx[Y_D]=2*w.U1*(q2*q3-q0*q1)/rb->m;
    
if(x[Z]>=0&&w.U1*(1-2*q1*q1-2*q2*q2)/rb->m<=G)
dx[Z_D]=0;

else
dx[Z_D]=G-w.U1*(1-2*q1*q1-2*q2*q2)/rb->m;
    dx[P]=(w.U2+(rb->Iy-rb->Iz)*q*r)/rb->Ix;
    dx[Q]=(w.U3+(rb->Iz-rb->Ix)*p*r)/rb->Iy;
    dx[R]=(w.U4+(rb->Ix-rb->Iy)*p*q)/rb->Iz;

    dx[0]=x_d;
    dx[1]=y_d;
    dx[2]=z_d;

    dx[3]=-(q1*p+q2*q+q3*r)/2;
    dx[4]=(q0*p+q2*r-q3*q)/2;
    dx[5]=(q0*q+q3*p-q1*r)/2;
    dx[6]=(q0*r+q1*q-q2*p)/2;
}


/*RK4*/
void rb_step(Rigid6DOF *rb,Actuator *a,double dt)
{
    double x0[17],k1[17],k2[17],k3[17],k4[17],xt[17];
    Rigid6DOF tmp;
    memcpy(x0,rb->x,sizeof(x0));
    rb_deriv(rb,a,k1);
    for(int i=0;i<17;i++)
        xt[i]=x0[i]+0.5*dt*k1[i];
    tmp=*rb;
    memcpy(tmp.x,xt,sizeof(xt));
    rb_deriv(&tmp,a,k2);
    for(int i=0;i<17;i++)
        xt[i]=x0[i]+0.5*dt*k2[i];
    memcpy(tmp.x,xt,sizeof(xt));
    rb_deriv(&tmp,a,k3);
    for(int i=0;i<17;i++)
        xt[i]=x0[i]+dt*k3[i];
    memcpy(tmp.x,xt,sizeof(xt));
    rb_deriv(&tmp,a,k4);
    for(int i=0;i<17;i++)
        rb->x[i]=x0[i]+dt/6.0*
        (k1[i]+2*k2[i]+2*k3[i]+k4[i]);

    /*quaternion normalization*/
    double mod=sqrt(rb->x[Q0]*rb->x[Q0]+rb->x[Q1]*rb->x[Q1]+rb->x[Q2]*rb->x[Q2]+rb->x[Q3]*rb->x[Q3]);
    for(int i=Q0;i<=Q3;i++)
        rb->x[i]/=mod;
}
double phi,theta,psi;
/*quaternion to euler angles*/
void Q2E(Rigid6DOF *rb)
{
    phi=atan2(
        2*(rb->x[3]*rb->x[4]+rb->x[5]*rb->x[6]),
        1-2*(rb->x[4]*rb->x[4]+rb->x[5]*rb->x[5])
    );

    theta=asin(2*(rb->x[3]*rb->x[5]-rb->x[6]*rb->x[4]));

    psi=atan2(
        2*(rb->x[3]*rb->x[6]+rb->x[4]*rb->x[5]),
        1-2*(rb->x[5]*rb->x[5]+rb->x[6]*rb->x[6])
    );
}


int main()
{

    Rigid6DOF rb;
    Actuator ac;

    rb_init(&rb);
    ac_init(&ac);

    const double dt=0.001;
    const double T_end=5.0;
    int N=(int)(T_end/dt);

    /*case 1:ascend,then hover*/
    printf("Hover\n");

    

    printf("  t      x       y       z      phi     theta   psi\n");

    for(int i=0;i<=N;i++)
    {
        double t=i*dt;
        double s;
        if(t<2.85)
        {s=0.58;}
        else if(t<3.0)
        {s=0.58-2/15*(t-2.85);}
        else
        {s=0.56;}
        for(int j=0;j<4;j++)
        {
          ac.sigma[j]=s;
        }
        Q2E(&rb);
        if(i%500==0)
        {
            printf("%5.2f  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
                   t,
                   rb.x[0],rb.x[1],rb.x[2],
                   phi*180/pi,theta*180/pi,psi*180/pi);
        }
        if(i<N)
            rb_step(&rb,&ac,dt);
    }

    /*case 2:constant roll torque*/
    rb_init(&rb);
    ac_init(&ac);

    printf("\nRoll\n");
    ac.sigma[0]=0.5604755;
    ac.sigma[1]=0.5604755;
    ac.sigma[2]=0.56;
    ac.sigma[3]=0.56;
    printf("  t      x       y       z      phi     theta   psi\n");

    for(int i=0;i<=N;i++)
    {
        double t=i*dt;
        Q2E(&rb);

        if(i%500==0)
        {
            printf("%5.2f  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
                   t,
                   rb.x[0],rb.x[1],rb.x[2],
                   phi*180/pi,theta*180/pi,psi*180/pi);
        }

        if(i<N)
            rb_step(&rb,&ac,dt);
    }

    return 0;
}
