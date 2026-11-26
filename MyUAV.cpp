#include <stdio.h>
#include <math.h>
#include <string.h>
#define pi 3.1415926
#define G 9.80
#define NX 12
/*12个状态量：x,y,z,phi,theta,psi,x_d,y_d,z_d,p,q,r */
/* 无人机参数输入 */ 
typedef struct
{
/* 质量 轴距 转动惯量 12个状态量 */ 
	double m;
    double l;
	double Ix,Iy,Iz;
	double x[NX] ;
} Rigid6DOF; /*最后命名*/ 

/*力的输入*/
typedef struct
{
	double U1;
	double U2,U3,U4;
 } Wrench;
 
 /*初始化函数*/
 void rb_init(Rigid6DOF *rb)  
 {
 	memset(rb,0,sizeof(*rb));
 	rb->m=1.0;
 	rb->Ix=0.0025;
 	rb->Iy=0.0025;
 	rb->Iz=0.0045;
 	rb->l=0.17;

 }
 
 /*定义导数*/ 
 void rb_deriv(const Rigid6DOF *rb,const Wrench *w,double dx[NX])
 {
 	const double *x=rb->x;
 	double phi=x[3],th=x[4],psi=x[5];
 	double x_d=x[6],y_d=x[7],z_d=x[8];
 	double p=x[9],q=x[10],r=x[11];
 	double spsi=sin(psi),cpsi=cos(psi);
 
    
    dx[6]=-(sin(th)*cpsi*cos(phi)+sin(phi)*spsi)*w->U1/rb->m; /*x的二阶导数*/ 
    dx[7]=-(cos(phi)*sin(th)*spsi-sin(phi)*cpsi)*w->U1/rb->m;
    dx[8]=G-w->U1*cos(th)*cos(phi)/rb->m;
    dx[9]=(w->U2+(rb->Iy-rb->Iz)*q*r)/rb->Ix;
    dx[10]=(w->U3+(rb->Iz-rb->Ix)*p*r)/rb->Iy;
    dx[11]=(w->U4+(rb->Ix-rb->Iy)*p*q)/rb->Iz;
    dx[0]=x_d;
    dx[1]=y_d;
    dx[2]=z_d;
    dx[3]=p+sin(phi)*tan(th)*q+cos(phi)*tan(th)*r;
    dx[4]=cos(phi)*q-sin(phi)*r;
    dx[5]=sin(phi)/cos(th)*q+cos(phi)/cos(th)*r;
}

/* RK4 数值积分*/ 
void rb_step( Rigid6DOF *rb,const Wrench *w,double dt)
{
	double x0[NX],k1[NX],k2[NX],k3[NX],k4[NX],xt[NX];
	Rigid6DOF tmp;
	memcpy(x0, rb->x, sizeof(x0)); 

    rb_deriv(rb, w, k1);

    for (int i = 0; i < NX; i++) 
	
	xt[i] = x0[i] + 0.5*dt*k1[i];
    tmp = *rb; 
	memcpy(tmp.x, xt, sizeof(xt));
    rb_deriv(&tmp, w, k2);

    for (int i = 0; i < NX; i++) xt[i] = x0[i] + 0.5*dt*k2[i];
    memcpy(tmp.x, xt, sizeof(xt));
    rb_deriv(&tmp, w, k3);

    for (int i = 0; i < NX; i++) xt[i] = x0[i] + dt*k3[i];
    memcpy(tmp.x, xt, sizeof(xt));
    rb_deriv(&tmp, w, k4);

    for (int i = 0; i < NX; i++)
        rb->x[i] = x0[i] + dt/6.0*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

}

int main(void) 
{
	Rigid6DOF rb;
	rb_init(&rb);
	const double dt=0.001;
	const double T_end=5.0;


	/*场景1悬停*/ 
	printf("悬停\n"); 
	rb_init(&rb);
	rb.x[3]=0.005 ;
    Wrench w_hover = {  rb.m * G,  0.0, 0.0, 0.0 }; 
    printf("  t      x       y       z    phi   theta   psi\n");
    for (double t = 0.0; t <= T_end + 1e-12; t += dt) {
        if (fmod(t, 0.5) < dt)
            printf("%5.2f  %7.4f %7.4f %7.4f  %5.2f %6.2f %6.2f\n",
                   t, rb.x[0], rb.x[1], rb.x[2],
                   rb.x[3], rb.x[4], rb.x[5]);
        rb_step(&rb, &w_hover, dt);
}
    /*场景2恒定 滚转力矩*/
    printf ("恒定滚转力矩\n");
    rb_init(&rb);
    Wrench w_roll = { rb.m * G,  0.01, 0.0, 0.0 };

    printf("  t     phi(rad)   phi_d(rad/s)   X       Z\n");
    for (double t = 0.0; t <= T_end + 1e-12; t += dt) 
	{
        if (fmod(t, 0.5) < dt)
            printf("%5.2f  %8.3f  %9.4f  %7.3f  %7.3f\n",
                   t, rb.x[3], rb.x[9], rb.x[0], rb.x[2]);
        rb_step(&rb, &w_roll, dt);
    }

    return 0;
}













#if 0
/*  限值函数  */
static double clamp(double v,double low,double high) 
{
	return v<low? low : (v>high? high : v);
}

/* PID 控制器*/
typedef struct
{
	double p_ref[3];
	double psi_ref;
	
	double kp_xy,ki_xy,kd_xy;
	double kp_z,ki_z,kd_z;
	
	double kp_att,ki_att,kd_att;
	double kp_yaw,ki_yaw,kd_yaw;
	
    /* 积分状态 */ 
    double int_ex,int_ey,int_ez;
	double int_ephi,int_eth,int_epsi;
	 
	/* 积分限幅 */
	double int_limit_pos;
	double int_limit_att;
 } controller;
 
 void ctrl_init(controller *c)
 {
 	c->p_ref[0]=20.0;
 	c->p_ref[1]=0.0;
 	c->p_ref[2]=-20.0;
 	c->psi_ref=0.0;
    /* 外环先关闭积分，调试好再打开*/ 
    c->kp_xy=1.0;c->ki_xy=0.0;c->kd_xy=2.0;
	c->kp_z=2.0;c->ki_z=0.0;c->kd_z=3.0;
	
	/*内环 wn≈10，zeta≈1 */
	c->kp_att=100.0;c->ki_att=0.0;c->kd_att=20.0;
	c->kp_yaw=50.0;c->ki_yaw=0.0;c->kd_yaw=10.0;

	c->int_ex=c->int_ey=c->int_ez=0.0;
	c->int_ephi=c->int_eth=c->int_epsi=0.0;
	
	c->int_limit_pos=2.0;
	c->int_limit_att=0.5;
}

void ctrl_update(controller *c,const Rigid6DOF *rb,Wrench *w,double dt)
{
	const double *x=rb->x;
	double phi=x[4],th=x[5],psi=x[6];
	
	double ex=c->p_ref[0]-x[1];
	double ey=c->p_ref[1]-x[2];
	double ez=c->p_ref[2]-x[3];
	
	c->int_ex+=ex*dt;
	c->int_ey+=ey*dt;
	c->int_ez+=ez*dt;
	
	/*误差的积分需要限幅*/ 
	c->int_ex=clamp(c->int_ex,-c->int_limit_pos/c->ki_xy-1e-9,c->int_limit_pos/c->ki_xy+1e-9);
	c->int_ey=clamp(c->int_ey,-c->int_limit_pos/c->ki_xy-1e-9,c->int_limit_pos/c->ki_xy+1e-9);
	c->int_ez=clamp(c->int_ez,-c->int_limit_pos/c->ki_z-1e-9,c->int_limit_pos/c->ki_z+1e-9);

    double ax_des=c->kp_xy*ex+c->ki_xy*c->int_ex-c->kd_xy*x[7];
    double ay_des=c->kp_xy*ey+c->ki_xy*c->int_ey-c->kd_xy*x[8];
    double az_des=c->kp_z*ez+c->ki_z*c->int_ez-c->kd_z*x[9];
    
    /*期望的加速度需要限幅*/ 
    const double a_xy=0.5*G;
    const double a_z=0.4*G;
    ax_des=clamp(ax_des,-a_xy,a_xy);
    ay_des=clamp(ay_des,-a_xy,a_xy);
    az_des=clamp(az_des,-a_z,a_z);
    /*预先输入一个U1，仅用于求解期待滚转、俯仰角*/ 
    double U1_est=rb->m*G;
    double phi_des=asin(clamp(rb->m*ay_des/U1_est,-0.99,0.99));
    double theta_des=asin(clamp(-rb->m*ax_des/(U1_est*cos(phi_des)),-0.99,0.99));
    
    /*期待的滚转、俯仰角需要限幅*/ 
    double const TILT=30.0*pi/180.0;
	phi_des = clamp(phi_des,-TILT,TILT) ;
	theta_des=clamp(theta_des,-TILT,TILT);
	/*求出期待的滚转、俯仰角后反解U1*/ 
	double U1 = rb->m * (G - az_des) / (cos(theta_des)*cos(phi_des));
	/*U1需要限幅*/ 
    U1 = clamp(U1, 0.0, 2.0 * rb->m * G);

	double ephi=phi_des-phi;
	double eth=theta_des-th;
	double epsi=c->psi_ref-psi;
	
	c->int_ephi += ephi*dt;
	c->int_eth += eth*dt;
	c->int_epsi += epsi*dt;
	
	c->int_ephi =clamp(c->int_ephi,-c->int_limit_att/c->ki_att-1e-9,c->int_limit_att/c->ki_att+1e-9);
	c->int_eth =clamp(c->int_eth,-c->int_limit_att/c->ki_att-1e-9,c->int_limit_att/c->ki_att+1e-9);
	c->int_epsi =clamp(c->int_epsi,-c->int_limit_att/c->ki_att-1e-9,c->int_limit_att/c->ki_att+1e-9);
    
    double phi_dd_des=c->kp_att*ephi+c->ki_att*c->int_ephi-c->kd_att*x[10];
    double theta_dd_des=c->kp_att*eth+c->ki_att*c->int_eth-c->kd_att*x[11];
    double psi_dd_des=c->kp_att*epsi+c->ki_att*c->int_epsi-c->kd_att*x[12];
    
    w->U1= U1 ;
    w->U2=rb->Ix*phi_dd_des/rb->l;
    w->U3=rb->Iy*theta_dd_des/rb->l;
    w->U4=rb->Iz*psi_dd_des/rb->l;
}
    
int main()
{
	Rigid6DOF rb;
	controller c;
	Wrench w;
	
	rb_init(&rb);
	ctrl_init(&c);
	memset(&w,0,sizeof(w));
	
	w.U1=rb.m*G;
	
	double const dt=0.001;
	double const T_end=20.0;
	
	printf("目标: (%.1f, %.1f, %.1f)\n\n",
           c.p_ref[0], c.p_ref[1], c.p_ref[2]);
    printf("  t      X       Y       Z     phi   theta  psi    U1\n");
    printf("----------------------------------------------------------\n");
    
    for (double t=0.0;t<=T_end+1e-9;t+=dt)
    {
    	if(fmod(t,1.0)<dt)
    	{
    		printf("%5.1f  %6.2f  %6.2f  %6.2f  %5.1f  %5.1f  %5.1f  %6.2f\n",
                   t, rb.x[1], rb.x[2], rb.x[3],
                   rb.x[4]*180/pi, rb.x[5]*180/pi, rb.x[6]*180/pi,
                   w.U1);
		}
		
		/* 先算控制，再积分 */
		ctrl_update(&c,&rb,&w,dt);
		rb_step(&rb,&w,dt);
	}
	printf("\n终点:  X=%.3f  Y=%.3f  Z=%.3f\n", rb.x[1], rb.x[2], rb.x[3]);
    printf("误差:  dX=%.3f  dY=%.3f  dZ=%.3f  m\n",
           rb.x[1]-c.p_ref[0],
           rb.x[2]-c.p_ref[1],
           rb.x[3]-c.p_ref[2]);
    return 0;
}
#endif
