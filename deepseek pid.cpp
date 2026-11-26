#include <stdio.h>
#include <math.h>
#include <string.h>

#define pi 3.141592653589793
#define G  9.8

/* 状态量 12个 */ 
enum {
    X=0, Y, Z,
    PHI, THETA, PSI,
    X_D, Y_D, Z_D,
    PHI_D, THETA_D, PSI_D,
    NX
};
/* 无人机的基本参数及状态量12个 */ 
typedef struct {
    double m, l, g;
    double Ix, Iy, Iz;
    double x[NX];
} Rigid6DOF;

/* 四个力 */ 
typedef struct {
    double U1, U2, U3, U4;
} Wrench;

/* 限幅工具 */
static double clamp(double v, double lo, double hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

/* ---------- 初始化 ---------- */
void rb_init(Rigid6DOF *rb) /*initial*/ 
{
    memset(rb, 0, sizeof(*rb)); /*memory set*/ 
    rb->m  = 1.0;
    rb->g  = G;
    rb->Ix = 0.0025;
    rb->Iy = 0.0025;
    rb->Iz = 0.0045;
    rb->l  = 0.17;
}

/* 运动模型 */ 
void rb_deriv(const Rigid6DOF *rb, const Wrench *w, double dx[NX])
{
    const double *x = rb->x;
    double phi = x[3], th = x[4], psi = x[5];
/* 可以填写x[3] 也可以填写x[PHI]*/ 
    double sphi = sin(phi), cphi = cos(phi);
    double sth  = sin(th),  cth  = cos(th);
    double spsi = sin(psi), cpsi = cos(psi);

    dx[X_D] = -(w->U1/rb->m) * (cpsi*sth*cphi + spsi*sphi);
    dx[Y_D] = -(w->U1/rb->m) * (spsi*sth*cphi - cpsi*sphi);
    dx[Z_D] =  rb->g - (w->U1/rb->m) * cth*cphi;

    dx[PHI_D]   = rb->l * w->U2 / rb->Ix;
    dx[THETA_D] = rb->l * w->U3 / rb->Iy;
    dx[PSI_D]   = rb->l * w->U4 / rb->Iz;

    dx[X] = x[X_D];   dx[Y] = x[Y_D];   dx[Z] = x[Z_D];
    dx[PHI]   = x[PHI_D];
    dx[THETA] = x[THETA_D];
    dx[PSI]   = x[PSI_D];
}

/* ---------- RK4 ---------- */
void rb_step(Rigid6DOF *rb, const Wrench *w, double dt)
{
    double x0[NX], k1[NX], k2[NX], k3[NX], k4[NX], xt[NX];
    Rigid6DOF tmp;
    memcpy(x0, rb->x, sizeof(x0));
    rb_deriv(rb, w, k1);
    for (int i = 0; i < NX; i++) xt[i] = x0[i] + 0.5*dt*k1[i];
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

/* ============================================================
 *  PID 控制器（级联：位置外环 + 姿态内环）
 * ============================================================ */
typedef struct {
    /*位置、偏航角参考值*/ 
	double p_ref[3];          /* 元素个数=3 */
    double psi_ref;

    /* 外环位置 PID */
    double kp_xy, ki_xy, kd_xy;
    double kp_z,  ki_z,  kd_z;

    /* 内环姿态 PID */
    double kp_att, ki_att, kd_att;
    double kp_yaw, ki_yaw, kd_yaw;

    /* 积分状态 */
    double int_ex, int_ey, int_ez; /* integrate */
    double int_ephi, int_eth, int_epsi;

    /* 积分限幅 */
    double int_limit_pos;     /* 位置积分限幅 */
    double int_limit_att;     /* 姿态积分限幅 */
} Controller;
/* 控制器初始化*/ 
void ctrl_init(Controller *c)
{
    c->p_ref[0] = 20.0;
    c->p_ref[1] =  0.0;
    c->p_ref[2] = -20.0;      /* 向上 20 m */
    c->psi_ref  =  0.0;

    /* 外环：先 PD，积分先关掉（ki=0），调好再开 */
    c->kp_xy = 1.0;   c->ki_xy = 0.0;   c->kd_xy = 2.0;
    c->kp_z  = 2.0;   c->ki_z  = 0.0;   c->kd_z  = 3.0;

    /* 内环：wn≈10, zeta≈1 */
    c->kp_att = 100.0; c->ki_att = 0.0; c->kd_att = 20.0;
    c->kp_yaw =  50.0; c->ki_yaw = 0.0; c->kd_yaw = 10.0;

    c->int_ex = c->int_ey = c->int_ez = 0.0;
    c->int_ephi = c->int_eth = c->int_epsi = 0.0;

    c->int_limit_pos = 2.0;    /* 位置积分项最多贡献 2 m/s2 */
    c->int_limit_att = 0.5;    /* 姿态积分项最多贡献 0.5 rad/s2 */
}
void w_init(Wrench *w,Rigid6DOF *rb) 
{
	w->U1=rb->m*rb->g;
	w->U2=0;
	w->U3=0;
	w->U4=0;
}
/* 控制器更新 */ 
void ctrl_update(Controller *c, const Rigid6DOF *rb, Wrench *w, double dt) 
{
    const double *x = rb->x;
    double g   = rb->g;
    double phi = x[PHI], th = x[THETA], psi = x[PSI];

    /* ================= 外环：位置 PID ================= */
    double ex = c->p_ref[0] - x[X];
    double ey = c->p_ref[1] - x[Y];
    double ez = c->p_ref[2] - x[Z];

    /* 积分（带限幅，防饱和） */
    c->int_ex += ex * dt;
    c->int_ey += ey * dt;
    c->int_ez += ez * dt;
    c->int_ex = clamp(c->int_ex, -c->int_limit_pos/c->ki_xy - 1e-9,
                                 c->int_limit_pos/c->ki_xy + 1e-9);
    c->int_ey = clamp(c->int_ey, -c->int_limit_pos/c->ki_xy - 1e-9,
                                 c->int_limit_pos/c->ki_xy + 1e-9);
    c->int_ez = clamp(c->int_ez, -c->int_limit_pos/c->ki_z  - 1e-9,
                                 c->int_limit_pos/c->ki_z  + 1e-9);

    /* 期望加速度：P + I + D（D 项对速度取负） */
    double ax_des = c->kp_xy*ex + c->ki_xy*c->int_ex - c->kd_xy*x[X_D];
    double ay_des = c->kp_xy*ey + c->ki_xy*c->int_ey - c->kd_xy*x[Y_D];
    double az_des = c->kp_z *ez + c->ki_z *c->int_ez - c->kd_z *x[Z_D];

    /* 加速度限幅（物理可行性） */
    const double a_xy = 0.5 * g;   /* 水平最大 ~4.9 m/s2 */
    const double a_z  = 0.4 * g;   /* 垂直最大 ~3.9 m/s2 */
    ax_des = clamp(ax_des, -a_xy, a_xy);
    ay_des = clamp(ay_des, -a_xy, a_xy);
    az_des = clamp(az_des, -a_z,  a_z);

    /* ============ 反解期望姿态 + 总推力 ============ */
    /* ### ψ≈0 时：
       ax = -(U1/m)·sθcφ
       ay =  (U1/m)·sφ
       az =  g - (U1/m)·cθcφ         */
    double U1_est = rb->m * g;
    double phi_des   = asin(clamp(rb->m*ay_des / U1_est, -0.99, 0.99));
    double theta_des = asin(clamp(-rb->m*ax_des / (U1_est*cos(phi_des)),
                                  -0.99, 0.99));

    /* 最大倾斜 30°，防止大误差时过度倾斜 */
    const double TILT = 30.0 * pi / 180.0;
    phi_des   = clamp(phi_des,   -TILT, TILT);
    theta_des = clamp(theta_des, -TILT, TILT);

    /* 用期望姿态精确计算 U1（高度补偿） */
    double U1 = rb->m * (g - az_des) / (cos(theta_des)*cos(phi_des));
    U1 = clamp(U1, 0.0, 2.0 * rb->m * g);

    /* ================= 内环：姿态 PID ================= */
    double ephi = phi_des   - phi;
    double eth  = theta_des - th;
    double epsi = c->psi_ref - psi;

    c->int_ephi += ephi * dt;
    c->int_eth  += eth  * dt;
    c->int_epsi += epsi * dt;
    c->int_ephi = clamp(c->int_ephi, -c->int_limit_att/c->ki_att - 1e-9,
                                     c->int_limit_att/c->ki_att + 1e-9);
    c->int_eth  = clamp(c->int_eth,  -c->int_limit_att/c->ki_att - 1e-9,
                                     c->int_limit_att/c->ki_att + 1e-9);
    c->int_epsi = clamp(c->int_epsi, -c->int_limit_att/c->ki_yaw - 1e-9,
                                     c->int_limit_att/c->ki_yaw + 1e-9);

    double phi_dd_des   = c->kp_att*ephi + c->ki_att*c->int_ephi
                        - c->kd_att*x[PHI_D];
    double theta_dd_des = c->kp_att*eth  + c->ki_att*c->int_eth
                        - c->kd_att*x[THETA_D];
    double psi_dd_des   = c->kp_yaw*epsi + c->ki_yaw*c->int_epsi
                        - c->kd_yaw*x[PSI_D];

    /* 角加速度 -> 力矩 -> 推力差（U = L/l） */
    w->U1 = U1;
    w->U2 = rb->Ix * phi_dd_des   / rb->l;
    w->U3 = rb->Iy * theta_dd_des / rb->l;
    w->U4 = rb->Iz * psi_dd_des   / rb->l;
}

/* ============================================================
 *  主函数：飞到 (20, 0, -20)
 * ============================================================ */
int main(void)
{
    Rigid6DOF rb;
    Controller ctrl;
    Wrench w;

    rb_init(&rb); /*应用初始化*/ 
    ctrl_init(&ctrl);
    w_init(&w,&rb); /*力和力矩清零*/ 
    

    const double dt    = 0.001;
    const double T_end = 30.0;

    printf("目标: (%.1f, %.1f, %.1f)\n\n",
           ctrl.p_ref[0], ctrl.p_ref[1], ctrl.p_ref[2]);
    printf("  t      X       Y       Z     phi   theta  psi    U1\n");
    printf("----------------------------------------------------------\n");

    for (double t = 0.0; t <= T_end + 1e-12; t += dt) {

        if (fmod(t, 1.0) < dt) {
            printf("%5.1f  %6.2f  %6.2f  %6.2f  %5.1f  %5.1f  %5.1f  %6.2f\n",
                   t, rb.x[X], rb.x[Y], rb.x[Z],
                   rb.x[PHI]*180/pi, rb.x[THETA]*180/pi, rb.x[PSI]*180/pi,
                   w.U1);
        }

        ctrl_update(&ctrl, &rb, &w, dt);   /* 先算控制，再积分 */
        rb_step(&rb, &w, dt);
    }

    printf("\n终点:  X=%.3f  Y=%.3f  Z=%.3f\n", rb.x[X], rb.x[Y], rb.x[Z]);
    printf("误差:  dX=%.3f  dY=%.3f  dZ=%.3f  m\n",
           rb.x[X]-ctrl.p_ref[0],
           rb.x[Y]-ctrl.p_ref[1],
           rb.x[Z]-ctrl.p_ref[2]);
    return 0;
}
