/* ============================================================
 *  理想六自由度刚体模型（无扰动、无模型不匹配）
 *  坐标系：NED，机体系 x 前 / y 右 / z 下
 *
 *  状态(12)：
 *    [0..2]  位置 pn, pe, pd          (m)
 *    [3..5]  欧拉角 phi, theta, psi   (rad)
 *    [6..8]  机体系线速度 u, v, w     (m/s)
 *    [9..11] 机体系角速度 p, q, r     (rad/s)
 *
 *  输入(6)：机体系合外力 Fx,Fy,Fz（不含重力） + 合外力矩 L,M,N
 * ============================================================ */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.1415926
#define G  9.8

/* ---------- 状态索引 ---------- */
enum {
    PN = 0, PE, PD,
    PHI, THETA, PSI,
    U, V, W,
    P, Q, R,
    NX                       /* = 12 */
};

/* ---------- 对象 ---------- */
typedef struct {
    double m;                /* 质量 kg */
    double g;                /* 重力加速度 */
    double Ixx, Iyy, Izz;    /* 主惯量 kg·m2 */
    double Ixz;              /* 惯量积（对称机体常为 0） */
    double x[NX];            /* 状态 */
} Rigid6DOF;

/* ---------- 输入 ---------- */
typedef struct {
    double Fx, Fy, Fz;       /* 机体系合外力（不含重力） N */
    double L,  M,  N;        /* 机体系合外力矩 N·m */
} Wrench;

/* ---------- 初始化 ---------- */
void rb_init(Rigid6DOF *rb)
{
    memset(rb, 0, sizeof(*rb));
    rb->m   = 1.0;
    rb->g   = G;
    rb->Ixx = 0.0025;
    rb->Iyy = 0.0025;
    rb->Izz = 0.0045;
    rb->Ixz = 0.0;
}

/* ============================================================
 *  状态导数 dx = f(x, u)
 * ============================================================ */
void rb_deriv(const Rigid6DOF *rb, const Wrench *w, double dx[NX])
{
    const double *x = rb->x;

    double phi = x[PHI],   th = x[THETA], psi = x[PSI];
    double u   = x[U],     v  = x[V],     ww  = x[W];
    double p   = x[P],     q  = x[Q],     r   = x[R];

    double sphi = sin(phi),   cphi = cos(phi);
    double sth  = sin(th),    cth  = cos(th);
    double spsi = sin(psi),   cpsi = cos(psi);

    /* ---------- 平动动力学（机体系） ---------- */
    /* m·(v?_b + ω×v_b) = F_b + 重力在机体系投影 */
    dx[U] = r*v  - q*ww + w->Fx/rb->m - rb->g*sth;
    dx[V] = p*ww - r*u  + w->Fy/rb->m + rb->g*sphi*cth;
    dx[W] = q*u  - p*v  + w->Fz/rb->m + rb->g*cphi*cth;

    /* ---------- 转动动力学（欧拉方程） ---------- */
    double Ixx = rb->Ixx, Iyy = rb->Iyy, Izz = rb->Izz, Ixz = rb->Ixz;
    double Gam = Ixx*Izz - Ixz*Ixz;

    dx[P] = ( Ixz*(Ixx - Iyy + Izz)*p*q
            - (Izz*(Izz - Iyy) + Ixz*Ixz)*q*r
            + Izz*w->L + Ixz*w->N ) / Gam;

    dx[Q] = ( (Izz - Ixx)*p*r + Ixz*(r*r - p*p) + w->M ) / Iyy;

    dx[R] = ( (Ixx*(Ixx - Iyy) + Ixz*Ixz)*p*q
            - Ixz*(Ixx - Iyy + Izz)*q*r
            + Ixz*w->L + Ixx*w->N ) / Gam;

    /* ---------- 运动学：位置 ---------- */
    dx[PN] = u*cth*cpsi + v*(sphi*sth*cpsi - cphi*spsi)
                         + ww*(cphi*sth*cpsi + sphi*spsi);
    dx[PE] = u*cth*spsi + v*(sphi*sth*spsi + cphi*cpsi)
                         + ww*(cphi*sth*spsi - sphi*cpsi);
    dx[PD] = -u*sth + v*sphi*cth + ww*cphi*cth;

    /* ---------- 运动学：欧拉角 ---------- */
    double cth_safe = cth;
    if (fabs(cth_safe) < 1e-4)
        cth_safe = (cth_safe >= 0.0) ? 1e-4 : -1e-4;

    dx[PHI]   = p + (q*sphi + r*cphi) * tan(th);
    dx[THETA] = q*cphi - r*sphi;
    dx[PSI]   = (q*sphi + r*cphi) / cth_safe;
}

/* ---------- RK4 ---------- */
void rb_step(Rigid6DOF *rb, const Wrench *w, double dt)
{
    double x0[NX], k1[NX], k2[NX], k3[NX], k4[NX], xt[NX];
    Rigid6DOF tmp;

    memcpy(x0, rb->x, sizeof(x0));

    rb_deriv(rb, w, k1);

    for (int i = 0; i < NX; i++) xt[i] = x0[i] + 0.5*dt*k1[i];
    tmp = *rb; memcpy(tmp.x, xt, sizeof(xt));
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
 *  主函数：两个理想场景
 *    A. 悬停配平：Fz = -mg，其余为 0，看它一动不动
 *    B. 恒定滚转力矩：L = 0.01 N·m，看滚转角速度线性增长
 * ============================================================ */
int main(void)
{
    Rigid6DOF rb;
    rb_init(&rb);

    const double dt    = 0.001;
    const double T_end = 2.0;

    /* ==================== 场景 A：悬停配平 ==================== */
    printf("==== 场景 A：悬停配平 ====\n");
    rb_init(&rb);

    Wrench w_hover = { 0.0, 0.0, -rb.m * rb.g,  0.0, 0.0, 0.0 };

    printf("  t     pn      pe      pd    phi   theta   psi\n");
    for (double t = 0.0; t <= T_end + 1e-12; t += dt) {
        if (fmod(t, 0.5) < dt)
            printf("%5.2f  %7.4f %7.4f %7.4f  %5.2f %6.2f %6.2f\n",
                   t, rb.x[PN], rb.x[PE], rb.x[PD],
                   rb.x[PHI]*180/PI, rb.x[THETA]*180/PI, rb.x[PSI]*180/PI);
        rb_step(&rb, &w_hover, dt);
    }

    /* ==================== 场景 B：恒定滚转力矩 ==================== */
    printf("\n==== 场景 B：恒定滚转力矩 L = 0.01 N·m ====\n");
    rb_init(&rb);

    /* 先让 Z 向推力配平重力，同时施加恒定滚转力矩 */
    Wrench w_roll = { 0.0, 0.0, -rb.m * rb.g,  0.01, 0.0, 0.0 };

    printf("  t     phi(deg)   p(rad/s)   pn       pd\n");
    for (double t = 0.0; t <= T_end + 1e-12; t += dt) {
        if (fmod(t, 0.5) < dt)
            printf("%5.2f  %8.3f  %9.4f  %7.3f  %7.3f\n",
                   t, rb.x[PHI]*180/PI, rb.x[P], rb.x[PN], rb.x[PD]);
        rb_step(&rb, &w_roll, dt);
    }

    return 0;
}
