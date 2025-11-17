//axi
#include "grid/multigrid.h"
#if AXIS
# include "axi.h" // fixme: does not run with -catch
#endif
#if MOMENTUM
# include "momentum.h"
#else
#include "navier-stokes/centered.h"
#if CASE2
# define FILTERED 1
#endif
#if CLSVOF
# include "two-phase-clsvof.h"
#elif LEVELSET
# include "two-phase-levelset.h"
#else
# include "two-phase.h"
#endif
#endif
#if LEVELSET
# include "integral.h"
#else
# include "tension.h"
#endif
#if REDUCED
# include "reduced.h"
#endif

#ifndef LEVEL
# define LEVEL 8
#endif

/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. */

#if MOMENTUM  //动量
q.t[right] = dirichlet(0); //q为动量密度
q.t[left]  = dirichlet(0);
#else
u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);
#endif

int main() {

  /**
  The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
  $256\times 64$ grid points. */

  dimensions (nx = 4);  //设置x：y=4：1
  size (2 [1]);
  DT = 1. [0,1];
  init_grid (1 << LEVEL);
  
  /**
  Hysing et al. consider two cases (1 and 2), with the densities, dynamic
  viscosities and surface tension of fluid 1 and 2 given below. */

  rho1 = 1000.[0], mu1 = 10.;  // works also with rho1 = 1000.[-3,0,1]
#if CASE2
  rho2 = 1., mu2 = 0.1;
#else
  rho2 = 100., mu2 = 1.;
#endif

#if LEVELSET
  #if CASE2
  const scalar sigma[] = 1.96;
  #else
  const scalar sigma[] = 24.5;
  #endif
  d.sigmaf = sigma;
#else // !LEVELSET
  #if CASE2
  f.sigma = 1.96;
  #else
  f.sigma = 24.5;
  #endif
#endif // !LEVELSET
  
  /**
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  TOLERANCE = 1e-4 [*];
#if REDUCED
  G.x = -0.98;
  Z.x = 1.;
#endif
  run();
}

event init (t = 0) {

  /**
  The bubble is centered on (0.5,0) and has a radius of 0.25. */

#if LEVELSET
  foreach()
    d[] = sqrt (sq(x - 0.5) + sq(y)) - 0.25;
#else
  fraction (f, sq(x - 0.5) + sq(y) - sq(0.25));
#endif
}

/**
We add the acceleration of gravity. */

#if !REDUCED
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
}
#endif

/**
A utility function to check the convergence of the multigrid
solvers. */

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}

/**
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

event logfile (i++) {
  ddouble xb = 0., yb = 0., vb = 0., Vb = 0., Ub = 0.;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:vb)
                    reduction(+:Vb) reduction(+:Ub)) {
    double dv = (1. - f[])*dv();  // 气泡体积分数
    xb += x*dv;
    yb += y*dv;
    Vb += dv;
#if MOMENTUM
    Ub += q.x[]*dv / rho(f[]);
#else
    Ub += u.x[]*dv;
#endif
  }

  if (Vb < 1e-10) return 0;  // 防止除零

  xb /= Vb;  // 质心 x
  yb /= Vb;  // 质心 y
  Ub /= Vb;  // 平均上升速度

  // 等效直径
  double d = cbrt(6.0*Vb/M_PI);

  // 浮力 = 阻力（稳态假设）
  double g = 0.98;
  double Fb = rho1 * g * Vb;           // 浮力（向上）
  double Fd = Fb;                      // 阻力（向下，数值相等）

  // 投影面积 A = π d² / 4
  double A = M_PI * d*d / 4.0;

  // 阻力系数 Cd
  double Usq = Ub*Ub + 1e-20;
  double Cd = 2.0 * Fd / (rho1 * Usq * A);

  // 雷诺数 Re
  double Re = rho1 * fabs(Ub) * d / mu1;

  // 初始体积
  static double V0 = 0.;
  if (i == 0) {
    V0 = Vb;
    printf("# t    V/V0     y_b     U_b     d       Re       Cd       Fb\n");
  }

  // 输出：时间、归一化体积、质心y、速度、直径、Re、Cd、浮力
  printf("%g %g %g %g %g %g %g %g\n",
         t, Vb/V0, yb, Ub, d, Re, Cd, Fb);
#if !MOMENTUM
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
#endif
  putchar ('\n');
  fflush (stdout);
}

/**
At $t=3$ we output the shape of the bubble. */
FILE *ffmpeg_gif = NULL;
FILE *ffmpeg_mp4 = NULL;

// 替换你的 movies 事件
event movies (i += 10) {
  static int nf = 0;
  if (nf++ == 0) {
    // 启动 ffmpeg 管道（实时接收 PPM 图像流）
    ffmpeg_gif = popen(
      "ffmpeg -y -f image2pipe -vcodec ppm -r 20 -i - "
      "-loop 0 -vf 'scale=512:trunc(iw*0.5/2)*2' bubble.gif", "w");
    ffmpeg_mp4 = popen(
      "ffmpeg -y -f image2pipe -vcodec ppm -r 20 -i - "
      "-c:v libx264 -pix_fmt yuv420p -vf 'scale=512:trunc(iw*0.5/2)*2' bubble.mp4", "w");
  }

  // 正确调用：前两个是位置参数，其余是命名参数
  output_ppm (f, ffmpeg_gif, n = 512, min = 0, max = 1, linear = true);
  output_ppm (f, ffmpeg_mp4, n = 512, min = 0, max = 1, linear = true);
}

event end (t = end) {
  if (ffmpeg_gif) pclose(ffmpeg_gif);
  if (ffmpeg_mp4) pclose(ffmpeg_mp4);
  fprintf(stderr, "Animation generated: bubble.gif and bubble.mp4\n");
}
#if _GPU && SHOW
event display (i++) {
  output_ppm (f, fp = NULL, fps = 30, spread = -1, n = 1024);
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, fp = NULL, fps = 30, spread = 6, n = 1024);
}
#endif

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){5e-4,1e-3,1e-3}, LEVEL);
}
#endif


event stop(t=50.){
  printf("simution stop at 3s");
}
