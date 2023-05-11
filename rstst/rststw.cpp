#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BM.h"

#define Np     1000 // 4の倍数であること;NP=4*r^2*lo
#define Nn     100
#define lowall 1
#define R   80.  //固定;
                //  //R=sqrt(np/4/lo);
// ,0.1より大きいこと;
#define tmax     5000 // 973.686//2*100たうとする;
#define tmaxlg   500 // 緩和時間は10たうとする;
#define dtlg     0.0001
#define dt       0.0001
#define temp     50.         // v0^2=2D/tau,ここではDを入れること;
#define dim      2           // 変えるときはEomを変えること;
#define cut      1.122462048 // 3.
#define skin     1.5
#define tau      50.
#define ensemble 1
// #define polydispersity 0.2 コードも変える;
#define folder_name "stwne4"
#define mgn         0.// Omega=omega/tau,ここではomegaを入れること;
//#define radios 1. //粒径の平均値を変えるときはヒストグラムの変え方も変えること:現在は1;
// v4:lohistをNpで割らなくした;
// v五:ディレクトリ名でRを先にした;
//stwv2:ティレク取りメイをloに、ヒストグラムの取り方を外からに;
//rststw:壁との相互作用をWCAに、countbound廃止、corrのカウント方法変更;
void tout_update(double *tout) { *tout *= 1.1; }

void ini_coord_circle(double (*x)[dim]) {
    double R2 = R - 0.5;
    double num_max = sqrt(Np / M_PI);
    double bitween = (R2 - 0.1) / num_max;
    int    n025 = Np * 0.25;
    int    i, j, k = 0;
    x[0][0] = 0.;
    x[0][1] = 0.;
    for (j = 1; j < num_max; j++) {
        x[k][0] = j * bitween;
        x[k][1] = 0.;
        x[k + n025][0] = -j * bitween;
        x[k + n025][1] = 0.;
        x[k + n025 * 2][0] = 0.;
        x[k + n025 * 2][1] = j * bitween;
        x[k + n025 * 3][0] = 0.;
        x[k + n025 * 3][1] = -j * bitween;
        k++;
    }
    for (j = 1; j < num_max; j++) {
        for (i = 1; i < num_max; i++) {
            if (i * bitween * i * bitween + j * j * bitween * bitween <
                R2 * R2) {
                x[k][0] = i * bitween;
                x[k][1] = j * bitween;
                x[k + n025][0] = -i * bitween;
                x[k + n025][1] = j * bitween;
                x[k + 2 * n025][0] = -i * bitween;
                x[k + 2 * n025][1] = -j * bitween;
                x[k + 3 * n025][0] = i * bitween;
                x[k + 3 * n025][1] = -j * bitween;
                k++;
                if (k >= Np * 0.25)
                    break;
            } else {
                continue;
            }
        }
        if (k >= Np * 0.25)
            break;
    }
}
void setcoord_cirwall(double (*x)[dim]) {
    double theta = 1. / (R * lowall);
    double Nwall = 2. * M_PI * R * lowall;
    for (int i = Np; i < Nwall; i++) {
        x[i][0] = R * cos(theta * i);
        x[i][1] = R * sin(theta * i);
    }
}
void set_diameter(double *a) {
    for (int i = 0; i < Np; i++)
        a[i] = 0.5;
}



void ini_array(double (*x)[dim]) {
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++)
            x[i][j] = 0.0;
}

void list_verlet(int (*list)[Nn], double (*x)[dim]) {
    double dx, dy, dr2;
    double thresh = cut + skin;
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < Nn; j++)
            list[i][j] = 0;

    for (int i = 0; i < Np; i++)
        for (int j = 0; j < Np; j++) {
            if (j > i) {
                dx = x[i][0] - x[j][0];
                dy = x[i][1] - x[j][1];
                //   dx-=L*floor((dx+0.5*L)/L);
                //   dy-=L*floor((dy+0.5*L)/L);
                dr2 = dx * dx + dy * dy;
                if (dr2 < thresh * thresh) {
                    list[i][0]++;
                    list[i][(int) list[i][0]] = j;
                }
            }
        }
}

void calc_force(double (*x)[dim], double (*f)[dim], double *a,
                int (*list)[Nn]) {
    double dx, dy, dr2, dUr, w2, w6, w12, aij;
    ini_array(f);

    for (int i = 0; i < Np; i++)
        for (int j = 1; j <= list[i][0]; j++) {
            dx = x[i][0] - x[list[i][j]][0];
            dy = x[i][1] - x[list[i][j]][1];

            dr2 = dx * dx + dy * dy;
            if (dr2 < cut * cut) {
                aij = (a[i] + a[list[i][j]]);
                w2 = aij * aij / dr2;
                w6 = w2 * w2 * w2;
                w12 = w6 * w6;
                dUr = (-48. * w12 + 24. * w6) / dr2 /* -12. * w12 / dr2*/;
                f[i][0] -= dUr * dx;
                f[list[i][j]][0] += dUr * dx;
                f[i][1] -= dUr * dy;
                f[list[i][j]][1] += dUr * dy;
            }
        }
}

void eom_aoup(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              double temp0, int (*list)[Nn], double (*F)[dim]) {
    double tauinv = dt / tau;
    double F0[dim],v0[2],fiw[2];
    double fluc = sqrt(temp0 / dt), ri,riw,aij,w2,w6,w12,dUr;
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {
        fiw[0]=0.;
        fiw[1]=0.;
////*force bitween wall;
        
        ri=sqrt(x[i][0]*x[i][0]+x[i][1]*x[i][1]);
        riw=R+0.5-ri;
        if(riw<cut){
            aij=0.5+a[i];
            w2=aij*aij/(riw*riw);
            w6=w2*w2*w2;
            w12=w6*w6;
            dUr=(-48.*w12+24.*w6)/(riw*ri);
            fiw[0]=dUr*x[i][0];
            fiw[1]=dUr*x[i][1];
        }
///till here;*/
        F0[0] = F[i][0];
        F0[1] = F[i][1];
        F[i][0] += (-F0[0] - F0[1] * mgn + fluc * gaussian_rand()) * tauinv;
        F[i][1] += (-F0[1] + F0[0] * mgn + fluc * gaussian_rand()) * tauinv;

        v0[0] = F[i][0] + f[i][0]+fiw[0];
        v0[1] = F[i][1] + f[i][1]+fiw[1];
        x[i][0] += v0[0] * dt;
        x[i][1] += v0[1] * dt;
        
    }
}

void ini_hist(double *hist, int Nhist) {
    for (int i = 0; i < Nhist; i++) {
        hist[i] = 0.;
    }
}

void make_v_thetahist(double (*x)[dim], double (*v)[dim], double(*hist),
                      double *hist2,double *lohist) {
    // lohist  と一緒に運用し、outputでv_theta[i]/lo[i];
    // v_thetaとomegaを算出、histがｖhist2がΩ;
    double v_t, dr,rint=(int)R,rsyou=R-rint,
                    bunbo = 1. / ( ensemble * floor(tmax / dt));
    int histint;
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        v_t = (x[i][0] * v[i][1] - x[i][1] * v[i][0]) / (dr*dr);
        if (dr < R) {
            histint=(int) floor(abs(dr-rsyou));
            hist[histint] += v_t * bunbo*dr;
            hist2[histint] += v_t * bunbo ;
            lohist[histint]+=bunbo;
        }
    }
}
void make_vt1hist(double (*x)[dim], double (*v)[dim], double(*hist)) {
    // lohist  と一緒に運用し、outputでv_theta[i]/lo[i];
    double v_t, dr, rint = (int) R, rsyou = R - rint,
                    bunbo = 1. / (ensemble * floor(tmax / dt));
    int histint;
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        v_t = (x[i][0] * v[i][1] - x[i][1] * v[i][0]) / dr;
        if (dr < R)
            hist[(int) floor(abs(dr - rsyou))] += v_t * bunbo;
    }
}

void make_lo1hist(double (*x)[dim], double(*hist)) {
    double v_t, dr, rint = (int) R, rsyou = R - rint,
                    bunbo = 1. / (ensemble * floor(tmax / dt));
    int histint;
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        if (dr < R)
            hist[(int) floor(abs(dr - rsyou))] += bunbo;
    }
}
void make_lohistx0(double (*x)[dim], double(*hist),
                   double (*x0)[dim]) { // 関数を使うときは使用回数で割ること
    double v_t, dr, rint = (int) R, rsyou = R - rint,
                    bunbo = 1. / (ensemble * floor(tmax / dt));
    int histint;
    for (int i = 0; i < Np; i++) {
        dr = sqrt((x[i][0] - x0[i][0]) * (x[i][0] - x0[i][0]) +
                  (x[i][1] - x0[i][1]) * (x[i][1] - x0[i][1]));
        if (dr <= R)
            hist[(int) floor(abs(dr - rsyou))] += 1.;
    }
}
void calc_diff(double (*x)[dim], double (*v)[dim], double *v_thetahist,
               double *lohist, double *diffv_theta, double *difflo) {
    double dr, v_theta, Rhistinv = 1., bub = 1. / Np, dif;
    double rint = (int) R, rsyou = R - rint;
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        if (dr < R) {
            v_theta = (x[i][0] * v[i][1] - x[i][1] * v[i][0]) / dr;
            dif = (dr - lohist[(int) floor(abs(dr - rsyou))]);
            difflo[(int) floor(abs(dr - rsyou))] += dif * dif / Np;
            dif = v_theta - v_thetahist[(int) floor(abs(dr - rsyou))] /
                                lohist[(int) floor(abs(dr - rsyou))];
            diffv_theta[(int) floor(abs(dr - rsyou))] +=
                dif * dif / lohist[(int) floor(abs(dr - rsyou))];
        }
    }
}
void out_gosahist(double (*x)[dim], double (*v)[dim]) {
    int Nphist = (int) (R  + 1.);
     double v_thetahist[Nphist],
        lohist[Nphist], diffv_theta[Nphist], difflo[Nphist];
    double Mgn = mgn / tau, v0 = temp / tau;

    ini_hist(lohist, Nphist);
    ini_hist(v_thetahist, Nphist);
    ini_hist(difflo, Nphist);
    ini_hist(diffv_theta, Nphist);
    make_lo1hist(x, lohist);
    make_vt1hist(x, v, v_thetahist);
    calc_diff(x, v, v_thetahist, lohist, diffv_theta, difflo);
    char          filename[128];
    double        bitthist = R / Nphist,ryou=R-(int)R;
    std::ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/lzdif_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Nphist; i++) {
        if (lohist[i] != 0.) {
            file << (i + 0.5+ryou)  << "\t" << (v_thetahist[i] / lohist[i])
                 << "\t" << diffv_theta[i] << std::endl;

        } else {
            file << i +0.5*ryou << "\t" << 0 << "\t" << 0 << std::endl;
        }
    }
    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/lodif_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename, std::ios::app);
    for (int i = 0; i < Nphist; i++) {
        file << (i + 0.5+ryou)   << "\t"
             << lohist[i] / (2. * M_PI * (i + 0.5+ryou) ) << "\t"
             << difflo[i] << std::endl;
    }
    file.close();
}

void output(int k, double (*v)[dim], double (*x)[dim], int l) {
    char   filename[128];
    double        Mgn = mgn / tau;
    std::ofstream file;
    sprintf(filename,
            "./%s_coorlo%.2fv0%.1ftau%.3fm%.3f/"
            "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), tau, Mgn, l);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; i++) {
        file << k * dt << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0]
             << "\t" << v[i][1]  << std::endl;
    }
    file.close();
}
void out_setup() {
    char          filename[128];
    double        Mgn = mgn / tau, v0 = temp / tau;
    std::ofstream file;
    sprintf(filename, "./%slo%.2fv0%.1ftau%.3fm%.3f/setupr%fm%f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, R, Mgn);
    file.open(filename, std::ios::app); // append

    file << "dt=" << dt << std::endl;
    file << "cut" << cut << std::endl;
    file << "skin" << skin << std::endl;
    file << "Nn" << Nn << std::endl;
    file << "Np=" << Np << std::endl;
    file << "tmax=" << tmax << std::endl;
    file << "tmaxlg=" << tmaxlg << std::endl;
    file << "temp=" << temp << std::endl;
    file << "ens=" << ensemble << std::endl;
    file << "type=" << 2 << std::endl;
    file<<"2DkaraD"<<std::endl;
    file<<"壁はWCA"<<std::endl;
    file.close();
}

void outputhist(double *hist, int counthistv_theta, double *lohist,
                double *hist2) {
    char          filename[128];
    double        v_theta = 0.;
    double        bitthist = 1.;
    int Nphist = (int) (R  + 1.);
     double Mgn = mgn / tau,
        v0 = temp / tau;
        double rsyou=R-(int)R;
    std::ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/v_thetahist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    
    if (lohist[0] != 0.) {
        file << ( rsyou + 1.)*0.5 << "\t" << (hist[0] / lohist[0]) << std::endl;

        v_theta += hist[0] / Np;
    } else {
        file << (rsyou + 1.) * 0.5 << "\t" << 0 << std::endl;
    }
    for (int i = 1; i < Nphist; i++) {
        if (lohist[i] != 0.) {
            file << i +rsyou+0.5 << "\t" << (hist[i] / lohist[i]) << std::endl;

            v_theta += hist[i] /Np;
        } else {
            file << i +rsyou+0.5 << "\t" << 0 << std::endl;
        }
    }
    file.close();
    sprintf(filename, "./%slo%.2fv0%.1ftau%.3fm%.3f/v_theta_lo%.3f_tau%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau);
    file.open(filename, std::ios::app); // append
    file << tau << "\t" << Mgn <<  "\t" << R << "\t" << v_theta
         << std::endl;
    file << tau << "\t" << Mgn <<  "\t" << R << "\t" << v_theta
         << std::endl;

    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/omegahist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn,
            Np * 0.25 / (R * R), tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    
    if (lohist[0] != 0.) {
        file << (rsyou + 1.) * 0.5 << "\t" << (hist2[0] / lohist[0])
             << std::endl;

    } else {
        file << (rsyou + 1.) * 0.5 << "\t" << 0 << std::endl;
    }
    for (int i = 1; i < Nphist; i++) {
        if (lohist[i] != 0.) {
            file << i +rsyou+0.5 << "\t" << (hist2[i] / lohist[i]) << std::endl;

        } else {
            file << i +rsyou+0.5 << "\t" << 0 << std::endl;
        }
    }
    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/lohist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    file << (rsyou + 1.) * 0.5 << "\t" << (lohist[0] / (4. * (rsyou + 1.)))
         << std::endl;
    for (int i = 1; i < Nphist; i++) {
        file << i + 0.5 +rsyou<< "\t"
             << (lohist[i] / (8. * (i + 0.5+rsyou))) << std::endl;
    }
    file.close();
}

void outtuibi(double (*x)[dim], double t, double (*v)[dim], int cu) {
    char   filename[128];
    double v_theta;
    double Mgn = mgn / tau;
    v_theta = x[cu][0] * v[cu][1] - x[cu][1] * v[cu][0];

    std::ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/coordtui_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), tau, Mgn);
    file.open(filename, std::ios::app); // append

    file << t << "\t" << x[cu][0] << "\t" << x[cu][1] << "\t" << v_theta
         << std::endl;

    file.close();
}

void update(double (*x_update)[dim], double (*x)[dim]) {
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++)
            x_update[i][j] = x[i][j];
}

void calc_disp_max(double *disp_max, double (*x)[dim],
                   double (*x_update)[dim]) {
    double dx, dy;
    double disp;
    for (int i = 0; i < Np; i++) {
        dx = x[i][0] - x_update[i][0];
        dy = x[i][1] - x_update[i][1];

        disp = dx * dx + dy * dy;
        if (disp > *disp_max)
            *disp_max = disp;
    }
}

void auto_list_update(double *disp_max, double (*x)[dim],
                      double (*x_update)[dim], int (*list)[Nn]) {
    static int count = 0;
    count++;
    calc_disp_max(&(*disp_max), x, x_update);
    if (*disp_max > skin * skin * 0.25) {
        list_verlet(list, x);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        *disp_max = 0.0;
        count = 0;
    }
}


void calc_corr(double (*x)[dim], double (*x0)[dim], double (*v0)[dim],
               double (*v)[dim], double *xcor, double *vcor, int k,double *msd) {
    double bunbo = 1. / (Np * ensemble),dr;
    for (int i = 0; i < Np; i++) {
        for (int j = 0; j < dim; j++) {
            xcor[k] += x0[i][j] * x[i][j] * bunbo;
            vcor[k] += v0[i][j] * v[i][j] * bunbo;
            dr = x[i][j] - x0[i][j];
            msd[k] += dr * dr * bunbo;
            
        }
    }
}

void outputcorr(double *msd, double *vcor, double *t, int countout,
                double *msd2) { // takusannkaeru
    char          filename[128];
    double        v_theta;
    double        Mgn = mgn / tau, v0 = temp / tau;
    std::ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/xcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; i++) {
        file << t[i] << "\t" << msd[i] << std::endl;
    }
    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/vcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; i++) {
        file << t[i] << "\t" << vcor[i] << std::endl;
    }
    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/msd_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; i++) {
        file << t[i] << "\t" << msd2[i] << std::endl;
    }
    file.close();
}

void outv_thetat(double *hist, int counthistv_theta, double tout) {
    char          filename[128];
    double        v_theta = 0.;
    double        Mgn = mgn / tau;
    int           Nphist = (int) (R + 1.);
    std::ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/v_thetatime_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), tau, Mgn);
    file.open(filename, std::ios::app); // append
    double bunbo = (ensemble * floor(tmax / dt)) / (counthistv_theta * Np);
    for (int i = 0; i < Nphist; i++) {

        v_theta += hist[i] * bunbo;
    }
    file << tout << "\t" << v_theta << std::endl;
    file.close();
}
void ini_test(double (*x)[dim], double (*v)[dim]) {
    for (int i = 0; i < Np; i++) {
        x[i][0] = -1.;
        x[i][1] = 0.;
        v[i][0] = 0.;
        v[i][1] = 2.;
    }
}

void calc_bounddiff(double *bound, double *answer) {
    double dif;
    for (int i = 0; i < Np; i++)
        answer[0] += bound[i] / Np;
    for (int i = 0; i < Np; i++) {
        dif = answer[0] - bound[i];
        answer[1] += dif * dif / Np;
    }
}

int main() {
    double x[Np][dim], x_update[Np][dim], v[Np][dim], f[Np][dim], a[Np],
        F[Np][dim], x0[Np][dim],v0[Np][dim];
    int    list[Np][Nn];
    int    counthistv_theta = 0, countout = 0;
    int    Nphist = (int) (R + 1.);
    double hist[Nphist], lohist[Nphist], hist2[Nphist];
    double tout = 0.01, toutcoord = 0, U, disp_max = 0.0;

    int j = 0, k = 0, kcoord = 0;
    set_diameter(a);
    ini_coord_circle(x);
    ini_array(v);
    ini_array(F);
    ini_array(f);
    ini_hist(hist, Nphist);
    ini_hist(lohist, Nphist);
    ini_hist(hist2, Nphist);
    // ini_test(x,v);
    double Mgn = mgn / tau;
    char   foldername[128];
    sprintf(foldername, "%slo%.2fv0%.1ftau%.3fm%.3f", folder_name,
            Np * 0.25 / (R * R), temp / tau, tau, Mgn);
    const char *fname = foldername;
    mkdir(fname, 0777);
    char foldername2[128];
    sprintf(foldername2, "%s_coorlo%.2fv0%.1ftau%.3fm%.3f", folder_name,
            Np * 0.25 / (R * R), temp / tau, tau, Mgn);
    const char *fname2 = foldername2;
    mkdir(fname2, 0777);
    out_setup();
    
    while (tout  < tmax) {

            tout_update(&tout);

            countout++;
            
    }
    countout+=5;
    double msd[countout], t[countout], vcor[countout], msd2[countout];
    ini_hist(msd, countout);
    ini_hist(t, countout);
    ini_hist(msd2, countout);
    ini_hist(vcor, countout);

    j = 0;
    double ttemp = 5. * temp;
    if (ttemp / tau < 5)
        ttemp = 5. * tau;
    j = 0;
    while (j * dtlg < 10) {
        j++;
        auto_list_update(&disp_max, x, x_update, list);
        eom_aoup(v, x, f, a, ttemp, list, F);
        
    }
 
    j = 0;
    while (j * dt < tmaxlg) {
        j++;
        auto_list_update(&disp_max, x, x_update, list);
        eom_aoup(v, x, f, a, temp, list, F);
    }

    double r2test, r2testmax = 0.;
    int    ituibi = 0,tauch=tau/dt,tmaxch=tmax/dt;
    for (int ko = 0; ko < Np; ko++) {
        r2test = x[ko][0] * x[ko][0] + x[ko][1] * x[ko][1];
        if (r2test > r2testmax)
            r2testmax = r2test;
        ituibi = ko;
    }
    for (int i = 0; i < ensemble; i++) {
        for (int xnp = 0; xnp < Np; xnp++) {

            for (int xdim = 0; xdim < dim; xdim++) {
                x0[xnp][xdim] = x[xnp][xdim];
                v0[xnp][xdim] = v[xnp][xdim];
            }
        }
        j = 0;
        tout = 0.01;
        toutcoord = 0.;
        k = 0;
        kcoord = 0;

        calc_corr(x, x0, v0, v, msd, vcor, kcoord,msd2);
        t[0] = 0.;
        kcoord++;

        out_gosahist(x, v);
        output(j, v, x, k);
        k++;
        while (j  < tmaxch) {
            j++;
            auto_list_update(&disp_max, x, x_update, list);
            eom_aoup(v, x, f, a, temp, list, F);
            make_v_thetahist(x, v, hist, hist2,lohist);

            if (j  >= toutcoord) {
                output(j, v, x, k);
                outtuibi(x, toutcoord, v, ituibi);
                toutcoord += tauch;
                k++;
            }
            if (j * dt >= tout) {
                calc_corr(x, x0, v0, v, msd, vcor, kcoord,msd2);
                t[kcoord] = j * dt;
                kcoord++;
                tout_update(&tout);
                // std::cout<<k;
            }
        }
    }

    int counthazure = 0;
    for (int i = 0; i < Np; i++) {
        if (x[i][0] * x[i][0] + x[i][1] * x[i][1] > R * R)
            counthazure++;
    }
    char filename[128];

    std::ofstream file;
    
    sprintf(filename, "./%slo%.2fv0%.1ftau%.3fm%.3f/kekkalo%.3fm%.3f.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), Mgn);
    file.open(filename, std::ios::app); // append

    file  << counthistv_theta << " " << counthazure << std::endl;
    file.close();
    
    outputhist(hist, counthistv_theta, lohist, hist2);
    outputcorr(msd, vcor, t, countout, msd2);
    std::cout << "done" << std::endl;
    return 0;
}
