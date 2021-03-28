//!===========================================================================================!
//!===========================================================================================!
//! Numerical Algorithms Collections In Fortran 77 Written by Liu Fei Yu & Li                 !
//! Translate into C++ for easy access by Bao Biwen                                           !
//! Created on 2021.03.25 with Version 1.0                                                    !
//! Chapter 3 Numerical intergration                                                          !
//!===========================================================================================!
//!===========================================================================================!

#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const double pi = 3.14159265358;

//矩形法
//为方便起见 被积函数为sin(x) 积分范围为0到pi 根据定积分和求导只是改值为2
//利用矩形法来求该段积分
void juxing_inter()
{
    const int n = 500;//区间数
    double lower_limit = 0.0;
    double upper_limit = pi;
    double dx = (upper_limit-lower_limit)/n;//如果是n x_arr 从0到3.13XXX 到不了3.1415 n-1可以确保积分上下限
    double x_arr[n] = {0.0};
    for (int j=0;j<n;j++)
        {
            x_arr[j]=lower_limit+j*dx;
        }
//    for (int j=0;j<n;j++)
//        {
//            cout << x_arr[j] << endl;
//        }
    double sum = 0.0;//积分值
    for (int j=0;j<n;j++)
        {
            sum = sum + dx*sin(x_arr[j]);//sin(x)的x是弧度
        }
    cout << "利用矩形法 积分sin(x)在区间0到pi之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}


//矩形法 函数指针版本
//为方便起见 被积函数为sin(x) 积分范围为0到pi 根据定积分和求导只是改值为2
//利用矩形法来求该段积分

double func(double x)
{
    return sin(x);
}

void juxing_inter_fp(double (*func)(double x), double lower_limit, double upper_limit)
{
    const int n = 500;//区间数
//    double lower_limit = 0.0;
//    double upper_limit = pi;
    double dx = (upper_limit-lower_limit)/n;//如果是n x_arr 从0到3.13XXX 到不了3.1415 n-1可以确保积分上下限
    double x_arr[n] = {0.0};
    for (int j=0;j<n;j++)
        {
            x_arr[j]=lower_limit+j*dx;
        }
//    for (int j=0;j<n;j++)
//        {
//            cout << x_arr[j] << endl;
//        }
    double sum = 0.0;//积分值
    for (int j=0;j<n;j++)
        {
            sum = sum + dx*(*func)(x_arr[j]);//sin(x)的x是弧度
        }
    cout << "利用矩形法(函数指针版) 积分sin(x)在区间0到pi之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}

//梯形法
//为方便起见 被积函数为sin(x) 积分范围为0到pi 根据定积分和求导只是改值为2
//利用矩形法来求该段积分
void tixing_inter()
{
    const int n = 500;//区间数
    double lower_limit = 0.0;
    double upper_limit = pi;
    double dx = (upper_limit-lower_limit)/(n-1);//如果是n x_arr 从0到3.13XXX 到不了3.1415 n-1可以确保积分上下限
    double x_arr[n] = {0.0};
    for (int j=0;j<n;j++)
        {
            x_arr[j]=lower_limit+j*dx;
        }
//    for (int j=0;j<n;j++)
//        {
//            cout << x_arr[j] << endl;
//        }
    double sum = 0.0;//积分值
    for (int j=0;j<n-1;j++)
        {
            double dsum = (sin(x_arr[j])+sin(x_arr[j+1]))*dx*0.5;
            sum = sum + dsum;//sin(x)的x是弧度
        }
    cout << "利用梯形法 积分sin(x)在区间0到pi之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}


//梯形法 函数指针版本
//为方便起见 被积函数为sin(x) 积分范围为0到pi 根据定积分和求导只是改值为2
//利用矩形法来求该段积分
void tixing_inter_fp(double (*func)(double x), double lower_limit, double upper_limit)
{
    const int n = 500;//区间数
//    double lower_limit = 0.0;
//    double upper_limit = pi;
    double dx = (upper_limit-lower_limit)/(n-1);//如果是n x_arr 从0到3.13XXX 到不了3.1415 n-1可以确保积分上下限
    double x_arr[n] = {0.0};
    for (int j=0;j<n;j++)
        {
            x_arr[j]=lower_limit+j*dx;
        }
//    for (int j=0;j<n;j++)
//        {
//            cout << x_arr[j] << endl;
//        }
    double sum = 0.0;//积分值
    for (int j=0;j<n-1;j++)
        {
            double dsum = ((*func)(x_arr[j])+(*func)(x_arr[j+1]))*dx*0.5;
            sum = sum + dsum;//sin(x)的x是弧度
        }
    cout << "利用梯形法(函数指针版本) 积分sin(x)在区间0到pi之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}



//定步长辛普森积分法
void xps_dbc_inter()
{
    const int n = 500;//区间数
    double lower_limit = 0.0;
    double upper_limit = 1.0;
    double h =0.0;
    h = (upper_limit-lower_limit)/(n-1);
    double x_arr[n] = {0.0};
    for (int j=0;j<n;j++)
        {
            x_arr[j]=lower_limit+j*h;
        }
//    for (int j=0;j<n;j++)
//        {
//            cout << x_arr[j] << endl;
//        }
    double sum=0.0;
    int flag = 1; //奇数偶数项标记
    for (int j=1;j<n-1;j++)
        {
            if(flag==1)
            {
                sum = sum + 4.0*cos(x_arr[j]);
                flag = 0;
                continue;//很重要 不然就会执行下面的flag==0
            }
            if(flag==0)
            {
                sum = sum + 2.0*cos(x_arr[j]);
                flag=1;
                continue;
            }

        }
    sum = sum + cos(lower_limit) + cos(upper_limit);
    sum = sum*h/3.0;
    cout << "利用梯形法 积分cos(x)在区间0到1之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}


//定步长辛普森积分法
//函数指针版本

double func2(double x)
{
    return cos(x);
}

void xps_dbc_inter_fp(double (*func)(double x), double lower_limit, double upper_limit)
{
    const int n = 500;//区间数
//    double lower_limit = 0.0;
//    double upper_limit = 1.0;
    double h =0.0;
    h = (upper_limit-lower_limit)/(n-1);
    double x_arr[n] = {0.0};
    for (int j=0;j<n;j++)
        {
            x_arr[j]=lower_limit+j*h;
        }
//    for (int j=0;j<n;j++)
//        {
//            cout << x_arr[j] << endl;
//        }
    double sum=0.0;
    int flag = 1; //奇数偶数项标记
    for (int j=1;j<n-1;j++)
        {
            if(flag==1)
            {
                sum = sum + 4.0*(*func)(x_arr[j]);
                flag = 0;
                continue;//很重要 不然就会执行下面的flag==0
            }
            if(flag==0)
            {
                sum = sum + 2.0*(*func)(x_arr[j]);
                flag=1;
                continue;
            }

        }
    sum = sum + (*func)(lower_limit) + (*func)(upper_limit);
    sum = sum*h/3.0;
    cout << "利用梯形法(函数指针版本) 积分cos(x)在区间0到1之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}


//变步长辛普森积分法
double f(double x)
{
    return log(1+x)/(1+x*x);//写要求辛普森积分的函数
}

double simpson(double L, double R)
{//三点辛普森积分法，要求f(x)是全局函数
    double mid = (L + R) / 2.0;
    return (f(L) + 4.0 * f(mid) + f(R)) * (R - L) / 6.0;
}

double xps_bbc_integral(double L, double R, double Eps)
{//自适应辛普森积分递归过程
    double mid = (L + R) / 2.0;
    double ST = simpson(L, R);
    double SL = simpson(L, mid);
    double SR = simpson(mid, R);
    if(fabs(SL + SR - ST) <= 15.0 * Eps)
    {
        return SL + SR + (SL + SR - ST) / 15.0;//直接返回结果
    }

    return xps_bbc_integral(L, mid, Eps/2.0) + xps_bbc_integral(mid, R, Eps/2.0);//对半划分区间
}



//梯形二重积分法
//积分e^(x^2+y^2)在区间x从0到1 y从-sqrt(1-x^2)到sqrt(1-x^2)之间的值 答案为2.699022172
//其实如果积分区间对称 直接2倍0到sqrt(1-x^2)的积分值更准确
double f1d(double low, double up)
{
    const int n = 300;//区间数
    double dx = (up-low)/(n-1);
    double x_arr[n] = {0.0};
    for (int j=0;j<n;j++)
        {
            x_arr[j]=low+j*dx;
        }
    double sum = 0.0;//积分值
    for (int j=0;j<n-1;j++)
        {
            double dsum = (exp(x_arr[j]*x_arr[j])+exp(x_arr[j+1]*x_arr[j+1]))*dx*0.5;
            sum = sum + dsum;//sin(x)的x是弧度
        }
    return sum;
}

void tixing_inter_2d()
{
    const int n = 300;//区间数
    double lower_limit = 0.0;
    double upper_limit = 1.0;
    double dx = (upper_limit-lower_limit)/(n-1);//如果是n x_arr 从0到3.13XXX 到不了3.1415 n-1可以确保积分上下限
    double x_arr[n] = {0.0};
    for (int j=0;j<n;j++)
        {
            x_arr[j]=lower_limit+j*dx;
        }
//    for (int j=0;j<n;j++)
//        {
//            cout << x_arr[j] << endl;
//        }
    double sum = 0.0;//积分值
    for (int j=0;j<n-1;j++)
        {
            double part1 = exp(x_arr[j]*x_arr[j])*f1d(-sqrt(1.0-x_arr[j]*x_arr[j]),sqrt(1.0-x_arr[j]*x_arr[j]));
            double part2 = exp(x_arr[j+1]*x_arr[j+1])*f1d(-sqrt(1.0-x_arr[j+1]*x_arr[j+1]),sqrt(1.0-x_arr[j+1]*x_arr[j+1]));
            double dsum = (part1+part2)*dx*0.5;
            sum = sum + dsum;//sin(x)的x是弧度
        }
    cout << "利用梯形法 积分二重e^(x^2+y^2)在区间x从0到1 y从-sqrt(1-x^2)到sqrt(1-x^2)之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}


//变步长辛普森积分法 二重积分版本
//积分e^(x^2+y^2)在区间x从0到1 y从-sqrt(1-x^2)到sqrt(1-x^2)之间的值 答案为2.699022172
//跟用普通梯形法比较 改办法可以控制Eps到很小范围 比较精准
double ff(double x)
{
    return exp(x*x);//内层积分
}

double simpson_in(double L, double R)
{//三点辛普森积分法，要求f(x)是全局函数
    double mid = (L + R) / 2.0;
    return (ff(L) + 4.0*ff(mid) + ff(R)) * (R-L)/6.0;
}


double f2d(double L2, double R2, double Eps)
{
//    double L2 = -sqrt(1-x*x);
//    double R2 = sqrt(1-x*x); //内层上下限积分
    double mid = (L2 + R2) / 2.0;
    double ST = simpson_in(L2, R2);
    double SL = simpson_in(L2, mid);
    double SR = simpson_in(mid, R2);
    if(fabs(SL + SR - ST) <= 15.0 * Eps)
    {
        return SL + SR + (SL + SR - ST) / 15.0;//直接返回结果
    }

    return f2d(L2, mid, Eps/2.0) + f2d(mid, R2, Eps/2.0);//对半划分区间
}

double simpson2d(double L, double R, double Eps)
{//三点辛普森积分法，要求f(x)是全局函数
    double mid = (L + R) / 2.0;
    double part1 = exp(L*L)*f2d(-sqrt(1-L*L),sqrt(1-L*L),Eps);
    double part2 = 4.0*exp(mid*mid)*f2d(-sqrt(1-mid*mid),sqrt(1-mid*mid),Eps);
    double part3 = exp(R*R)*f2d(-sqrt(1-R*R),sqrt(1-R*R),Eps);
    return (part1 + part2 + part3)*(R-L)/6.0;

}

double xps_bbc_integral2d(double L, double R,double Eps) //L1 R1为第一层积分上下限
{//自适应辛普森积分递归过程
    double mid = (L + R) / 2.0;
    double ST = simpson2d(L, R, Eps);
    double SL = simpson2d(L, mid, Eps);
    double SR = simpson2d(mid, R, Eps);
    if(fabs(SL + SR - ST) <= 15.0 * Eps)
    {
        return SL + SR + (SL + SR - ST) / 15.0;//直接返回结果
    }

    return xps_bbc_integral2d(L, mid, Eps/2.0) + xps_bbc_integral2d(mid, R, Eps/2.0);//对半划分区间
}


//高斯方法求积
void gauss_inter_fp(double (*func)(double x), double lower_limit, double upper_limit)
{
    const int n = 500;//区间数
//    double lower_limit = 0.0;
//    double upper_limit = 1.0;
    int m = 0.5*n-1;
    double h =0.0;
    h = (upper_limit-lower_limit)/n;
    double sum = 0.0;
    for (int j=0; j<=m; j++)
        {
            double part1 = lower_limit + h*(1.0-1.0/sqrt(3.0)+2.0*j);
            double part2 = lower_limit + h*(1.0+1.0/sqrt(3.0)+2.0*j);
            double part3 = (*func)(part1)+(*func)(part2);
            sum = sum + part3*h;
        }
    cout << "利用高斯法 积分cos(x)在区间0到1之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}


//高斯勒让德求积法
//n=2 x0=-0.5773502692 x1=0.5773502692 A0=A1=1.0
//n=3 x0=-0.7745966692 x1=0.0 x2=0.7745966692 A0=A2=0.5555555556 A1=0.8888888889
//n=4 x0=-0.8611363116 x1=-0.3399810436 x2=0.3399810436 x3=0.8611363116 A0=A3=0.3478548541 A1=A2=0.6521451549
//n=5 x0=-0.9061798459 x1=-0.5384693101 x3=0.0 x4=0.5384693101 x5=0.9061798459 A0=A4=0.2369268851 A1=A3=0.4786286705 A2=0.5688888889
//高斯方法求积
//这里采用n=5

double func3(double x)
{
    return x*x+sin(x);
}

void gauss_lerande_inter_fp(double (*func)(double x), double lower_limit, double upper_limit)
{
    double Ak[5]={0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851};
    double xk[5]={-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459};

    double sum = 0.0;
    for (int j=0; j<5; j++)
        {
            double part1 = 0.5*(upper_limit-lower_limit)*xk[j]+0.5*(upper_limit+lower_limit);
            sum = sum + Ak[j]*(*func)(part1);
        }
    sum = sum*0.5*(upper_limit-lower_limit);
    cout << "利用高斯勒让德法 积分x^2+sin(x)在区间2.5到8.4之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}

//切比雪夫求积公式
//xk=cos[(2k+1)/(2n+2)*pi] k=0,1,...n
//Ak=pi/(n+1)

double func4(double x)
{
    return exp(x);
}

//这里默认的上下限是-1.0 到 1.0
//如果实际的是a到b 需要进行变换到-1.0到1.0再进行积分
//类似之前高斯积分的变换
void qiebi_inter_fp(double (*func)(double x), double lower_limit, double upper_limit)
{
    const int n = 20;
    double Ak[n]={0.0};
    double xk[n]={0.0};
    for (int j=0;j<n;j++)
        {
            Ak[j]=pi/(n+1);
            xk[j]=cos((2.0*j+1.0)/(2.0*n+2.0)*pi);//注意整除问题 要带个1.0 保证结果是double
        }

//        for (int j=0;j<n;j++)
//        {
//           cout << Ak[j] << '\t' << xk[j] << endl;
//        }

    double sum = 0.0;
    for (int j=0;j<n;j++)
        {
            double part1 = sqrt(1.0-xk[j]*xk[j]);
            sum = sum + Ak[j]*part1*(*func)(xk[j]);
        }

    cout << "利用切比雪夫 积分e^x在区间-1到1之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}

//这里是变换的版本
void qiebi_inter_fp1(double (*func)(double x), double lower_limit, double upper_limit)
{
    const int n = 20;
    double Ak[n]={0.0};
    double xk[n]={0.0};
    for (int j=0;j<n;j++)
        {
            Ak[j]=pi/(n+1);
            xk[j]=cos((2.0*j+1.0)/(2.0*n+2.0)*pi);//注意整除问题 要带个1.0 保证结果是double
        }

//        for (int j=0;j<n;j++)
//        {
//           cout << Ak[j] << '\t' << xk[j] << endl;
//        }

    double sum = 0.0;
    for (int j=0;j<n;j++)
        {
            double part1 = 0.5*(upper_limit-lower_limit)*xk[j]+0.5*(upper_limit+lower_limit);
            double part2 = sqrt(1.0-part1*part1);
            sum = sum + Ak[j]*part2*(*func)(part1);
        }
    sum = sum*0.5*(upper_limit-lower_limit);
    cout << "利用切比雪夫 积分e^x在区间-1到1之间的值为：" << sum << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    return;
}


int main()
{

    juxing_inter();
    juxing_inter_fp(func,0.0,pi);
    tixing_inter();
    tixing_inter_fp(func,0.0,pi);
    xps_dbc_inter();
    xps_dbc_inter_fp(func2,0.0,1.0);
    double ret = xps_bbc_integral(0.0,1.0,0.000001);
    cout << "利用变步长辛普森法 积分log(1+x)/(1+x^2)在区间0到1之间的值为：" << setprecision(8) << ret << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    tixing_inter_2d();
    double ret1 = xps_bbc_integral2d(0.0,1.0,0.00005);
    cout << "利用变步长辛普森法 积分二重e^(x^2+y^2)在区间x从0到1 y从-sqrt(1-x^2)到sqrt(1-x^2)之间的值为：" << setprecision(8) << ret1 << endl;
    cout << "-------------------------------------------------------------------"<<endl;
    gauss_inter_fp(func2,0.0,1.0);
    gauss_lerande_inter_fp(func3,2.5,8.4);
    qiebi_inter_fp(func4, -1.0, 1.0);
    qiebi_inter_fp1(func4, -1.0, 1.0);
    return 0;
}





