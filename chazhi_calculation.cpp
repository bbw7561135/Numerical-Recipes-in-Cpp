//!===========================================================================================!
//!===========================================================================================!
//! Numerical Algorithms Collections In Fortran 77 Written by Liu Fei Yu & Li                 !
//! Translate into C++ for easy access by Bao Biwen                                           !
//! Created on 2020.06.26 with Version 1.0                                                    !
//! Chapter 2 Interpolation and numerical differentials                                         !
//!===========================================================================================!
//!===========================================================================================!

#include <iostream>
#include <cmath>
using namespace std;

const double pi = 3.14159265358;

//一元全区间不等距插值
//对给定的n个插值节点x1-xn及对应的函数值y1-yn 利用拉格朗日插值公式 计算x点的函数值y(x)
//Large_chazhi where y_want should be pass in pointer as it will change the value
//the formula is in page32
void Large_chazhi(const double *x_arr, const double *y_arr, const int arr_len, double x_want, double* y_want)
{

    double sum=0.0;
    for(int i=0; i<arr_len; i++)
    {
        double lian_cheng = 1.0; //multiply
        for(int j=0; j<arr_len; j++)
        {
            if(j==i)
            {
                continue;
            }
            else
            {
                lian_cheng = lian_cheng * (x_want - x_arr[j])/(x_arr[i] - x_arr[j]);
            }
        }
        sum = sum + y_arr[i]*lian_cheng;
    }
    *y_want=sum;
    return;
}


//一元三点等距插值和微分 等距=dx相等 三点=参与计算的仅为相邻三点 既然仅用三点就存在选取点的判断方法
//利用拉格朗日三点插值公式对一元等距插值表选取相邻三点进行插值或求微分
void lg_czwf_3dian_dengju(const double *x_arr, const double *y_arr, const int arr_len,int intend, double x_want,double* y_want)
{
    //intend=0 cha zhi intend=1 wei fen
    double h = (x_arr[arr_len-1]-x_arr[0])/(arr_len-1.0);
    double i_1 = (x_want-x_arr[0])/h+0.5;
    if(i_1<1.0)
    {
        i_1=1.0; //in fortran the first element index is 1 while in c/c++ is 0
    }
    else if(i_1 <= arr_len-2.0)
    {
        i_1=i_1;
    }
    else
    {
        i_1= arr_len - 2.0;//fortran is from 1 to N while in c it is 0 to N-1 for N points N-1 intervals
    }

    double t = (x_want-x_arr[0])/h-(int)(i_1); //should be careful i_1 is an int
    if(intend==0) //chazhi
    {
        int ii = (int)(i_1);
        double c1 = 0.5*(t*t-t);
        double c2 = 1.0-t*t;
        double c3 = 0.5*(t*t+t);
        *y_want = c1*y_arr[ii-1]+c2*y_arr[ii]+c3*y_arr[ii+1];
        //fortran is from 1 to N while in c it is 0 to N-1 for N points N-1 intervals
    }
    if(intend==1) //weifen
    {
        int ii = (int)(i_1);
        double c1 = 0.5*(2.0*t-1.0);
        double c2 = -2.0*t;
        double c3 = 0.5*(2.0*t+1.0);
        *y_want = (c1*y_arr[ii-1]+c2*y_arr[ii]+c3*y_arr[ii+1])/h;
    }
    return;
}

//一元三点不等距插值
//选取最靠近插值点的相邻的三个插值点 用拉格朗日三点插值公式进行插值
void lg_cz_3dian_budengju(const double *x_arr, const double *y_arr, const int arr_len, double x_want,double* y_want)
{
    int i=0;
    if(x_want >= 0.5*(x_arr[arr_len-3]+x_arr[arr_len-2]))
    {        //fortran is from 1 to N while in c it is 0 to N-1 for N points N-1 intervals
        i = arr_len - 3;
    }
    else
    {
        for(int j=1;j<arr_len;j++)
        {
            if(x_want < 0.5*(x_arr[j]+x_arr[j+1]))
            {
                i = j-1;
                break;
            }
        }
    }

    double sum=0.0;
    for(int k=i;k<=i+2;k++)
    {
        double lian_cheng = 1.0;
        for(int j=i; j<=i+2; j++)
        {
            if(j==k)
            {
                continue;
            }
            else
            {
                lian_cheng = lian_cheng * (x_want-x_arr[j])/(x_arr[k]-x_arr[j]);
            }
        }
        sum = sum + y_arr[k]*lian_cheng;
    }
    *y_want = sum;
    return;
}

//一元三点不等距成组插值
//选取最靠近插值点的相邻的三个插值点 用拉格朗日三点插值公式进行插值
//和一元三点不等距插值不同的在于判断i 也就是第一个插值点有所不同 其他皆相同
//所谓成组 就是一个x数组 对应多个y数组 其实没必要 分开算也行 因为y组之间是独立的 只是x公用罢了
void lg_cz_3dian_czu_budengju(const double *x_arr, const double *y_arr, const int arr_len, double x_want,double* y_want)
{
    int i=0;
    if(x_want > x_arr[arr_len-2])
    {        //fortran is from 1 to N while in c it is 0 to N-1 for N points N-1 intervals
        i = arr_len - 3;
    }
    else if(x_want < x_arr[1])
    {
        i = 0;
    }
    else
    {
        for(int j=1;j<arr_len;j++)
        {
            if(x_want > x_arr[j] && x_want < x_arr[j+1])
            {
                if(fabs(x_want-x_arr[j])>fabs(x_arr[j+1]-x_want))
                {
                    i = j;
                    break;
                }
                else
                {
                    i=j-1;
                    break;
                }
            }
        }
    }

    double sum=0.0;
    for(int k=i;k<=i+2;k++)
    {
        double lian_cheng = 1.0;
        for(int j=i; j<=i+2; j++)
        {
            if(j==k)
            {
                continue;
            }
            else
            {
                lian_cheng = lian_cheng * (x_want-x_arr[j])/(x_arr[k]-x_arr[j]);
            }
        }
        sum = sum + y_arr[k]*lian_cheng;
    }
    *y_want = sum;
    return;
}

//埃特金插值要反复用到两点线性插值
double linear_2pcz(double x0, double y0, double x1, double y1, double xwant)
{
    double ywant=0.0;
    ywant = (xwant-x1)/(x0-x1)*y0+(xwant-x0)/(x1-x0)*y1;
    return ywant;
}

//埃特金插值
//从给定的n个插值节点中选取最靠近插值点x的相邻的m(m<=n)个插值点 尽量使x位于m个插值点中心
//用埃特金反复线性插值公式对一元函数进行插值 反复线性插值逐步产生一次 二次 至m-1次拉格朗日多项式
//比如m=3 可以产生2次拉格朗日插值多项式 这里默认m=3
//原书中的i是选取的m点的起始插值点 i是在n个插值点中的位置
void atj_cz(const double *x_arr, const double *y_arr, const int arr_len,double x_want,double* y_want)
{
    int i = 0;
    for(int j=0;j<arr_len-1;j++)
    {
        if(x_want < x_arr[0])
        {
            i = 0;
            break;
        }
        if(x_want > x_arr[arr_len-1])
        {
            i = arr_len-3;
            break;
        }
        if(x_want>x_arr[j] && x_want<x_arr[j+1])
        {
            if(arr_len-1-j==1)//要是j是倒数第二个插值点 那么只有两个插值点不够 要往前一个点
            {
                i=j-1;
                break;
            }
            else
            {
                i=j;
                break;
            }
        }
    }
    //m=3 3个插值点 x0y0 x1y1 x2y2 x为待插值点
    //利用x0y0和x1y1构造一个一次多项式 并计算x的值 得P01(x)
    //利用x0y0和x2y2构造一个一次多项式 并计算x的值 得P02(x)
    //利用x1P01(x)和x2P02(x)构造一个一次多项式 计算P012(x)即为所求值
    double P01x = linear_2pcz(x_arr[i],y_arr[i],x_arr[i+1],y_arr[i+1],x_want);
    double P02x = linear_2pcz(x_arr[i],y_arr[i],x_arr[i+2],y_arr[i+2],x_want);
    double P012x = linear_2pcz(x_arr[i+1],P01x,x_arr[i+2],P02x,x_want);
    *y_want = P012x;
    return;
}

//埃尔密特插值
//给定n个插值点x1-xn及对应的函数值y1-yn和一阶导数y1'-yn'
// 利用埃尔密特公式对一元函数进行插值
void amt_cz(const double *x_arr, const double *y_arr, const double *dy_arr,const int arr_len,double x_want,double* y_want)
{
    double sum=0.0;
    for(int i=0; i<arr_len; i++)
    {
        double lian_cheng_hi = 1.0; //累乘因子hi
        double lian_jia_ai=0.0;//累加因子ai
        for(int j=0; j<arr_len; j++)
        {
            if(j==i)
            {
                continue;
            }
            else
            {
                lian_cheng_hi= lian_cheng_hi * (x_want - x_arr[j])/(x_arr[i] - x_arr[j])* (x_want - x_arr[j])/(x_arr[i] - x_arr[j]);
                lian_jia_ai = lian_jia_ai + 1.0/(x_arr[i] - x_arr[j]);
            }
        }
        sum = sum + lian_cheng_hi*(y_arr[i]+(x_arr[i]-x_want)*(2.0*lian_jia_ai*y_arr[i]-dy_arr[i]));
    }
    *y_want=sum;
    return;
}


int main()
{
    const int n = 6;
    double x[n] = {0.0,0.1,0.195,0.3,0.401,0.5};
    double y[n] = {0.39894,0.39695,0.39142,0.38138,0.36812,0.35206};
    double x_want = 0.15;
    double x_want1 = 0.3;
    double x_want2 = 0.45;
    double y_want = 0.0;
    //一元全区间不等距插值 一元=仅x为因变量 全区间=所有x都参与插值计算 不等距=dx在区间内有不同 插值=利用周围点计算x处的函数值y(x)
    //对给定的n个插值节点x1-xn及对应的函数值y1-yn 利用拉格朗日插值公式 计算x点的函数值y(x)
    //Large_chazhi where y_want should be pass in pointer as it will change the value
    //the formula is in page32
    cout<< "--------------------一元全区间不等距插值--------------------------" << endl;
    Large_chazhi(x,y,n,x_want,&y_want);
    cout << "x(0.15) = "<<y_want << endl;
    cout<< "-----------------------------------------------" << endl;
    Large_chazhi(x,y,n,x_want1,&y_want);
    cout << "x(0.3) = "<<y_want << endl;
    cout<< "-----------------------------------------------" << endl;
    Large_chazhi(x,y,n,x_want2,&y_want);
    cout << "x(0.45) = "<<y_want << endl;
    cout<< "-----------------------------------------------" << endl;

    int intend = 0;//0 for chazhi 1 for wei fen
    double x1[n] = {0.0,0.1,0.2,0.3,0.4,0.5};
    double y1[n] = {0.39894,0.39695,0.39104,0.38138,0.36827,0.35206};
    double x_want3 = 0.04;
    double y_want3 = 0.0;
    //一元三点等距插值和微分 等距=dx相等 三点=参与计算的仅为相邻三点 既然仅用三点就存在选取点的判断方法
    //利用拉格朗日三点插值公式对一元等距插值表选取相邻三点进行插值或求微分
    cout<< "--------------------一元三点等距插值和微分------------------------" << endl;
    lg_czwf_3dian_dengju(x1,y1,n,intend,x_want3,&y_want3);
    cout <<"x(0.04) =  "<<y_want3 << endl;
    cout<< "-----------------------------------------------" << endl;
    double x_want4 = 0.2;
    lg_czwf_3dian_dengju(x1,y1,n,intend,x_want4,&y_want3);
    cout <<"x(0.2) =  "<<y_want3 << endl;
    cout<< "-----------------------------------------------" << endl;
    double x_want5 = 0.24;
    lg_czwf_3dian_dengju(x1,y1,n,intend,x_want5,&y_want3);
    cout <<"x(0.24) =  "<<y_want3 << endl;
    cout<< "-----------------------------------------------" << endl;
    double x_want55 = 0.48;
    lg_czwf_3dian_dengju(x1,y1,n,intend,x_want55,&y_want3);
    cout <<"x(0.48) =  "<<y_want3 << endl;
    cout<< "-----------------------------------------------" << endl;

    double x_want6 = 0.0;
    double y_want6 = 0.0;
    //一元三点不等距插值
    //选取最靠近插值点的相邻的三个插值点 用拉格朗日三点插值公式进行插值
    cout<< "---------------------一元三点不等距插值---------------------------" << endl;
    lg_cz_3dian_budengju(x,y,n,x_want6,&y_want6);
    cout <<"x(0.0) = "<< y_want6 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_cz_3dian_budengju(x,y,n,x_want,&y_want6);
    cout << "x(0.15) = "<< y_want6 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_cz_3dian_budengju(x,y,n,x_want1,&y_want6);
    cout << "x(0.3) = "<< y_want6 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_cz_3dian_budengju(x,y,n,x_want2,&y_want6);
    cout << "x(0.45) = "<< y_want6 << endl;
    cout<< "-----------------------------------------------" << endl;



    double x2[10] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};
    double y2[10] = {0.0,0.3096,0.5712,0.7546,0.8561,0.8954,0.9028,0.9036,0.9110,0.9248};
    double x_want7 = 3.25;
    double y_want7 = 0.0;
    intend = 0;//0 for chazhi 1 for wei fen
    //一元三点分段等距插值和微分
    cout<< "---------------------一元三点分段等距插值和微分---------------------------" << endl;
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7,&y_want7);
    cout <<"x(3.25) = "<< y_want7 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7+1.0,&y_want7);
    cout <<"x(4.25) = "<< y_want7 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7+1.5,&y_want7);
    cout <<"x(4.75) = "<< y_want7 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7+1.75,&y_want7);
    cout <<"x(5.0) = "<< y_want7 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7+5.0,&y_want7);
    cout <<"x(8.25) = "<< y_want7 << endl;
    cout<< "-----------------------------------------------" << endl;
    double x3[4] = {0.0,1.0,2.0,3.0};
    double y3[4] = {0.0,0.3096,0.5712,0.7546};
    lg_czwf_3dian_dengju(x3,y3,4,intend,3.25,&y_want7);
    cout <<"x(3.25) = "<< y_want7 << endl;
    cout<< "-----------------------------------------------" << endl;

    double x4[7] = {3.0,4.0,5.0,6.0,7.0,8.0,9.0};
    double y4[7] = {0.7546,0.8561,0.8954,0.9028,0.9036,0.9110,0.9248};
    lg_czwf_3dian_dengju(x4,y4,7,intend,3.25,&y_want7);
    cout <<"x(3.25) = "<< y_want7 << endl;
    cout<< "-----------------------------------------------" << endl;
    //the book uses two batches of chazhi dian
    //for x < 3.0 use x3 y3 for x > 3.0 use x4 y4

    double x5[10] = {0.0,0.8,1.2,2.0,2.6,3.4,4.2,5.0,5.6,5.8};
    double y5[10] = {0.0,0.2502,0.3671,0.5712,0.6915,0.8043,0.8681,0.8954,0.9019,0.9025};
    double x_want8 = 5.5;
    double y_want8 = 0.0;
    //一元三点不等距分段插值
    cout<< "---------------------一元三点不等距分段插值---------------------------" << endl;
    lg_cz_3dian_budengju(x5,y5,10,x_want8,&y_want8);
    cout <<"x(5.5) = "<< y_want8 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_cz_3dian_budengju(x5,y5,10,x_want8-2.5,&y_want8);
    cout <<"x(3.0) = "<< y_want8 << endl;
    cout<< "-----------------------------------------------" << endl;

    cout<< "---------------------一元三点不等距成组插值---------------------------" << endl;
    //一元三点不等距成组插值
    //选取最靠近插值点的相邻的三个插值点 用拉格朗日三点插值公式进行插值
    //和一元三点不等距插值不同的在于判断i 也就是第一个插值点有所不同 其他皆相同
    double x6[6]={0.20,0.24,0.28,0.32,0.36,0.40};
    double y6a[6]={0.19867,0.23770,0.27636,0.31457,0.35227,0.38942};
    double y6b[6]={0.80400,0.76576,0.72784,0.69024,0.65296,0.61600};
    double x_want9 = 0.22;
    double y_want9 = 0.0;
    lg_cz_3dian_czu_budengju(x6,y6a,6,x_want9,&y_want9);
    cout <<"x(0.22) = "<< y_want9 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_cz_3dian_czu_budengju(x6,y6a,6,x_want9+0.16,&y_want9);
    cout <<"x(0.38) = "<< y_want9 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_cz_3dian_czu_budengju(x6,y6b,6,x_want9,&y_want9);
    cout <<"x(0.22) = "<< y_want9 << endl;
    cout<< "-----------------------------------------------" << endl;
    lg_cz_3dian_czu_budengju(x6,y6b,6,x_want9+0.16,&y_want9);
    cout <<"x(0.38) = "<< y_want9 << endl;
    cout<< "-----------------------------------------------" << endl;

    double x7[4]={0.5,0.65,0.8,1.0};
    double y7[4]={0.4794,0.6052,0.7174,0.8415};
    //埃特金插值
    //从给定的n个插值节点中选取最靠近插值点x的相邻的m(m<=n)个插值点 尽量使x位于m个插值点中心
    //用埃特金反复线性插值公式对一元函数进行插值 反复线性插值逐步产生一次 二次 至m-1次拉格朗日多项式
    //比如m=3 可以产生2次拉格朗日插值多项式 这里默认m=3
    cout<< "----------------------------埃特金插值--------------------------------" << endl;
    double x_want10 = 0.1;
    double y_want10 = 0.0;
    atj_cz(x7,y7,4,x_want10, &y_want10);
    cout <<"x(0.1) = "<< y_want10 << endl;
    cout<< "-----------------------------------------------" << endl;
    atj_cz(x7,y7,4,x_want10+0.4, &y_want10);
    cout <<"x(0.5) = "<< y_want10 << endl;
    cout<< "-----------------------------------------------" << endl;
    atj_cz(x7,y7,4,x_want10+0.5, &y_want10);
    cout <<"x(0.6) = "<< y_want10 << endl;
    cout<< "-----------------------------------------------" << endl;
    atj_cz(x7,y7,4,x_want10+1.2, &y_want10);
    cout <<"x(1.3) = "<< y_want10 << endl;
    cout<< "-----------------------------------------------" << endl;

    cout<< "----------------------------埃尔密特插值--------------------------------" << endl;
    double x_want11 = 0.25;
    double y_want11 = 0.0;
    double x8[3]={0.1,0.3,0.5};
    double y8[3]={0.099833,0.295520,0.479426};
    double dydx8[3]={0.995004,0.955336,0.877583};
    amt_cz(x8,y8,dydx8,3,x_want11,&y_want11);
    cout <<"x(0.25) = "<< y_want11<< endl;
    cout<< "-----------------------------------------------" << endl;

    return 0;
}





