//!===========================================================================================!
//!===========================================================================================!
//! Numerical Algorithms Collections In Fortran 77 Written by Liu Fei Yu & Li                 !
//! Translate into C++ for easy access by Bao Biwen                                           !
//! Created on 2020.06.26 with Version 1.0                                                    !
//! Chapter 2 Interpolation and numerical differeting                                         !
//!===========================================================================================!
//!===========================================================================================!

#include <iostream>
#include <cmath>
using namespace std;

const double pi = 3.14159265358;

//quan qu jian bu deng ju cha zhi
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
}


//san dian deng ju cha zhi he wei fen
void lg_czwf_3dian_dengju(const double *x_arr, const double *y_arr, const int arr_len,int intend, double x_want,double* y_want)
{
    //intend=0 cha zhi intend=1 wei fen
    double h = (x_arr[arr_len-1]-x_arr[0])/(arr_len-1);
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
    }
    if(intend==1) //weifen
    {
        int ii = (int)(i_1);
        double c1 = 0.5*(2.0*t-1.0);
        double c2 = -2.0*t;
        double c3 = 0.5*(2.0*t+1.0);
        *y_want = (c1*y_arr[ii-1]+c2*y_arr[ii]+c3*y_arr[ii+1])/h;
    }
}

//san dian bu deng ju cha zhi
void lg_cz_3dian_budengju(const double *x_arr, const double *y_arr, const int arr_len, double x_want,double* y_want)
{
    int i=0;
    if(x_want >= 0.5*(x_arr[arr_len-3]+x_arr[arr_len-2]))
    {
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
    Large_chazhi(x,y,n,x_want,&y_want);
    cout << y_want << endl;
    Large_chazhi(x,y,n,x_want1,&y_want);
    cout << y_want << endl;
    Large_chazhi(x,y,n,x_want2,&y_want);
    cout << y_want << endl;
    cout<< "-------------------" << endl;
    cout<< "-------------------" << endl;
    int intend = 0;//0 for chazhi 1 for wei fen
    double x1[n] = {0.0,0.1,0.2,0.3,0.4,0.5};
    double y1[n] = {0.39894,0.39695,0.39104,0.38138,0.36827,0.35206};
    double x_want3 = 0.04;
    double y_want3 = 0.0;
    lg_czwf_3dian_dengju(x1,y1,n,intend,x_want3,&y_want3);
    cout << y_want3 << endl;
    double x_want4 = 0.2;
    lg_czwf_3dian_dengju(x1,y1,n,intend,x_want4,&y_want3);
    cout << y_want3 << endl;
    double x_want5 = 0.24;
    lg_czwf_3dian_dengju(x1,y1,n,intend,x_want5,&y_want3);
    cout << y_want3 << endl;
    cout<< "-------------------" << endl;
    cout<< "-------------------" << endl;
    double x_want6 = 0.0;
    double y_want6 = 0.0;
    lg_cz_3dian_budengju(x,y,n,x_want6,&y_want6);
    cout << y_want6 << endl;
    lg_cz_3dian_budengju(x,y,n,x_want,&y_want6);
    cout << y_want6 << endl;
    lg_cz_3dian_budengju(x,y,n,x_want1,&y_want6);
    cout << y_want6 << endl;
    lg_cz_3dian_budengju(x,y,n,x_want2,&y_want6);
    cout << y_want6 << endl;

    cout<< "-------------------" << endl;
    cout<< "-------------------" << endl;

    double x2[10] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};
    double y2[10] = {0.0,0.3096,0.5712,0.7546,0.8561,0.8954,0.9028,0.9036,0.9110,0.9248};
    double x_want7 = 3.25;
    double y_want7 = 0.0;
    intend = 0;//0 for chazhi 1 for wei fen
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7,&y_want7);
    cout << y_want7 << endl;
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7+1.0,&y_want7);
    cout << y_want7 << endl;
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7+1.5,&y_want7);
    cout << y_want7 << endl;
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7+1.75,&y_want7);
    cout << y_want7 << endl;
    lg_czwf_3dian_dengju(x2,y2,10,intend,x_want7+5.0,&y_want7);
    cout << y_want7 << endl;
    cout<< "-------------------" << endl;
    cout<< "-------------------" << endl;
    double x3[4] = {0.0,1.0,2.0,3.0};
    double y3[4] = {0.0,0.3096,0.5712,0.7546};
    lg_czwf_3dian_dengju(x3,y3,4,intend,3.25,&y_want7);
    cout << y_want7 << endl;

    double x4[7] = {3.0,4.0,5.0,6.0,7.0,8.0,9.0};
    double y4[7] = {0.7546,0.8561,0.8954,0.9028,0.9036,0.9110,0.9248};
    lg_czwf_3dian_dengju(x4,y4,7,intend,3.25,&y_want7);
    cout << y_want7 << endl; //the book uses two batches of chazhi dian
    //for x < 3.0 use x3 y3 for x > 3.0 use x4 y4
    
    double x5[10] = {0.0,0.8,1.2,2.0,2.6,3.4,4.2,5.0,5.6,5.8};
    double y5[10] = {0.0,0.2502,0.3671,0.5712,0.6915,0.8043,0.8681,0.8954,0.9019,0.9025};

    double x_want8 = 5.5;
    double y_want8 = 0.0;
    lg_cz_3dian_budengju(x5,y5,10,x_want8,&y_want8);
    cout << y_want8 << endl;
    return 0;
}





