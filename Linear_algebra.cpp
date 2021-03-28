//!===========================================================================================!
//!===========================================================================================!
//! Numerical Algorithms Collections In Fortran 77 Written by Liu Fei Yu & Li                 !
//! Translate into C++ for easy access by Bao Biwen                                           !
//! Created on 2021.03.25 with Version 1.0                                                    !
//! Chapter 4  Linear algebra                                                                 !
//!===========================================================================================!
//!===========================================================================================!

#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const double pi = 3.14159265358;


void Gauss_calculation(void)//Gauss消去法解线性方程组
{
    //给定矩阵
    const int n=3;
    double x[n][n]={1.0,-1.0,1.0,5.0,-4.0,3.0,2.0,1.0,1.0};
    double b[n]={-4.0,-12.0,11.0};
    double xans[n]={0.0};
    cout << "原始增广矩阵" <<endl;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
         cout << x[i][j] << '\t';
        }
        cout << "|||" << b[i] << endl;
    }
    cout<<"利用以上A与b组成的增广阵进行高斯消去法计算方程组"<<endl;
    double bili = 0.0; //比例系数
    for (int i=1;i<n;i++)
        {
            for(int j=i;j<n;j++)
            {
                bili = x[j][i-1]/x[i-1][i-1];//确定本次消元的比例某行的比例系数
                b[j]=b[j]-bili*b[i-1];//先处理b
                for (int k=i-1;k<n;k++)
                    {
                        x[j][k]=x[j][k]-bili*x[i-1][k];
                    }

            }
        }
    //高斯消元法存在失真问题 比如答案是0 可能出现1.23e-16这样的答案
    //还有就是消元的过程中 可能出现对角线元素为0 即使不为0 如果x_kk很小 也会影响数值稳定性
    //这种程序还是直接去参考别人编好的程序比较稳妥
    //其他衍生方法 列主元 全主元 列主元-约当消去法 等其他高级方法
    cout <<"顺序消元后的上三角系数增广矩阵如下" <<endl;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
         cout << x[i][j] << '\t';
        }
        cout << "|||" << b[i] << endl;
    }
    cout << "利用回代法求解上三角方程组得到解" << endl;
    for (int j=n-1;j>=0;j--)//这里注意j>=0 而不是j<=0;
        {
            for (int k=0;k<n;k++)
                {
                    if(k==j)
                    {
                        continue;
                    }
                    b[j] = b[j] - xans[k]*x[j][k];
                }
            xans[j]=b[j]/x[j][j];
        }

    for(int i=0;i<n;i++)
    {
         cout << xans[i] << endl;
    }

}

int main()
{
    Gauss_calculation();
    return 0;
}
