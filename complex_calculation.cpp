//!===========================================================================================!
//!===========================================================================================!
//! Numerical Algorithms Collections In Fortran 77 Written by Liu Fei Yu & Li                 !
//! Translate into C++ for easy access by Bao Biwen                                           !
//! Created on 2020.06.19 with Version 1.0                                                    !
//! Chapter 1 complex calculation                                                             !
//!===========================================================================================!
//!===========================================================================================!
#include <iostream>
#include <cmath>
using namespace std;

const double pi = 3.14159265358;

//(a+ib)/(c+id) the result is e+if
//(a+ib)/(c+id)=[(ac+bd)+i(bc-cd)]/(c^2+d^2)
void complex_divide(double a, double b, double c, double d, double* e, double* f)
{
    if(abs(c)>=abs(d))
    {
        *e=(a+d/c*b)/(c+d/c*d);
        *f=(b-d/c*a)/(c+d/c*d);
    }
    else
    {
        *e=(c/d*a+b)/(c/d*c+d);
        *f=(c/d*b-a)/(c/d*c+d);
    }

    return;
}

//e^z where z=a+ib the result is c+df
//e^z=e^a cos(b)+ie^a sin(b)
void complex_ez(double a, double b, double* c, double* d)
{
    *c=exp(a)*cos(b);
    *d=exp(a)*sin(b);
    return;
}


//z=a+ib z!=0
//ln(z) = ln|z|e^(i argz) = ln|z|+i*argz
//argz=augment(z) -pi<argz<=pi
//ln(z)=c+id where c=ln|z| d=argz
void complex_lnz(double a, double b, double* c, double* d)
{
    //ln|z|=ln(sqrt(a^2+b^2))=0.5*ln(a^2+b^2)
    //the following procedure is to avoid overflow as in early years the computer is not so powerful
    if(abs(a)<1.0 && abs(b)<1.0)
    {
        double part1 = log(2.0*abs(a)+2.0*abs(b));
        double part2 = log(8.0*a*a/(2.0*abs(a)+2.0*abs(b))+8.0*b*b/(2.0*abs(a)+2.0*abs(b)));
        *c = 0.5*(part1+part2)-0.5*log(8.0);
    }
    else if(abs(a)>=1.0 && abs(b)>=1.0)
    {
        double part1 = log(0.25*abs(a)+0.25*abs(b));
        double part2 = log(0.5*0.25*a*a/(0.25*abs(a)+0.25*abs(b))+0.5*0.25*b*b/(0.25*abs(a)+0.25*abs(b)));
        *c = 0.5*(part1+part2)+0.5*log(8.0);
    }
    else
    {
        *c = 0.5*log(a*a+b*b);
    }

    if(abs(a)>=abs(b) && a != 0.0)
    {
        if(a>0.0)
        {
            *d=atan(b/a);
        }
        else if(b>=0.0)
        {
            *d=atan(b/a)+pi;
        }
        else
        {
            *d=atan(b/a)-pi;
        }
    }
    else
    {
        double signb = (b>0.0) ? 1.0:-1.0;
        *d=-atan(a/b)+ signb*0.5*pi;
    }
        return;
}

//z=a+ib |z|=sqrt(a^2+b^2)
void complex_mode(double a, double b, double* c)
{
    if(abs(a)>abs(b))
    {
        *c = abs(a)*sqrt(1.0+pow(b/a,2.0));
    }
    else
    {
        *c = abs(b)*sqrt(1.0+pow(a/b,2.0));
    }
    return;
}


//z=x+iy sqrt(z)=a+ib (where a>0 is wanted)
void complex_sqrt(double a, double b, double *c, double *d)
{
    double mode = 0.0;
    if(a>=0.0)
    {
        complex_mode(a,b,&mode);
        *c = sqrt(0.5*a+0.5*mode);
        *d = 0.5*b/(*c);
    }
    else
    {
        complex_mode(a,b,&mode);
        *d = sqrt(0.5*a+0.5*mode);
        *c = 0.5*b/(*d);
    }

    return;
}


//sin(z) where z = x + iy the result in c + id
//sinh and cosh is in cmath
void complex_sinz(double a, double b, double *c, double *d)
{
    *c = sin(a)*cosh(b);
    *d = cos(a)*sinh(b);
    return;
}

//cos(z) where z = x + iy the result is in c + id
//sinh and cosh is in cmath
void complex_cosz(double a, double b, double *c, double *d)
{
    *c = cos(a)*cosh(b);
    *d = -sin(a)*sinh(b);
    return;
}


//tan(z) where z = x + iy the result is in c + id
void complex_tanz(double a, double b, double *c, double *d)
{
    double sin_re = 0.0;
    double sin_im = 0.0;
    double cos_re = 0.0;
    double cos_im = 0.0;
    complex_sinz(a,b,&sin_re,&sin_im);
    complex_cosz(a,b,&cos_re,&cos_im);
    complex_divide(sin_re,sin_im,cos_re,cos_im,c,d);
    return;
}

//u^t where u=a+ib t=c+id u!=0
//u^t=e^tln(u)=e^t(ln(u)+i2npi) when n=0 it is the principal value //i2npi is due to augment(z)= argz + 2npi
//e^t(ln(u)+i2npi)=e^(c+id)(ln(u)+i2npi) where ln(u)=ln|u|+iarg(u)
//then u^t=e^(cln|u|-d(arg(u)+2npi)) *(cos(V)+isin(V))
//where V=d*ln|u|+c(arg(u)+2npi)
//the final result is e+if
void complex_zz(double a, double b, double c, double d, double* e, double *f)
{
    double lnu = 0.0; //mode ln|u|
    double argu = 0.0;
    double n = 0.0; //n=0 to get the principal value
    complex_lnz(a,b,&lnu,&argu);
    double part1 = c*lnu-d*(argu+2.0*n*pi);
    double part2 = exp(part1);
    double part3 = d*lnu+c*(argu+2.0*n*pi);
    *e = part2*cos(part3);
    *f = part2*sin(part3);
    return;
}

int main()
{
    double a = 1.5e10;
    double b = 1.0e20;
    double c = 2.0e38;
    double d = 1.0e30;
    double e = 0.0;
    double f = 0.0;
    complex_divide(a,b,c,d,&e,&f);
    //(a+ib)/(c+id) the result is e+if
    cout << "(a+ib)/(c+id) the result is e+if" <<endl;
    cout << e << " + " << f << "i" << endl;
    cout << "----------------------------------------" << endl;

    double aa = 1.0;
    double bb = pi/4.0; //in rad unit
    double cc = 0.0;
    double dd = 0.0; //cos、sin、asin、acos in cmath is in rad unit
    complex_ez(aa,bb,&cc,&dd);
    //e^z where z=a+ib the result is c+df
    cout << "e^z where z=a+ib the result is c+df" << endl;
    cout << cc << " + " << dd << "i" << endl;
    cout << "----------------------------------------" << endl;

    double aaa = 1.922115512;
    double bbb = 1.922115512;
    double ccc = 0.0;
    double ddd = 0.0;
    //z=a+ib z!=0
    //ln(z) = ln|z|e^(i argz) = ln|z|+i*argz
    //argz=augment(z) -pi<argz<=pi
    //ln(z)=c+id where c=ln|z| d=argz
    complex_lnz(aaa,bbb,&ccc,&ddd);
    cout << "ln(z) where z=a+ib and z !=0" << endl;
    cout << ccc << " + " << ddd << "i" << endl;
    cout << "----------------------------------------" << endl;

    double a1 = 1.264e38;
    double b1 = 1.548e38;
    double c1 = 0.0;
    //z=a+ib |z|=sqrt(a^2+b^2)
    complex_mode(a1,b1,&c1);
    cout << "|z|=sqrt(a^2+b^2) where z=a+ib" << endl;
    cout << c1 << endl;
    cout << "----------------------------------------" << endl;

    double a2 = 1.264e38;
    double b2 = 1.548e38;
    double c2 = 0.0;
    double d2 = 0.0;
    //z=x+iy sqrt(z)=a+ib
    complex_sqrt(a2,b2,&c2,&d2);
    cout << "z=x+iy sqrt(z)=a+ib" << endl;
    cout << c2 << " + " << d2 << "i" << endl;
    cout << "----------------------------------------" << endl;

    double a3 = 0.25;
    double b3 = 0.25;
    double c3 = 0.0;
    double d3 = 0.0;
    double c4 = 0.0;
    double d4 = 0.0;
    //sin(z) where z = x + iy the result in c + id
    //sinh and cosh is in cmath
    complex_sinz(a3,b3,&c3,&d3);
    cout << "sin(z) where z = x + iy the result in c + id" << endl;
    cout << c3 << " + " << d3 << "i" << endl;
    cout << "----------------------------------------" << endl;
    //cos(z) where z = x + iy the result in c + id
    //sinh and cosh is in cmath
    cout << "cos(z) where z = x + iy the result in c + id" << endl;
    complex_cosz(a3,b3,&c4,&d4);
    cout << c4 << d4 << "i" << endl;
    cout << "----------------------------------------" << endl;
    //tan(z) where z = x + iy the result is in c + id
    //tan(z)=sin(z)/cos(z)
    double c5 = 0.0;
    double d5 = 0.0;
    cout << "tan(z) where z = x + iy the result is in c + id" << endl;
    complex_tanz(a3,b3,&c5,&d5);
    cout << c5 << "+" << d5 << "i" << endl;
    cout << "----------------------------------------" << endl;

    //sinh and cosh is in cmath
    //cosh(y)=(e^y+e^-y)/2
    //sinh(y)=(e^y-e^-y)/2
    cout << "sinh(y) and cosh(y) where y is real number" << endl;
    cout << "sinh(0.25)=" << sinh(b3) << endl;
    cout << "cosh(0.25)=" << cosh(b3) << endl;
    cout << "----------------------------------------" << endl;

    double aa1 = 1.0;
    double bb1 = 1.0;
    double cc1 = 1.0;
    double dd1 = 1.0;
    double ee1 = 0.0;
    double ff1 = 0.0;
    complex_zz(aa1,bb1,cc1,dd1,&ee1,&ff1);
    //u^t where u=a+ib t=c+id u!=0
    cout << "u^t where u=a+ib t=c+id u!=0" << endl;
    cout << ee1 << "+" << ff1 << "i" << endl;
    cout << "----------------------------------------" << endl;
    return 0;
}
