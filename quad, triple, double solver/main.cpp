#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cstdarg>
#include <complex>
using namespace std;

void complex_sqrt(double &a, double &b)
{
	double re = a, im = b, r  = sqrt(re*re+im*im);
    if(im == 0) 
    { 
        r = sqrt(r);
        if(re >= 0) {a=r; b=0;} 
        else {a=0; b=r;}
    } 
    else 
    {// im != 0
        a = sqrt(0.5*(re+r));
        b = 0.5*im/a;
    }
}

bool too_small(double a, double e)
{	
	return fabs(a) < e;
}
int sol2(double *x, double a, double b, double c, double eps)
/*
2 : x[0], x[1]
1 : x[0] ± x[1]i
*/
{
	if(too_small(a, eps))
		a = eps;
	double D = b*b -4*a*c, a_t=0.5*a;
	if(D >= 0)
	{
		double D_1=sqrt(D);
		x[0] = (D_1 - b)*a_t;
		x[1] = (-D_1 - b)*a_t;
		//cout<<"x: \n"<<endl<<x[0]<<endl<<x[1]<<endl;
		return 2;
	}
	else
	{
		x[0] = -b*a_t;
		x[1] = sqrt(-D)*a_t;
		return 1;
	}
}
int complex_sol2(double *x, double a, double b_re, double b_im, double c_re, double c_im)
{
	complex<double> b(b_re, b_im), c(c_re,c_im), D = b*b -4*a*c, x1=-b + sqrt(D), x2=-b -sqrt(D);
	x1.real(x1.real()*0.5);
	x1.imag(x1.imag()*0.5);
	x2.real(x2.real()*0.5);
	x2.imag(x2.imag()*0.5);
	x[0]=x1.real();
	x[1]=x1.imag();
	x[2]=x2.real();
	x[3]=x2.imag();
	if(x[1] == 0)
		return 1;
	if(x[3] == 0)
	{
		double temp = x[2];
		x[2]=x[0];
		x[0]=temp;
		x[3]=x[1];
		x[1] = 0;
		return 1;
	}
	else
		return 0;
}
int sol3(double *x, double a, double b, double c, double d, double eps)
/*
ax^3 + bx^2 + cx + d = 0
3 : x[0], x[1], x[2]
1 : x[0], x[1] ± x[2]i
*/
{
	if(too_small(a, eps))
		a = eps;
	//y^3 + py + q = 0
	double p, q, two_seven = 1./27., a_3=1./(a*a*a), one_third = 1./3., a_1=1./a, root_3 = sqrt(3);
	p=(3*a*c - b*b)*one_third*a_3*a;
	q=(2*b*b*b - 9*b*c*a + 27 *a*a*d)*two_seven*a_3;
	//y=z-p/(3z) -> z^6 + qz^3 - p^3/27 = 0 || D = b^2 - 4ac -> z^3 = (+-sqrt(D) - q)/2
	double z[2];
	int solves_z = sol2(z, 1, q, -p*p*p*two_seven, eps);
	if(too_small(z[0], eps))
		z[0] = eps;
	if(too_small(z[1], eps))
		z[1] = eps;
	if(solves_z == 2)
	{
		double y[2];
		z[0] = cbrt(z[0]);
		z[1] = cbrt(z[1]);
		y[0] = z[0] - p/z[0]*one_third;
		y[1] = z[1] - p/z[1]*one_third;
		x[0] = y[0] - b*a_1*one_third;
		x[1] = y[1] - b*a-1*one_third;
		double
		z_re=z[0]*-0.5,
		z_im=z[0]*root_3*0.5,
		y_re=z_re - p*z_re/(z_re*z_re + z_im*z_im)*one_third,
		y_im=z_im + p*z_im/(z_re*z_re + z_im*z_im)*one_third;

		x[1]=y_re - b*one_third*a_1;
		x[2]=y_im;
		return 1;
	}
	else if(solves_z = 1)
	{
		if(too_small(p, eps))
		p = eps;
		double cosing = 3*root_3*q/(p*sqrt(-p)*0.5),
		ang = acos(cosing)*180./M_PI,
		psy[3] = {ang*one_third, ang*one_third + 120, ang*one_third +240},
		r = sqrt(abs(p)*one_third);
		for(int i = 0; i < 3; ++i)
			x[i] = (r - p/r*one_third)*cos(psy[i]*M_PI/180);
		return 3;
	}
	return 0;

}
int sol4(double *x, double a, double b, double c, double d, double e, double eps)
/*
ax^4 + bx^3 +cx^2 + dx + e => y^4 + py^2 + r
8w^3 - 4pw^2  - 8rw + (4pr - q^2) = 0
4 :	x[0], x[1], x[2], x[3] 
2 : x[0], x[1], x[2] ± x[3]i
0 : x[0] ± x[1], x[2] ± x[3]i
*/
{
	if(too_small(a, eps))
		a = eps;
	//y^4 + py^2 + r
	double p,q,r, a_4=1./(a*a*a*a), one_eight = 1./8, a_1=1./a;
	p = (8*a*c - 3*b*b)*a*a*a_4*one_eight;
	q=(b*b*b - 4*b*c*a + 8*d*a*a)*one_eight*a_4*a;
	r = (256*a*a*a*e - 3*b*b*b*b + 16*a*b*b*c - 64*a*a*b*d)*one_eight*one_eight*one_eight*2*a_4;
	double w[3];
	int solves_w = sol3(w, 8, -4*p, -8*r, 4*p*r - q*q, eps), solves_1, solves_2;

	//cout<<solves_w<<endl;
	//cout<<w[0]<<" "<<w[1]<<" "<<w[2]<<endl;
	//cout<<2*w[0] - p<<endl;
	if(solves_w == 3)
	{
		int i =	0;
		//for(; i < 3; ++i)
		//cout<<i<<" "<<w[i]<<" "<<2*w[i] - p<<endl;
		double f = sqrt(2*w[i] - p);
		if(too_small(f, eps))
			f = eps;
		double g = -q/f*0.5, y[2];
		solves_1 = sol2(y, 1, -f, w[i] - g, eps);
		//cout<<y[0]<<" "<<y[1]<<endl;
		x[0] = y[0] - b*0.25*a_1;
		x[1] = y[1];
		if (solves_1  == 2) x[1]-= b*0.25*a_1;
		solves_2 = sol2(y, 1, f, w[i] + g, eps);
		x[2] = y[0] - b*0.25*a_1;
		x[3] = y[1];
	}
	else
	{
		double f = sqrt(2*w[0] - p);
		if(too_small(f, eps))
		f=eps; 
		double g = -q/f*0.5, y[2];
		solves_1 = sol2(y, 1, -f, w[0] - g, eps);
		//cout<<y[0]<<" "<<y[1]<<endl;
		x[0] = y[0] - b*0.25*a_1;
		x[1] = y[1];
		if (solves_1  == 2) x[1]-= b*0.25*a_1;
		solves_2 = sol2(y, 1, f, w[0] + g, eps);
		x[2] = y[0] - b*0.25*a_1;
		x[3] = y[1];
	}
	if(solves_2 == 2) x[3]-=b*0.25*a_1;

	if(solves_1 == 2 && solves_2 == 2)
		return 4;
	else if(solves_1 == 2 && solves_2 == 1)
		return 2;
	else if(solves_1 == 1 && solves_2 == 1)
		return 0;
	else
	{
		double temp;
		temp=x[0];
		x[0]=x[2];
		x[2]=temp;

		temp=x[1];
		x[1]=x[3];
		x[3]=temp;
	}
	//complex_sol2(x, 1, 2,
}
string quadric_poly(double a, double b, double c, double eps)
{
	string res="There are ";
	double x[2];
	int solution = sol2(x, a, b, c, eps);
	if (solution == 2)
	{
		res+="2 real sloutions:\n";
		res+=to_string(x[0]) + "\n" + to_string(x[1]) + "\n"; 
	}
	else
		res+="2 complex solutions:\n" + to_string(x[0]) + " + " + to_string(fabs(x[1])) +
		"\n" + to_string(x[0]) + " - " + to_string(fabs(x[0])) + "i\n";
}
string qubic_poly(double a, double b, double c, double d, double eps)
{

	string res="There are ";
	double x[3];
	int solution = sol3(x,a,b,c,d, eps);
	if(solution == 3)
	{
		res += "3 real solutions:\n";
		for(int i = 0; i < 3; ++i)
			res+=to_string(x[i]) + "\n";
	}
	if(solution == 1)
	{
		res+="1 real, 2 complex solutions:\n";
		res+=to_string(x[0]) + "\n";
		res+=to_string(x[1]) + " + " + to_string(fabs(x[2])) + "i\n" 
		+ to_string(x[1]) + " - " + to_string(fabs(x[2])) + "i\n"; 
	}
	return res;
}
string quatic_poly(double a, double b, double c, double d, double e, double eps)
{
	string res="Threre are ";
	double x[4];
	int solution = sol4(x, a, b, c, d, e, eps);
	if(solution == 4)
	{
		res+="4 real solutions:\n";
		for(int i = 0; i < 4; ++i)
			res+=to_string(x[i])+"\n";
	}
	else if(solution == 2)
	{
		res+="2 real, 2 complex solutions:\n";
		res+=to_string(x[0]) + "\n";
		res+=to_string(x[1]) + "\n";
		res+= to_string(x[2]) + " + " + to_string(fabs(x[3])) + "i\n" + 
		to_string(x[2]) + " - " + to_string(fabs(x[3])) + "i\n"; 
	}
	else 
	{
		res+="4 complex solutions:\n";
		res+= to_string(x[0]) + " + " + to_string(fabs(x[1])) + "i\n" +
		to_string(x[0]) + " - " + to_string(fabs(x[1])) + "i\n"; 
		res+= to_string(x[2]) + " + " + to_string(fabs(x[3])) + "i\n" +
		to_string(x[2]) + " - " + to_string(fabs(x[3])) +"i\n"; 
	}
	return res;
}
int main()
{

	double x[4], eps=0.000001;
	/*
	double a=0, b =-9;
	complex_sqrt(a,b);
	cout<<a<<" "<<b;
	cout<<complex_sol2(x, 1, 0, 2, -1, 2)<<endl;
	cout<<x[0]<<" "<<x[1]<<endl<<x[2]<<" "<<x[3]<<endl;
	*/
	cout << quatic_poly(2, 2, -3, -3, -4, eps);
}