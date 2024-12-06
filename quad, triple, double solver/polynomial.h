#include <fstream>
#include <string>
#include <math.h>
#include <cstdarg>
#include <complex>
#include <set>
#include <vector>
#include <iostream>
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
double safe_acos(double a) 
{
    if (a <= -1.0) 
        return M_PI;
    else if (a >= 1.0) 
        return 0;
    else 
        return acos(a);
}
int correct(double *x, int size)
{
	int new_size=0, i;
	set <double> res;
	for(i = 0; i < size; ++i)
		res.insert(x[i]);
	i=0;
	for(set <double> :: iterator it = res.begin(); it != res.end(); ++it, new_size++, ++i)
		x[i]=*it;
	return new_size;
}
int sol2(double *x, double a, double b, double c, double eps)
/*
2 : x[0], x[1]
1 : x[0] ± x[1]i
*/
{
	if(too_small(a, eps))
		a = eps;
	double D = b*b -4*a*c, a_t=0.5/a;
	//cout<<"D :"<<D<<endl;
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
	if(x[2] == -0)
		x[2] =0;
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
int complex_solve_2(double *x, complex <double> b, complex <double> c)
{
	complex <double> D_1 = sqrt(b*b*0.25 - c), x1 = D_1 - b*0.5, x2 = -D_1 - b*0.5;
	x[0]=x1.real();
	x[1]=x1.imag();
	x[2]=x2.real();
	x[3]=x2.imag();
	if(x[1] == 0 && x[3] == 0)
		return 2;
	else if(x[1] == 0)
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
	if(d == 0)
	{
		double y[2], res1=sol2(y, a, b, c, eps);
		//cout<<y[0]<<" "<<y[1]<<endl;
		x[0]=0;
		x[1]=y[0];
		x[2]=y[1];
		if(res1 == 2)
			return 3;
		else
			return 1;
	}
	if(too_small(a, eps))
		a = eps;
	//y^3 + py + q = 0
	double D, p, q, two_seven = 1./27., a_3=1./(a*a*a), one_third = 1./3., a_1=1./a, root_3 = sqrt(3);
	p=(3.0*a*c - b*b)*one_third*a_3*a;
	q=(2.0*b*b*b - 9*b*c*a + 27 *a*a*d)*two_seven*a_3;
	//D = q*q*0.25 + p*p*p*one_third*one_third*one_third;
	//cout<<"D_1 :"<<D<<endl;
	/*
	p = -b*b/(3*a*a) + c/a; 
	q = 2*b*b*b/(27*a*a*a) - b*c/(3*a*a) + d/a;
	*/
	//cout<<"pq "<<p<<" "<<q<<endl;
	if(p == 0 && q == 0)
	{	
		x[0]= -b*one_third*a_1;
		x[2]=x[0];
		x[1]=x[2];
		return -1;
	}
	else if(p == 0 && q != 0)
	{
		x[0] = cbrt(-q) - b*one_third*a_1;
		x[1] = cbrt(q)*0.5;
		x[2] = cbrt(q)*root_3*0.5;
		return 1;
	}
	else if(p != 0 && q == 0)
	{
		x[0] = -b*one_third*a_1;
		if(p <= 0)
		{
			x[1]=sqrt(-p) - b*one_third*a_1;
			x[2]=-sqrt(-p) -b*one_third*a_1;
			return 3;
		}
		else
		{
			x[1]=0;
			x[2]=sqrt(p);
			return 1;
		}
	}
	else
	{
		//y=z-p/(3z) -> z^6 + qz^3 - p^3/27 = 0 || D = b^2 - 4ac -> z^3 = (+-sqrt(D) - q)/2
		double z[2];
		int solves_z = sol2(z, 1, q, -p*p*p*two_seven, eps);
		//cout<<q<<" "<<-p*p*p*two_seven<<"|||||||\n";
		//cout<<z[0]<<" "<<z[1]<<endl;	
			if(solves_z == 2)
		{
			//cout<<z[0]<<" "<<z[1]<<endl;
			double y[2];
			z[0] = cbrt(z[0]);
			z[1] = cbrt(z[1]);
			y[0] = z[0] - p/z[0]*one_third;
			y[1] = z[1] - p/z[1]*one_third;
			x[0] = y[0] - b*a_1*one_third;
			x[1] = y[1] - b*a_1*one_third;
			double
			z_re=z[0]*-0.5,
			z_im=z[0]*root_3*0.5,
			y_re=z_re - p*z_re/(z_re*z_re + z_im*z_im)*one_third,
			y_im=z_im + p*z_im/(z_re*z_re + z_im*z_im)*one_third;

			x[1]=y_re - b*one_third*a_1;
			x[2]=y_im;
			if(fabs(x[2]) <= eps)
				return -2;
			return 1;
		}
		else if(solves_z = 1)//x= a +_ bi
		{
			//if(too_small(p, eps))
				//p = eps;
			double cosing = (1.5)*root_3*q/(p*sqrt(-p)),

			ang = acos(cosing)*(180./M_PI),
			psy[3] = {ang*one_third, ang*one_third + 120., ang*one_third +240.},
			r = sqrt(fabs(p)*one_third);
			//cout<<"cos, acos ang r "<<cosing<< " "<<ang<<" "<<r<<endl;
			//cout<<"cosing "<<cosing<<" ang"<<ang<<" "<<"r "<<r<<endl;;
			for(int i = 0; i < 3; ++i)
				x[i] = (r - one_third*p/r)*cos(psy[i]*(M_PI/180)) - b*one_third*a_1;
			return 3;
			
		}
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
	double a_1=1./a;
	if (b == 0 && d == 0)//ax^4 + cx^2 + e = 0 => (x^2 + b_1x + c_1)(x^2 + b_2x + c_2)
	{	
		c*=a_1;
		e*=a_1;
		complex <double> b_1, c_1(e, 0);
		c_1=sqrt(c_1);
		b_1=c_1;
		b_1.real(2*b_1.real() - c);
		b_1.imag(2*b_1.imag());
		b_1=sqrt(b_1);
		double complex_x1[4]={}, complex_x2[4]={};
		int res1 = complex_solve_2(complex_x1, b_1, c_1), res2 = complex_solve_2(complex_x2, -b_1, c_1);
		if(res1 == 2)
		{
			if(res2 == 2)
			{
				x[0]=complex_x1[0];
				x[1]=complex_x1[2];
				x[2]=complex_x2[0];
				x[3]=complex_x2[2];
				return 4;
			}
			else
			{
				x[0]=complex_x1[0];
				x[1]=complex_x1[2];
				x[2]=complex_x2[0];
				x[3]=complex_x2[1];
				return 2;
			}
		}
		else if(res1 == 1)
		{
			if(res2 == 2)
			{
				x[0]=complex_x1[0];
				x[1]=complex_x2[0];
				x[2]=complex_x1[2];
				x[3]=complex_x1[3];
				return 2;
			}
			else
			{
				x[0]=complex_x1[0];
				x[1]=complex_x2[0];
				x[2]=complex_x1[2];
				x[3]=complex_x1[3];
				return 2;
			}
		}
		else
		{
			x[0]=complex_x1[0];
			x[1]=complex_x1[1];
			x[2]=complex_x2[0];
			x[3]=complex_x2[1];
			return 0;
		}
	}
	else if (d == 0 && e == 0)
	{
		double y[2]={};
		int res = sol2(y, a, b, c, eps);
		x[0] = 0;
		x[1] = y[0];
		x[2] = y[1];
		if(res == 2)
			return 3;
		else if(res == 1)
			return 1;
	}
	else if (e == 0)
	{
		bool is_zero=false;
		int res = sol3(x, a, b, c, d, eps);
		switch(res)
		{
			case 3:
			{
				for(int i = 0; i < 3; ++i)
				{
					if(x[i] == 0)
						is_zero=true;
				}
				if(is_zero)
					return 3;
				x[3] = 0;
				if(!is_zero) 
					return 4;
				break;
			}
			case -2:
			{
				is_zero = false;
				for(int i = 0; i < 2; ++i)
				{
					if(x[i] == 0)
						is_zero=true;
				}
				if(is_zero)
					return -2;
				x[2] = 0;
				if(!is_zero) 
					return 3;
				break;
			} 
			case 1:
			{
				if(x[0] == 0)
					return 1;
				else
				{
					x[3]=x[2];
					x[2]=x[1];
					x[1]=x[0];
					x[0]=0;
					return 2;
				}
				break;
			}
			case -1:
			{
				if (x[0] == 0)
					return -1;
				else
				{
					x[1] = x[0];
					x[0] = 0;
					return -2;
				}
				break;
			}
			default:
				break;
		}
	}
	
	//y^4 + py^2 + r
	double p,q,r, a_4=1./(a*a*a*a), one_eight = 0.125, one_third=1./3;
	
	p = (8*a*c - 3*b*b)*a*a*a_4*one_eight;
	q=(b*b*b - 4*b*c*a + 8*d*a*a)*one_eight*a_4*a;
	r = (256*a*a*a*e - 3*b*b*b*b + 16*a*b*b*c - 64*a*a*b*d)*one_eight*one_eight*one_eight*2*a_4;
	//cout<<"pqr "<<p<<" "<<q<<" "<<r<<endl;
	//cout<<-4*p<<" "<<-8*r<<" "<<4*p*r - q*q<<"||||||||||\n";
	if(p == 0 && q == 0 && r == 0)
	{
		x[0] = -b*0.25*a_1;
		return -1;
	}
	else if(q == 0)
	{
		double y[4]={};
		int res = sol4(y, 1, 0, p, 0, r, eps);
		if(res == 4)
		{
			for(int i = 0; i < 4; ++i)
				x[i] = y[i] -b*0.25*a_1;
			return 4;
		}
		else if (res == 2)
		{
			x[0] = y[0] - b*0.25*a_1;
			x[1] = y[1] - b*0.25*a_1;
			x[2] = y[2] - b*0.25*a_1;
			x[3]=y[3];
			return 2;
		}
		else if (res == 0)
		{
			x[0] = y[0] - b*0.25*a_1;
			x[1] = y[1];
			x[2] = y[2] - b*0.25*a_1;
			x[3] = y[3];
			return 0;
		}


	}
	double w[3];
	int solves_w = sol3(w, 8, -4*p, -8*r, 4*p*r - q*q, eps), solves_1, solves_2;

	//cout<<solves_w<<endl;
	//cout<<w[0]<<" "<<w[1]<<" "<<w[2]<<endl;
	//cout<<2*w[0] - p<<endl;
	if(solves_w == 3)
	{
		//cout<<111111<<endl;
		if(2*w[0] - p >= 0)
		{
			//for(; i < 3; ++i)
			//cout<<i<<" "<<w[i]<<" "<<2*w[i] - p<<endl;
			//cout<<w[0]<<endl<<w[1]<<endl<<w[2]<<endl;
			//cout<<2*w[0] - p<<endl;
			double f = sqrt(2*w[0] - p);
			if(too_small(f, eps))
				f = eps;
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
		else
		{
			int i =1;
			while(i < 3 && 2*w[i] - p < 0)
				++i;
			if(i == 3)
			{
				i=0;
				double y[4];
				solves_1 = complex_sol2(y, 1, 0, -sqrt(p - 2*w[i]), 0, 0.5*q/sqrt(p-2*w[i]));
			}
			else
			{
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

			//cout<<-sqrt(p - 2*w[0])<<" "<< 0.5*q/sqrt(p-2*w[0])<<endl;
			//for(int i = 0; i < 4; ++i)
				//cout<<y[i]<<",";
			//cout<<endl;
		}
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
	{
		int res = correct(x, 4);
		if(res == 1 || res == 2)
			return -res;
		return res;
	}
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
		return  2;
	}
}
string get_str_W(const vector <double> a)
{
    string res="";
    for(int i = 0; i < a.size(); ++i)
    {
        if(i == 0 && a[i] < 0)
        {    
            res+="-";
            if(fabs(a[i]) != 1)
                res+=to_string(fabs(a[i]));
            res+="x^" + to_string(a.size()-1-i);
        }
        else if(i == 0 && a[i] > 0)
        {
            if(fabs(a[i]) != 1)
                        res+=to_string(fabs(a[i]));
                    res+="x^" + to_string(a.size()-1-i);
        }
        else if(i == a.size() - 1)
        {
            if(a[i] != 0)
            {
                    if(a[i] > 0)
                {
                    res+=" + ";
                    res+=to_string(fabs(a[i]));
                }
                else if(a[i] < 0)
                {
                    res+=" - ";
                    res+=to_string(fabs(a[i]));
                }
            }
        }
        else if(a[i] != 0)
        {
            if(a[i] > 0)
            {
                res+=" + ";
                if(fabs(a[i]) != 1)
                    res+=to_string(fabs(a[i]));
                res+="x";
                if(i != a.size() - 2)
                	res+="^" + to_string(a.size()-1-i);
            }
            else if(a[i] < 0)
            {
                res+=" - ";
                if(fabs(a[i]) != 1)
                    res+=to_string(fabs(a[i]));
                res+="x";
                if(i != a.size() - 2)
                	res+="^" + to_string(a.size()-1-i);
            }
        }
    }
    res+=" = 0\n";
    return res;
}
string quadric_poly(double a, double b, double c, double eps)
{
	vector <double> coef(3);
	coef[0] = a;
	coef[1] = b;
	coef[2] = c;
	string res = "For equation:\n";
	res+=get_str_W(coef);
	res+="There are ";
	double x[2];
	int solution = sol2(x, a, b, c, eps);
	if (solution == 2)
	{
		res+="2 real sloutions:\n";
		res+=to_string(x[0]) + "\n" + to_string(x[1]) + "\n"; 
	}
	else
		res+="2 complex solutions:\n" + to_string(x[0]) + " + " + to_string(fabs(x[1])) +
		"\n" + to_string(x[0]) + " - " + to_string(fabs(x[1])) + "i\n";
	return res;
}
string cubic_poly(double a, double b, double c, double d, double eps)
{
	vector <double> coef(4);
	coef[0] = a;
	coef[1] = b;
	coef[2] = c;
	coef[3] = d;
	string res = "For equation:\n";
	res+=get_str_W(coef);
	res+="There ";
	double x[3];
	int solution = sol3(x,a,b,c,d, eps);
	if(solution == 3)
	{
		res += "are 3 real solutions:\n";
		for(int i = 0; i < 3; ++i)
			res+=to_string(x[i]) + "\n";
	}
	else if(solution == 1)
	{
		res+="are 1 real, 2 complex solutions:\n";
		res+=to_string(x[0]) + "\n";
		res+=to_string(x[1]) + " + " + to_string(fabs(x[2])) + "i\n" 
		+ to_string(x[1]) + " - " + to_string(fabs(x[2])) + "i\n"; 
	}
	else if(solution == -1)
	{
		res+="is 1 real solution:\n";
		res+=to_string(x[0]) + "\n";
	}
	else if(solution == -2)
	{
		res+="are 2 real solutions:\n";
		res+=to_string(x[0]) + "\n" + to_string(x[1]) + "\n";
	}
	return res;
}
string quatic_poly(double a, double b, double c, double d, double e, double eps)
{
	vector <double> coef(5);
	coef[0] = a;
	coef[1] = b;
	coef[2] = c;
	coef[3] = d;
	coef[4] = e;
	string res = "For equation:\n";
	res+=get_str_W(coef);
	res+="Threre ";
	double x[4];
	int solution = sol4(x, a, b, c, d, e, eps);
	cout<<solution<<endl;
	switch(solution)
	{
		case 4:
		{
			res+="are 4 real solutions:\n";
			for(int i = 0; i < 4; ++i)
				res+=to_string(x[i])+"\n";
			break;
		}
		case 3:
		{
			res+="are 3 real solutions:\n";
			res+= to_string(x[0]) + "\n" + to_string(x[1]) + "\n" + to_string(x[2]) + "\n";
			break;
		}
		case 2:
		{
			res+="are 2 real, 2 complex solutions:\n";
			res+=to_string(x[0]) + "\n";
			res+=to_string(x[1]) + "\n";
			res+= to_string(x[2]) + " + " + to_string(fabs(x[3])) + "i\n" + 
			to_string(x[2]) + " - " + to_string(fabs(x[3])) + "i\n"; 
			break;
		}
		case 1:
		{
			res+="are 1 real, 2 complex solutions:\n";
			res+=to_string(x[0]) + "\n" + to_string(x[1]) + " + " + to_string(fabs(x[2])) +"\n"
			+ to_string(x[1]) + " - " + to_string(fabs(x[2])) + "\n";
			break;
		}
		case 0:
		{
			res+="are 4 complex solutions:\n";
			res+= to_string(x[0]) + " + " + to_string(fabs(x[1])) + "i\n" +
			to_string(x[0]) + " - " + to_string(fabs(x[1])) + "i\n"; 
			res+= to_string(x[2]) + " + " + to_string(fabs(x[3])) + "i\n" +
			to_string(x[2]) + " - " + to_string(fabs(x[3])) +"i\n"; 
			break;
		}
		case -1:
		{
			res+="is 1 real solution:\n";
			res+=to_string(x[0])+ "\n";
			break;
		}
		case -2:
		{
			res+="are 2 real solutions:\n";
			res+=to_string(x[0]) + "\n" + to_string(x[1]) + "\n";
			break;
		}
		default:
		{
			res+="No solutions found\n";
			break;
		}
	}
	return res;
}