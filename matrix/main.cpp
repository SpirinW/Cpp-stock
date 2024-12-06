#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdarg.h>
using namespace std;
class Matrix
{
	vector <vector <double>> a;
public:
	Matrix(unsigned int n, unsigned int m, ...)
	{
		va_list args;
		va_start(args, m);
		for(int i = 0; i < n; ++i)
		{
			vector <double> temp;
			for(int j  = 0; j < m; ++j)
				cout<<va_arg(args, double;
				//temp.push_back(va_arg(args, double));
			a.push_back(temp);
		}
	}
	void print()
	{
		for(int i = 0; i < a.size(); ++i)
		{
			for(int j = 0; j < a[i].size(); ++j)
				cout<<a[i][j];
			cout<<endl;
		}
	}
	ostream& operator <<(ostream& os)
	{
		for(int i = 0; i < a.size(); ++i)
		{
			for(int j = 0; j < a[i].size(); ++j)
				os<<a[i][j]<<" ";
			os<<endl;
		}
		return os;
	}
};

int main()
{
	Matrix a(2, 3, 1,2,3,4,5,6);
	//a.print();
	//cout<<a;

	return 0;
}