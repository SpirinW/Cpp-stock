#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include "polynomial.h"
using namespace std;
bool check_eq(double a, double b, double e)
{
    return fabs(a-b) < e;
}
bool test_eq(const vector <double> a, const vector <double> ans, double eps)
{
    vector <bool> count;
    bool result=true;
    if(a.size() == 4)
    {//x^3
        double x[3], solve = sol3(x, a[0], a[1], a[2], a[3], eps);
        if(solve == 3)//3 real
        {

            for(int i = 0; i < ans.size(); ++i)
            {
                bool found = false;
                for(int k = 0; k < 3; ++k)
                {
                    if(check_eq(x[k], ans[i], eps))
                        found =true;
                }
                if(found)
                    count.push_back(true);
                else
                    count.push_back(false);
            }
        }
        else if(solve == 1)//1 real, 2 complex
        {

            count.push_back(check_eq(x[0], ans[0], eps));
            count.push_back(check_eq(x[1], ans[1], eps));
            count.push_back(check_eq(x[2], ans[2], eps));
        }
        for(int i = 0; i < count.size(); ++i)
            result*=count[i];
        return result;
    }
    else if(a.size() == 5)
    {//x^4
        double x[4], solve = sol4(x, a[0], a[1], a[2], a[3], a[4], eps);
        if(solve == 4)//4 real
        {

            for(int i = 0; i < ans.size(); ++i)
            {
                bool found = false;
                for(int k = 0; k < 4; ++k)
                {
                    if(check_eq(x[k], ans[i], eps))
                        found =true;
                }
                if(found)
                    count.push_back(true);
                else
                    count.push_back(false);
            }
        }
        else if(solve == 2)//
        for(int i = 0; i < count.size(); ++i)
            result*=count[i];
        return result;
    }
}

void run_test_write(double eps)
{
    ifstream file;
    file.open("input.txt");
    ofstream res;
    res.open("test_input.txt");
    int type;
    int solution;
    vector <double> a(5);
    for(int test = 1; file>>type; test++)
    {
        cout<<"Test "<<test<<":\n";
        res<<type<<"| ";
        if(type == 4)
        {
            for (int i = 0; i < 5; ++i)
                file>>a[i];
            double x[4];
            solution = sol4(x, a[0], a[1], a[2], a[3], a[4], eps);
            cout<<quatic_poly(a[0], a[1], a[2], a[3], a[4], eps);
            for(int i = 0; i < 5; ++i)
                res<<a[i]<<" ";
            res<<"| ";
            switch(solution)
            {
                case 4: case 3: case -1: case -2:
                {
                    res<<abs(solution)<<" : ";
                    for(int i = 0; i < abs(solution); ++i)
                        res<<x[i]<<" ";
                    break;
                }
                case 1: case 2:
                {
                    res<<solution<<";2 : ";
                    for(int i = 0; i < solution; ++i)
                        res<<x[i]<<" ";

                    break;
                }   
                case 0:
                    res<<solution<<";4 : ";
                    for(int i = 0; i < 4; ++i)
                        res<<x[i]<<" ";
                    break;
                default:
                    break;
            }
            res<<endl;
        }
        else if(type == 3)
        {
            for (int i = 0; i < 4; ++i)
                file>>a[i];
            double x[3]={};
            solution = sol3(x, a[0], a[1], a[2], a[3], eps);
            cout<<cubic_poly(a[0], a[1], a[2], a[3], eps);
            for(int i = 0; i < 4; ++i)
                res<<a[i]<<" ";
            res<<"| ";
            
            if (solution == 3)
            {
                res<<"3 : ";
                for(int i = 0; i < 3; ++i)
                    res<<x[i]<<" ";
            }
            else if (solution == 1)
            {
                res<<"1;2 : ";
                for(int i = 0; i < 3; ++i)
                    res<<x[i]<<" ";
            }
            else if (solution == -1 || solution == -2)
            {
                res<<abs(solution)<<" : ";
                for(int i = 0; i < abs(solution); ++i)
                    res<<x[i]<<" ";
            }
            res<<endl;
        }
        cout<<endl;
    }
    file.close();
    res.close();
}
void check_tests(double eps)
{
    ifstream file;
    file.open("input.txt");
    ofstream res;
    res.open("test_results.txt");
    int type;
    double solution;
    vector <double> a(5);
    for(int test = 1; file>>type; test++)
    {
        cout<<"Test "<<test<<":\n";
        if(type == 4)
        {
            for (int i = 0; i < 5; ++i)
                file>>a[i];
            double x[4];
            //solution = sol4(x, a[0], a[1], a[2], a[3], a[4], eps);
            cout<<quatic_poly(a[0], a[1], a[2], a[3], a[4], eps);
            //res<<quatic_poly(a[0], a[1], a[2], a[3], a[4], eps);

        }
        else if(type == 3)
        {
            for (int i = 0; i < 4; ++i)
                file>>a[i];
            double x[3]={};
            //solution = sol3(x, a[0], a[1], a[2], a[3], eps);
            cout<<cubic_poly(a[0], a[1], a[2], a[3], eps);
        }
    }
    file.close();
    res.close();

}
int main(int argc, char const *argv[])
{
    double eps=0.000001, x[4]={};
    /*
    inaccurate

    cout<<"===============\n";
    cout<< cubic_poly(4, 3, 2, 1, eps);
    cout<<cubic_poly(1, 0, -2, 1, eps);
    cout<<cubic_poly(1, 3, -24, 1, eps);
    cout <<cubic_poly(1, 0, 0, 0, eps);
    cout<<cubic_poly(2, -3,-3,2,eps);
    ------
    cout<<cubic_poly(8, -18, -38.5, 86.625, eps);
    -----------
    
    cout<<quadric_poly(3, 4, 2, eps);
    run_test_3(eps);
    cout<<quadric_poly(1, 2, -3, eps);
    cout<<cubic_poly(5, -2, -1, 1, eps);
    cout<<quatic_poly(1, 0, 0, -3, 1, eps);
    vector <double> a={-1, 2, -3, 0, 5};
    cout<<get_str_W(a);
    cout<<cubic_poly(1, -2, -1, 1, eps);
    cout<<cubic_poly(5, -8, -8, 5, eps);
    cout<<quatic_poly(1, 2, 0, -2, -1, eps);
    cout<<quatic_poly(1, 4, 6, 4, 1, eps);
    cout<<quatic_poly(1, -2,  3, -2, 2, eps);
    check_tests(eps);
    run_test_write(eps);
    cout<<cubic_poly(5, -1, -20, 4, eps);
    cout<<"======================\n";
    cout<<cubic_poly(5, -8, -8, 5, eps);
    cout<<"======================\n";
    cout<<cubic_poly(1, -2, -1, 1, eps);
    cout<<"======================\n";
    cout<<quatic_poly(1, 1 ,1 ,1 ,1, eps);
    cout<<"======================\n";
    cout<<quatic_poly(7, 10, 4, 8, 10, eps);
    cout<<"======================\n";
    cout<<cubic_poly(5, -1, 20, 4, eps);
    cout<<cubic_poly(1, -2, -1, 1, eps);
    complex <double> b(1, 0), c(-1, 0);
    cout<<complex_solve_2(x, -b, c);
    for (int i = 0; i < 4; ++i)
        cout<<"| "<<x[i]<<endl;
    cout<<cubic_poly(1, -6, 11, -6, eps)<<endl;
    cout<<quatic_poly(2, 0, 5, 0, -3, eps);
    */
    cout<<cubic_poly(1, -2, -1, 1, eps);
    cout<<cubic_poly(5, -8, -8, 5, eps);
    cout << quatic_poly(2, 2, -3, -3, -4, eps);
    cout<< quatic_poly(1, 2, 6, 5, 6, eps);
    return 0;
}