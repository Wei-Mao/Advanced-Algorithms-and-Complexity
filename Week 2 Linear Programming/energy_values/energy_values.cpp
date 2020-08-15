#include <iostream>
#include <algorithm>  // std::swap
#include <vector>
#include <limits>
using namespace std;

const double EPS = 1e-6;
const int PRECISION = 20;
const int int_min = std::numeric_limits<int>::min();

typedef std::vector<long double> Column;
typedef std::vector<long double> Row;
typedef std::vector<Row> Matrix;

struct Position
{
  int column;
  int row;
};

class Equation {
public:
    Matrix a;
    Column b;
    /* Equation(const Matrix &a, const Column &b): */
    /*     a(a), */
    /*     b(b) */
    /* {} */
};

class Gaussian_Elimination_New
{

public:
    /*
       Coefficient Matrix is a square matrix.
    */
    Equation Linear_Eqn;
    int size;

public:
    void ReadEquation()
    {
        std::cin >> size;
        Matrix a(size, std::vector <long double> (size, 0.0));
        Column b(size, 0.0);
        for (int row = 0; row < size; ++row) {
            for (int column = 0; column < size; ++column)
                std::cin >> a[row][column];
            std::cin >> b[row];
        }
        Linear_Eqn.a = a;
        Linear_Eqn.b = b;
    }


    void Row_Exchange(int i, int j)
    {
        // swap Coefficient matrix
        swap(Linear_Eqn.a[i], Linear_Eqn.a[j]);

        // swap right-hand side of the linear equation.
        swap(Linear_Eqn.b[i], Linear_Eqn.b[j]);
    }

    int pivoting(Position &pivot_element, vector<bool> &used_row, vector<bool> &used_col)
    {
        /*
	   Assume all rows of zeros lie at the bottom of augmented matrix.
	   A = [1 2 3 4
		1 2 3 3
		0 0 0 0
		0 0 0 0]
	   b = [1
		2
		0
		0]
        */

	// Find the pivot position
	pivot_element.row = 0;
	pivot_element.column = 0;

	while(used_row[pivot_element.row])
	  pivot_element.row++;
	
	while(used_col[pivot_element.column])
	  pivot_element.column++;

        // Put the largest entry in the pivot column on the pivot position
        int max_value = -1;
        int max_idx = -1;

        for(int j(pivot_element.row); j < Linear_Eqn.a.size(); j++)
        {
            if(abs(Linear_Eqn.a[j][pivot_element.column]) > max_value)
            {
                max_value = abs(Linear_Eqn.a[j][pivot_element.column]);
                max_idx = j;
            }
        }
	int max_row = max_idx;
	return max_row;
	// Selecting pivot in such away still gives the possibility of zero pivot
        // Swap the current row with the row with largest entry in the pivot column.
    }

    void elimination(int pivot_row, int pivot_col)
    {
        /*
       Subtracting multiplier of pivot_row from current_row.
        */

        long double pivot_value = Linear_Eqn.a[pivot_row][pivot_col];

        for (int i(pivot_row + 1); i < Linear_Eqn.a.size(); i++)
        {
            long double multiplier = Linear_Eqn.a[i][pivot_col] / pivot_value;
            Linear_Eqn.b[i] -= multiplier * Linear_Eqn.b[pivot_row];
            for (int j(pivot_col); j < Linear_Eqn.a[0].size(); j++)
            {
                Linear_Eqn.a[i][j] -= multiplier * Linear_Eqn.a[pivot_row][j];
            }
        }
    }

    int back_substitution(vector<long double> &sol, vector<Position> pivots)
    {
	int status = 1; // 1 -unique solution. 0- no solution. 2: infinitely many solutions


	if (pivots.size() == 0)
	{
	  for(int k(0); k < Linear_Eqn.b.size() ; k++)
	  {
	    if (abs(Linear_Eqn.b[k]) > EPS)
	    {
	      return 0;
	    }
	  }

	  for (int k(0); k < Linear_Eqn.b.size(); k++)
	  {
	    sol.push_back(1e9);
	  }

	}
	if (pivots.size() >0)
	{
	int bottom_pivot = pivots[pivots.size()-1].row;

	if (bottom_pivot < Linear_Eqn.a.size()-1) {status = 2;}
	for (int k(Linear_Eqn.a.size()-1); k > bottom_pivot; k--)
	{
	  if (abs(Linear_Eqn.b[k]) > EPS)
	  {
	    status = 0;
	  }
	}

	if (status == 0)
	{
	  // No solution. It makes no sense to do back_substitution.
	  return status;
	}

        sol.resize(Linear_Eqn.a[0].size(), 1);

        for(int i(pivots.size() - 1); i >=0; i--)
        {
            int p_row = pivots[i].row;
            int p_col = pivots[i].column;
            sol[p_col] = Linear_Eqn.b[p_row];
            for(int j(p_col + 1); j < Linear_Eqn.a[0].size(); j++)
            {
		// Free variables are set to zero
                sol[p_col] -= Linear_Eqn.a[p_row][j] * sol[j];
            }

	      sol[p_col] = sol[p_col] / Linear_Eqn.a[p_row][p_col];
        }
	}

	return status; // Indicates the execution of the back_substitution step.
    }

    vector<long double> Gaussian_PP()
    {
        if (size == 0)
        {
            return {}; 
	    // return {}; indicates "return an object of the function's return type initialized with an empty list-initializer"
	    // Reference at: https://stackoverflow.com/questions/39487065/what-does-return-statement-mean-in-c11
        }

        vector<long double> sol;
	vector<bool> used_row(Linear_Eqn.a[0].size(), false);
	vector<bool> used_col(Linear_Eqn.a[0].size(), false);
	vector<Position> pivots;

        for(int i(0); i < Linear_Eqn.a.size(); i++)
        {
            // partial pivot
	  Position pivot_element;
            int max_row = pivoting(pivot_element, used_row, used_col);
	    if (Linear_Eqn.a[max_row][pivot_element.column] == 0)
	    {
	      // Mark the pivot_element.column as used.
	      used_col[pivot_element.column] = true;

	      // enter next iteration
	      continue;
	    }

	    pivots.push_back(pivot_element);   
            Row_Exchange(max_row, pivot_element.row);
	    used_row[pivot_element.row] = true;
	    used_col[pivot_element.column] = true;

	    /* cout << "--------------------" << endl; */
	    /* for(int i(0); i < Linear_Eqn.a.size(); i++) */
	    /* { */
	    /*   for(int j(0); j < Linear_Eqn.a[0].size(); j++) */
	    /*   { */
		/* cout << Linear_Eqn.a[i][j] << " "; */
	    /*   } */
	    /*   cout << Linear_Eqn.b[i] << endl; */
	    /* } */
	    
	    /* cout << "--------------------" << endl; */

            // Clear out the entries below the pivot
            elimination(pivot_element.row, pivot_element.column);
        }

        int status = back_substitution(sol, pivots);
	/* cout << "solution status " << status << endl; */
        /* return make_pair(code, sol); */
	return sol;
    }

};


void PrintColumn(const Column &column) {
    int size = column.size();
    std::cout.precision(PRECISION);
    for (int row = 0; row < size; ++row)
        std::cout << column[row] << " ";
    std::cout << std::endl;
}

int main() {
    Gaussian_Elimination_New Linear_Eqn;
    Linear_Eqn.ReadEquation();
    vector<long double> solution = Linear_Eqn.Gaussian_PP();
    PrintColumn(solution);
    return 0;
}
