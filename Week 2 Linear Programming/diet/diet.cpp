#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <utility>  // pair
#include <limits>
#include <bitset>
#include <numeric>
using namespace std;

const long double EPS = 1e-6;
const int PRECISION = 20;
const long double INF = 1.0e+9;      
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
};


/**
    Solve the System of Linear Equations.
    Support arbitrary square matrix of Coefficients. 
**/
class Gaussian_Elimination
{

/*

*/
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

	// Extreme Case: No pivot at all
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
	    status = 2;
	    sol.push_back(INF); // The whole R^n space is feasible.
	  }
	}

	if (pivots.size() > 0)
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

        vector<long double> sol;  // Initally, set sol to empty.
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

class LP
{
  private:
    Matrix A;
    vector<long double> b;
    vector<long double> c;

  public:
    void read_data()
    {
      int n, m;   // n restrictions and m variables.
      cin >> n >> m;
      A.resize(n, vector<long double>(m));
      for (int i = 0; i < n; i++) {
	for (int j = 0; j < m; j++) {
	  cin >> A[i][j];
	}
      }
      b.resize(n);
      for (int i = 0; i < n; i++) {
	cin >> b[i];
      }
      c.resize(m);
      for (int i = 0; i < m; i++) {
	cin >> c[i];
      }
    }

    void prepare()
    {
      int n = A.size();  // number of constraints
      int m = A[0].size(); // number of variables
      
      // Add constraints to check for infinity solution.
      A.emplace_back(vector<long double>(m, 1));
      b.emplace_back(INF);
      ++n;

      // Add Positive Amount Constraints
      for (int i(0); i < m; i++)
      {
	vector<long double> vec(m, 0.0);
	vec[i] = -1;  // Ensure to get -x <= 0;
	A.emplace_back(std::move(vec));
	b.emplace_back(0.0);
	n++;
      }
    }

    
    void get_k_subsets(vector<int> &superset, int n, int k, int next_idx,
	vector<int> set, vector<vector<int>> &k_subsets, int i)
    {
      // base case: already k elements in set
      if (next_idx == k)
      {
	k_subsets.push_back(set);
	return;
      }

      if (i >= n)
      {
	return;
      }

      // Size of subset does not reach k
      set[next_idx] = superset[i];
      // i-th element is included 
      get_k_subsets(superset, n, k, next_idx + 1, set, k_subsets, i + 1);
      
      // i-th element is excluded
      // set[next_idx] need to be replaced.
      get_k_subsets(superset, n, k, next_idx, set, k_subsets, i + 1);
      
    }
    
    vector<Column> solve_all_equations()
    {

      vector<Column> all_sols;
      int n = A.size(); // number of the constraints;
      int m = A[0].size(); // number of variables; 
      // n >= m after preparing.
    
      vector<int> superset;
      for (int i(0); i < n; i++)
      {
	superset.push_back(i);
      }
      vector<int> set(m);

      vector<vector<int>> k_subsets;
      get_k_subsets(superset, n, m, 0, set, k_subsets, 0);
      /* cout << "size of k_subsets" << k_subsets.size(); */

      for (const auto& sub : k_subsets) {

        Matrix mat;
        Column col;

        for (auto j : sub) {
            mat.push_back(A[j]);
            col.push_back(b[j]);
        }

	Gaussian_Elimination Linear_Solver;
	Linear_Solver.Linear_Eqn.a = mat;
	Linear_Solver.Linear_Eqn.b = col;
	Linear_Solver.size = mat.size();
	auto sol = Linear_Solver.Gaussian_PP();

        if (!sol.empty()) {
	  // Solution Exists.
	  all_sols.push_back(sol);
        }
      }
      return all_sols;
    }
   
    pair<int, vector<long double>>
    solve_diet_problem()
    {
        prepare();
	/* cout << "Preprocessed A and b " << endl; */
	/* print_mat(A); */
	/* print_vec(b); */
	int n = A.size();
	int m = A[0].size();
	vector<Column> solutions = solve_all_equations();
	/* cout << "n " << n <<endl; */
	/* cout << "m " << m << endl; */
	/* cout << "solution size " << solutions.size() << endl; */

	if (solutions.size() == 0) {
	    return { -1, {} }; // no solution
	}

	int sol_index = -1;
	long double largest_pleasure = -(std::numeric_limits<long double>::max() / 2);

	// check solutions
	for (int i = 0; i < solutions.size(); ++i) 
	{

	    auto& sol = solutions[i];
	    /* cout << "One Vertex Solution" << endl; */
	    /* print_vec(sol); */
	    bool satisfied{ true };

	    for (int j = 0; j < n; ++j) 
	    {

		long double sum{ 0.0 };

		for (int k = 0; k < m; ++k) {
		    sum += A[j][k] * sol[k];
		}

		if ((sum - EPS) >  b[j]) {
		    satisfied = false;
		    break;
		}
	    }

	    long double pleasure{ 0.0 };
	    for (int k = 0; k < m; ++k) {
		pleasure += sol[k] * c[k];
	    }

	    if (satisfied && pleasure >= largest_pleasure) {
		largest_pleasure = pleasure;
		sol_index = i;
	    }
	}

	if (sol_index == -1) {
	    return { -1, {} }; // no solution
	}

	/* cout << "Optimal Solution " << sol_index << endl; */
	/* print_vec(solutions[sol_index]); */
	auto& sol = solutions[sol_index];  // solution to the linear programming
	long double  sum_of_vars = 0;
	for (int k(0); k < sol.size(); k++)
	{
	  sum_of_vars += sol[k];
	}
	if (sum_of_vars + EPS >= INF) {
	    return { 1, {} }; // infinity
	}

	return { 0, sol};

    }

    void print_mat(vector<vector<long double>> mat)
    {
      for(int i(0); i < mat.size(); i++)
      {
	for(int j(0); j < mat[0].size(); j++)
	{
	  cout << mat[i][j] << " ";
	}
	cout << endl;
      }
    }

    void print_vec(vector<long double> vec)
    {
      for (int k(0); k < vec.size(); k++)
      {
	cout << vec[k] << " ";
      }
      cout << endl;
    }
};

int main(){

  /* int n, m;   // n restrictions and m variables. */
  /* cin >> n >> m; */
  /* Matrix A(n, vector<long double>(m)); */
  /* for (int i = 0; i < n; i++) { */
  /*   for (int j = 0; j < m; j++) { */
  /*     cin >> A[i][j]; */
  /*   } */
  /* } */
  /* vector<long double> b(n); */
  /* for (int i = 0; i < n; i++) { */
  /*   cin >> b[i]; */
  /* } */
  /* vector<long double> c(m); */
  /* for (int i = 0; i < m; i++) { */
  /*   cin >> c[i]; */
  /* } */

  LP linear_solver;
  linear_solver.read_data();
  pair<int, vector<long double>> ans = linear_solver.solve_diet_problem();
  int m = ans.second.size();

  switch (ans.first) {
    case -1: 
      printf("No solution\n");
      break;
    case 0: 
      printf("Bounded solution\n");
      for (int i = 0; i < m; i++) {
        printf("%.18Lf%c", ans.second[i], " \n"[i + 1 == m]);
      }
      break;
    case 1:
      printf("Infinity\n");
      break;      
  }
  return 0;
}
