#include <iostream>
#include <cstdio>  // freopen
#include <fstream> // std::ifstream
#include <sstream> // std::istringstream

#include <vector>
#include <queue>
#include <stack>
#include <limits>
using namespace std;
const int int_max = std::numeric_limits<int>::max();

/*
   To index the vertex and edge by the index of the vertex quickly, we adopt the capacity  matrix as the approach to represent Network.
   Zero capacity on the edge means no link between vertices.
   Optimization: No need to compute a new residual Network at each iteration.
*/

class Network
{
  private:
    vector<vector<int>> capacity_matrix;
    vector<vector<int>> flow_matrix;
    vector<vector<int>> flow_matrix_res;
    vector<vector<int>> stock_data;
    int num_of_flights;
    int num_of_crews;
    int source = 0;
    int sink;

  public:
    void Data_Read()
    {
      /*
	 stock_data[i][j] means the price of i-th stock at j-th time point
	 Make sure to do unit test for Data_Read() procedure.
      */
      int num_of_stocks, num_of_points;
      cin >> num_of_stocks >> num_of_points;
      stock_data.resize(num_of_stocks, vector<int>(num_of_points, 0));
      
      for(int i(0); i < num_of_stocks; i++)
	for(int j(0); j < num_of_points; j++)
	{
	  cin >> stock_data[i][j];
	}

    }

    void graph_construct()
    {
      /*
	 Reference at http://mradwan.github.io/algorithms/2014/05/02/flows-cuts-and-matchings/
	 Create a DAG from the given charts, where each vertex corresponds to a chart, 
	 vertex V_i has an edge going to V_j if each kth point in the ith chart is less than the kth point in the jth chart, 
	 the answer then is the minimum path cover in that DAG.
      */
      num_of_flights = stock_data.size();
      num_of_crews = stock_data.size();
      int num_left = num_of_flights;
      int num_right = num_of_crews;

      int vertex_count = num_left + num_right + 2;
      // index 0 indicates the source, index vertex_count-1 indicates the sink.
      // 1 to num_left represent the left vertices.
      // num_left + 1 to (num_left + num_right) represent the right vertices.
      source = 0;
      sink = num_left + num_right + 1;
      flow_matrix.resize(vertex_count, vector<int>(vertex_count, 0));
      flow_matrix_res.resize(vertex_count, vector<int>(vertex_count, 0));
      capacity_matrix.resize(vertex_count, vector<int>(vertex_count, 0)); 

      // Connect left vertices to the acceptable right vertices with an directed edge of capacity 1
      for (int i = 1; i <= num_left; ++i)
	for (int j = 1; j <= num_right; ++j) 
	{
	  if (i == j)
	  {
	    // avoid loop.
	    continue;
	  }
	  
	  // check whether each k-th point in the left charts is strictly less than each k-th point in the right charts.
	  bool is_less = true;
	  for(int k(0); k < stock_data[0].size(); k++)
	  {
	    if (stock_data[i-1][k] >= stock_data[j-1][k])
	    {
	      is_less = false;
	      break;
	    }
	  }

	  if (is_less)
	  {
	    // Right vertex is accepted by the left vertex
	    capacity_matrix[i][num_left + j] = 1;
	    flow_matrix_res[i][num_left + j] = 1;
	  }
	}

      // Connect source 0 to each left vertex with a directed edge of capacity 1.
      for(int i(1); i <= num_left; i++)
      {
	capacity_matrix[0][i] = 1;
	flow_matrix_res[0][i] = 1;
      }

      // Connect each of the right vertices to sink with a directed edge of capacity 1 .
      for(int i(num_left +1); i <= (num_left + num_right); i++)
      {
	capacity_matrix[i][sink] = 1; 
	flow_matrix_res[i][sink] = 1;
      } 
    }

    bool BFS(int s, int t, vector<vector<int>> &flow_mat, vector<int> &prev)
    {
      /* Find a path from s to t.
	 Input:
	      s: the source; 
	      t: the sink.
	      flow_mat: The matrix of flows of the residual network.
	      prev: Record the parents of each vertex in the resulting shortest path tree.
	 Output:
		true if there is a path from s to t, false otherwise.
      */
      
      int num_of_vertices = flow_mat.size();
      vector<int> dist(num_of_vertices, int_max);

      dist[s] = 0;
      queue<int> queue_vertex;
      queue_vertex.push(s);

      while(!queue_vertex.empty())
      {
	int u = queue_vertex.front();
	queue_vertex.pop();

	for(int i(0); i < num_of_vertices; i++)
	{
	  if (flow_mat[u][i]> 0)
	  {
	    // there is edge from u to i.
	    // i is adjacent to u.
	    int v = i;
	    if (dist[v] == int_max)
	    {
	      // v has not been discovered
	      // enqueue v
	      queue_vertex.push(v);
	      dist[v] = dist[u] + 1;
	      prev[v] = u;
	    }
	  }
	}
      }

      return dist[t] != int_max;
    }
   
    void update_flow_matrix(int &total_flow, vector<int> &prev)
    {
      /*
	 Input:
	    max_flow(int): The total flow in the network.
	    prev: shortest path tree represented in the form of parent array. source is the root.
	 Remarks:
		Given current flow matrix and residual network, update the flow in the original network and
		residual network.
		And update the total flow;
      */

      // Find the min flow along the path
      int min_flow = int_max;
      for(int v(sink); v != source; v = prev[v])
      { 
	int u = prev[v]; // from u to v
	if (flow_matrix_res[u][v] < min_flow)
	{
	  min_flow = flow_matrix_res[u][v];
	}	
      }

      // update the flow in  the original network.
      for (int v(sink); v != source; v = prev[v])
      {
	int u = prev[v]; // from u to v

	// update the flow matrix of the residual network.
	// total sum of the flow between two vertices in the residual network equals the capacity.
	flow_matrix_res[u][v] -= min_flow;
	flow_matrix_res[v][u] += min_flow; 
	
	// update flow matrix of the original network.
	if (capacity_matrix[u][v] > 0)
	{
	  flow_matrix[u][v] += min_flow;
	}
	else
	{
	  flow_matrix[v][u] -= min_flow;
	}
      }

      // We add same flow min_flow to each edge in the selected path to ensure the conservation of flow.
      total_flow += min_flow;
    }
   
    int Edmonds_Karp()
    {
      /* Return the Maxflow
      */
      int total_flow = 0;
      int num_of_vertices = capacity_matrix.size();

      while(true)
      {
	/* vector<vector<int>> flow_matrix_res(num_of_vertices, vector<int>(num_of_vertices, 0)); */
	/* compute_residual_net(flow_matrix_res); */
	vector<int> prev(num_of_vertices, -1);
	bool is_path = BFS(source, sink, flow_matrix_res, prev);

	if(!is_path)
	{
	  return total_flow;
	}

	update_flow_matrix(total_flow, prev);
      }
    }

    int print_result()
    {
      graph_construct(); //  construct augmented bipartite graph.
      Edmonds_Karp();    //  Find the maxflow
      vector<int> matched(flow_matrix.size(), 0);
      int min_num_charts = num_of_flights;
      for(int i(1); i <= num_of_flights; i++)
	for(int j(num_of_flights + 1); j <= (num_of_flights + num_of_crews); j++)
      {
	if(flow_matrix[i][j] == 1 && matched[j] == 0)
	{
	  matched[j] = 1;
	  min_num_charts--;  // One match decreases the total number of charts by 1.
	}
      }

      return min_num_charts;
      /* cout << min_num_charts << endl; */
    } 
};


int main()
{
  // For Test
  /*
    Reference:
    convert int to string: https://stackoverflow.com/questions/5590381/easiest-way-to-convert-int-to-string-in-c
    convert string to file name: https://stackoverflow.com/questions/36824225/use-string-or-array-of-char-for-filename
    Redirect Standard Input and Output: https://stackoverflow.com/questions/5257509/freopen-equivalent-for-c-streams
    C++ Files and Streams: https://www.tutorialspoint.com/cplusplus/cpp_files_streams.htm#:~:text=A%20file%20must%20be%20opened,file%20for%20reading%20purpose%20only.
    getline: https://www.geeksforgeeks.org/getline-string-c/
    Read from file line by line: https://stackoverflow.com/questions/7868936/read-file-line-by-line-using-ifstream-in-c
    https://stackoverflow.com/questions/16777451/c-how-is-istream-is-converted-to-bool-inside-a-conditional-expression
  */

  int num_of_test;
  for (int i(1); i <= num_of_test; i++)
  {
    string test_file;
    string result_file;
    int min_chart;
    if(i <10)
    {
      test_file = "tests/0" + to_string(i);
    }
    else
    {
      test_file = "tests/" + to_string(i);
    }
    result_file = test_file + ".a";

    freopen(test_file.c_str(), "r", stdin);  

    std::ifstream infile(result_file.c_str());
    std::string line;
    while(std::getline(infile, line))
    {
      std::istringstream iss(line);
      if (!(iss >> min_chart))
      {
	break; //error
      }
    }
    Network net;
    net.Data_Read(); 
    int res = net.print_result();
    if (res != min_chart)
    {
      cout << "Fails on test case " << i << endl;
    }
    else
    {
      cout << "Succeed!" << endl;
    }
  }


  // For submission
  /* Network net; */
  /* net.Data_Read(); */ 
  /* cout << net.print_result() << endl; */
}
