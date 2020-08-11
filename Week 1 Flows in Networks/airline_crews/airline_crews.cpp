#include <iostream>
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
    int num_of_flights;
    int num_of_crews;
    int source = 0;
    int sink;

  public:
    void Data_Read()
    {
     // Remember to write a few lines of code to remove the loop while reading in the data.
     // The vertex indices are numbered from 1 instead of zero.
      int num_left, num_right;
      cin >> num_left >> num_right;
      num_of_flights = num_left;
      num_of_crews = num_right;
      int vertex_count = num_left + num_right + 2;
      // index 0 indicates the source, index vertex_count-1 indicates the sink.
      // 1 to num_left represent the left vertices.
      // num_left + 1 to (num_left + num_right) represent the right vertices.
      source = 0;
      sink = num_left + num_right + 1;
      flow_matrix.resize(vertex_count, vector<int>(vertex_count, 0));
      flow_matrix_res.resize(vertex_count, vector<int>(vertex_count, 0));
      capacity_matrix.resize(vertex_count, vector<int>(vertex_count, 0)); 

      for (int i = 1; i <= num_left; ++i)
	for (int j = 0; j < num_right; ++j) 
	{
	  int bit;
	  cin >> bit;
	  capacity_matrix[i][num_left + j + 1] = bit;
	  flow_matrix_res[i][num_left + j + 1] = bit;
	}

      // Add edges with weight 1 from source 0 to the left vertices
      for(int i(1); i <= num_left; i++)
      {
	capacity_matrix[0][i] = 1;
	flow_matrix_res[0][i] = 1;
      }

      // Add edges with capacity 1 form right vertices to sink 
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

    void print_result()
    {
      Edmonds_Karp();
      vector<int> matched(flow_matrix.size(), 0);
      vector<int> matched_crews(num_of_flights, -1);
      for(int i(1); i <= num_of_flights; i++)
	for(int j(num_of_flights + 1); j <= (num_of_flights + num_of_crews); j++)
      {
	if(flow_matrix[i][j] == 1 && matched[j] == 0)
	{
	  matched[j] = 1;
	  matched_crews[i - 1] = j - num_of_flights;
	}
      }

      for(int i(0); i < matched_crews.size(); i++)
      {
	cout << matched_crews[i] << " ";
      }
    }
    
};


int main()
{
  Network net;
  net.Data_Read(); 
  net.print_result();
}
