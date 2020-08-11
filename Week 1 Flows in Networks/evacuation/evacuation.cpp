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
    int source = 0;
    int sink;

  public:
    void Data_Read()
    {
     // Remember to write a few lines of code to remove the loop while reading in the data.
     // The vertex indices are numbered from 1 instead of zero.
      int vertex_count, edge_count;
      std::cin >> vertex_count >> edge_count;
      sink = vertex_count - 1;
      flow_matrix.resize(vertex_count, vector<int>(vertex_count, 0));
      flow_matrix_res.resize(vertex_count, vector<int>(vertex_count, 0));
      capacity_matrix.resize(vertex_count, vector<int>(vertex_count, 0));
      for (int i = 0; i < edge_count; ++i) {
	  int u, v, capacity;
	  std::cin >> u >> v >> capacity;
	  if (u == v)
	  {
	  // To avoid loop in the Network, we ignore the road from a city to itself.
	    continue;
	  }
	  capacity_matrix[u-1][v-1] += capacity;  // merge same edge into one by adding capacities
	  flow_matrix_res[u-1][v-1] += capacity;
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
};


int main()
{
  Network net;
  net.Data_Read(); 
  cout << net.Edmonds_Karp() << endl;
}
