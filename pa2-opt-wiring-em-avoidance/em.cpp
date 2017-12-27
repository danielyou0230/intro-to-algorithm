#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <vector>
#include <limits.h>

#define SORC 0
#define SINK 1
#define DISTANCE 0
#define PARENT   1 
#define SRC_SINK 0
#define SINK_SRC 1
using namespace std;

typedef struct 
{
	int x;
	int y;
	int capacity;
	int remain;
	int UID;
} node;

typedef struct {
	int f_dist; // forward
	int b_dist; // backward
	int f_cap;
	int b_cap;
	int src;
	int dest;
	// int type;
} edge_r;

typedef struct {
	int distance;
	int capacity;
} edge;

int N_SORC;
int N_SINK;
int n_node;

int** Graph;
int** dist_table;
edge** Edges;
edge_r** rGraph;

int** disMatrix;
vector<node> sources;
vector<node> sinks;

int distance(node*, node*);
void show_nodes(void);
void show_all_graph(void);
void show_distanceTable(void);
void update_edge (int src, int snk);
void update_rGraph (int src, int snk);
int find_bottleneck(vector<int> path);
vector<int> bellmanFord(void);

int main(int argc, char const *argv[])
{
	int tmp_cap;
	ifstream f_in;
	f_in.open(argv[1], ifstream::in);
	f_in >> n_node;
	cout << "Total number of nodes: " << n_node << endl;
	//
	for (int i = 0; i < n_node; ++i) {
		node tmp;
		f_in >> tmp.x;
		f_in >> tmp.y;
		f_in >> tmp_cap;
		tmp.capacity = abs(tmp_cap);
		tmp.remain = tmp.capacity;
		if (tmp_cap > 0) {
			sources.push_back(tmp);
			// cout << "[SOUR] ";
		}
		else {
			sinks.push_back(tmp);
		}
	}
	f_in.close();
	N_SORC = sources.size();
	N_SINK = sinks.size();
	//
	for (int i = 0; i < N_SORC; ++i)
		sources[i].UID = i;
	for (int i = 0; i < N_SINK; ++i) 
		sinks[i].UID = i + N_SORC;
	//
	disMatrix = new int*[N_SORC];
	for (int i = 0; i < N_SORC; ++i) {
		disMatrix[i] = new int[N_SINK];
		for (int j = 0; j < N_SINK; ++j) {
			disMatrix[i][j] = distance(&sources[i], &sinks[j]);
		}
	}
	//
	show_nodes();
	cout << endl;
	// initialise graph
	Graph = new int*[N_SORC];     // graph
	Edges = new edge*[N_SORC];     // capacity of edges
	rGraph = new edge_r*[N_SORC];
	for (int i = 0; i < N_SORC; ++i) {
		Graph[i] = new int[N_SINK]();
		Edges[i] = new edge[N_SINK]();
		rGraph[i] = new edge_r[N_SINK];
		for (int j = 0; j < N_SINK; ++j) {
			int tmp_dis = distance(&sources[i], &sinks[j]);
			int tmp_cap = min(sources[i].capacity, sinks[j].capacity);
			//
			Edges[i][j].distance = tmp_dis;
			Edges[i][j].capacity = tmp_cap;
			rGraph[i][j].src  = i;
			rGraph[i][j].dest = j;
		}
	}
	// show_all_graph();
	// init flow
	cout << "Initialising graph..." << endl;
	for (int i = 0; i < N_SORC; ++i) {
		for (int j = 0; j < N_SINK; ++j) { 
			// (1) source has no flow available
			if (sources[i].remain == 0) {
				// cout << "[SKIP] source " << i << " -> sink " << j << ": no remain flow-out." << endl;
				goto SKIP;
			}
			// (2) sink has no capacity to accept more
			if (sinks[j].remain == 0) {
				// cout << "[SKIP] source " << i << " -> sink " << j << ": no remain flow-in." << endl;
				goto SKIP;
			}
			// (3) the edge is fulfilled (might not happen)
			// if (Graph[i][j] == Edges[i][j]) {
			if (Graph[i][j] == Edges[i][j].capacity) {
				// cout << "[SKIP] source " << i << " -> sink " << j << ": edge reaches limit." << endl;
				goto SKIP;
			}
			//
			// if (sources[i].remain <= Edges[i][j] - Graph[i][j]) {
			if (sources[i].remain <= Edges[i][j].capacity - Graph[i][j]) {
				// sources[i].remain <= sinks[j].remain
				if (sources[i].remain <= sinks[j].remain) {
					// cout << "remain && edge enough" << endl;
					Graph[i][j] += sources[i].remain;
				}
				// sources[i].remain > sinks[j].remain
				else { 
					Graph[i][j] += sinks[j].remain;
				}
			}
			// else if (sources[i].remain > Edges[i][j] - Graph[i][j]) {
			else if (sources[i].remain > Edges[i][j].capacity - Graph[i][j]) {
				// sinks[j].remain <= Edges[i][j] - Graph[i][j]
				// if (sinks[j].remain <= Edges[i][j] - Graph[i][j]) {
				if (sinks[j].remain <= Edges[i][j].capacity - Graph[i][j]) {
					Graph[i][j] += sinks[j].remain;
				}
				// sinks[j].remain > Edges[i][j] - Graph[i][j]
				else { 
					// Graph[i][j] += Edges[i][j] - Graph[i][j];
					// Graph[i][j] = Edges[i][j];
					Graph[i][j] = Edges[i][j].capacity;
				}
			}
			else {
				cout << "Exception occurred, exiting..." << endl;
				return -1;
			}
			update_edge(i, j);
			//
			SKIP:
			update_rGraph(i, j);
		}
	}
	//
	cout << endl << "Graph initialised:" << endl;
	show_all_graph();
	//
	vector<int> negative_cycle = bellmanFord();
	int bottleneck = find_bottleneck(negative_cycle);
	return 0;
}

int distance(node* a, node* b) {
	int x = abs(a->x - b->x);
	int y = abs(a->y - b->y);
	// cout << x << " " << " " << y << " " << sqrt(pow(x, 2) + pow(y, 2)) << endl;
	return x + y;
}

void show_nodes(void) {
	cout << "[SOURCES]" << endl;
	for (int i = 0; i < N_SORC; ++i) {
		cout << "Node #" << i + 1 << ": (";
		cout << setw(3) << sources[i].x << ", ";
		cout << setw(3) << sources[i].y << ")" ;
		cout << " capacity = " << setw(3) << sources[i].capacity;
		// cout << ", distance: ";
		for (int j = 0; j < sinks.size(); ++j) {
			cout << setw(3) << disMatrix[i][j] << " | ";
		}
		cout << endl;
	}
	cout << "[SINKS]" << endl;
	for (int i = 0; i < N_SINK; ++i) {
		cout << "Node #" << i + 1 << ": (";
		cout << setw(3) << sinks[i].x << ", ";
		cout << setw(3) << sinks[i].y << ")" ;
		cout << " capacity = " << setw(3) << sinks[i].capacity << endl;
	}
}

void show_all_graph(void) {
	cout << endl;
	cout << " TYPE | Node# | UID | Coordinate | Flow | Remain |" << endl;
	cout << "==================================================" << endl;
	for (int i = 0; i < N_SORC; ++i) {
		cout << " SORC | " << setw(5) << i << " | " << setw(3) << sources[i].UID << " | ";
		cout << "(" << setw(3) << sources[i].x << ", " << setw(3) << sources[i].y << ")";
		cout << " | " << setw(4) << sources[i].capacity;
		cout << " |  " << setw(5) << sources[i].remain << " | " << endl;
	}
	// cout << endl;
	for (int i = 0; i < N_SINK; ++i) {
		cout << " SINK | " << setw(5) << i << " | " << setw(3) << sinks[i].UID << " | ";
		cout << "(" << setw(3) << sinks[i].x << ", " << setw(3) << sinks[i].y << ")";
		cout << " | " << setw(4) << sinks[i].capacity;
		cout << " |  " << setw(5) << sinks[i].remain << " |" << endl;
	}
	cout << endl;
	cout << "Residual Graph" << endl;
	cout << "    Node     Coordinate   UID         Node    Coordinate   UID  Capacity   Length"<< endl;
	for (int i = 0; i < N_SORC; ++i) {
		for (int j = 0; j < N_SINK; ++j) {
			cout << " SORC #" << setw(3) << i << "   ";
			cout << "(" << setw(3) << sources[i].x << ", " << setw(3) << sources[i].y << ")   ";
			cout << setw(3) << sources[i].UID << "  ->  SINK #" << setw(3) << j;
			cout << "  (" << setw(3) << sinks[j].x << ", " << setw(3) << sinks[j].y << ")   ";
			cout << setw(3) << sinks[j].UID << "     " ;
			cout << setw(5) << rGraph[i][j].f_cap << "      ";
			cout << setw(3) << rGraph[i][j].f_dist << " " << endl;
			//
			cout << " SORC #" << setw(3) << i << "   ";
			cout << "(" << setw(3) << sources[i].x << ", " << setw(3) << sources[i].y << ")   ";
			cout << setw(3) << sources[i].UID << "  <-  SINK #" << setw(3) << j;
			cout << "  (" << setw(3) << sinks[j].x << ", " << setw(3) << sinks[j].y << ")   ";
			cout << setw(3) << sinks[j].UID << "     " ;
			cout << setw(5) << rGraph[i][j].b_cap << "      ";
			cout << setw(3) << rGraph[i][j].b_dist << " ";
			cout << endl;
		}
	}
	cout << endl;
}

void show_distanceTable(void) {
	cout << "Distance Table:" << endl;
	cout << "     Node    |  Coordinate  |  UID  | PARENT | DISTANCE |" << endl;
	int z = 0;
	for (int i = 0; i < n_node; ++i) {
		if (i < N_SORC) {
			cout << " source #" << setw(3) << i << " | (";
			cout << setw(4) << sources[i].x << ", ";
			cout << setw(4) << sources[i].y << ") | ";
		}
		else {
			int j = i - N_SORC;
			cout << " sink   #" << setw(3) << j << " | (";
			cout << setw(4) << sinks[j].x << ", ";
			cout << setw(4) << sinks[j].y << ") | ";
		}
		//
		cout << setw(5) << z << " |  ";
		cout << setw(5) << dist_table[i][PARENT] << " | " << setw(8) << dist_table[i][DISTANCE] << " |"<< endl;
		// if (n_node % N_SINK == N_SINK - 1)
		// 	cout << endl;
		++z;
	}
	cout << endl;
}

void update_edge (int src, int snk) {
	int flow = 0;
	// cout << "Updating graph" << endl;
	cout << endl;
	// collect all outbound flows (KCL)
	for (int k = 0; k < N_SINK; ++k)
		flow += Graph[src][k];
	// update remain
	if (flow > sources[src].capacity) {
		cout << "[UPDATE] Flow out of sources " << src << " is invalid: ";
	}
	else {
		sources[src].remain = sources[src].capacity - flow;
		cout << "[UPDATE] Flow out of sources " << src << " succeed: ";
	}
	cout << flow  << "/" << sources[src].capacity << endl;
	// collect all inbound flows (KCL)
	flow = 0;
	for (int k = 0; k < N_SORC; ++k)
		flow += Graph[k][snk];
	// update remain
	if (flow > sinks[snk].capacity) {
		cout << "[UPDATE] Flow into sink " << snk << " is invalid: ";
	}
	else {
		sinks[snk].remain = sinks[snk].capacity - flow;
		cout << "[UPDATE] Flow into sink " << snk << " succeed: ";
	}
	cout << flow  << "/" << sinks[snk].capacity << endl;
}

void update_rGraph (int src, int snk) {
	rGraph[src][snk].f_cap = Edges[src][snk].capacity - Graph[src][snk];
	rGraph[src][snk].b_cap = Graph[src][snk];
	rGraph[src][snk].f_dist = rGraph[src][snk].f_cap == 0 ? 0 :  disMatrix[src][snk];
	rGraph[src][snk].b_dist = rGraph[src][snk].b_cap == 0 ? 0 : -disMatrix[src][snk];
}

int find_bottleneck(vector<int> path) {
	int capacity = 0;
	int btk_u, btk_v;
	for (int i = 0; i < path.size(); ++i) {
		// Edge (u, v)
		int u = path[(i + 1) % path.size()]; // uid
		int v = path[i];
		//cout << u << " -> " << v << ": ";
		//
		if (u > N_SORC) { // u -> v (sink -> src)
			int _u = u - N_SORC;
			if (capacity < rGraph[v][_u].b_cap) {
				capacity = rGraph[v][_u].b_cap;
				btk_u = u;
				btk_v = v;
			}
		}
		if (v > N_SORC) { // u -> v (src -> sink) 
			int _v = v - N_SORC;
			if (capacity < rGraph[_v][u].f_cap) {
				capacity = rGraph[_v][u].f_cap;
				btk_u = u;
				btk_v = v;
			}
		}
	}
	cout << "Bottleneck at " << btk_u << " -> " << btk_v << " | Capacity = " << setw(3) << capacity << endl;
	return capacity;
}

vector<int> bellmanFord(void) {
	dist_table = new int*[n_node];
	vector<int> negative_cycle;
	// dist_table = [sources] [sinks]
	// for v in V:
	// 	v.distance = infinity
	// 	v.p = None
	for (int i = 0; i < n_node; ++i) {
		dist_table[i] = new int[2];
		if (i < N_SORC)
			dist_table[i][DISTANCE] = 0;
		else
			dist_table[i][DISTANCE] = INT_MAX / 2;
		dist_table[i][PARENT]   = -1;
	}
	/*/ Step 2: Relax all edges |V| - 1 times. A simple shortest path from src
	// to any other vertex can have at-most |V| - 1 edges
	for (i = 1; i <= V-1; i++) {
		for (j = 0; j < E; j++) {
			int u = graph->edge[j].src;
			int v = graph->edge[j].dest;
			int weight = graph->edge[j].weight;
			if (dist[u] + weight < dist[v])
				dist[v] = dist[u] + weight;
		}
	}
	*/
	// cout << "Initialised Distance Table" << endl;
	// show_distanceTable();
	// 
	bool update;
	for (int i = 0; i < n_node; ++i) { // |V| - 1
		// cout << "Iteration #" << setw(5) << i + 1 << "/"<< setw(5) << n_node << endl;
		update = false;
		// check all node
		// cout << "Total  = " << n_node << endl;
		for (int itr_node = 0; itr_node < n_node; ++itr_node) {
			// check all outbound edges of the node
			if (itr_node < N_SORC) { // for sources
				// first half indices belong to sources
				int sorc_idx = itr_node; 
				//
				for (int itr = 0; itr < N_SINK; ++itr) {
					int weight = rGraph[sorc_idx][itr].f_dist;
					int idx_sink = N_SORC + itr;
					//
					if (weight == 0) {
						// cout << itr_node << " -> " << idx_sink << " " << weight << endl;
						continue;
					}
					// update all connected sinks
					if (dist_table[itr_node][DISTANCE] + weight < dist_table[idx_sink][DISTANCE]) {
						dist_table[idx_sink][DISTANCE] = dist_table[itr_node][DISTANCE] + weight;
						// dist_table[idx_sink][PARENT] = idx_sink;
						// cout << "[sorc " << itr_node << " PARENT = " << dist_table[idx_sink][PARENT] << endl;
						dist_table[idx_sink][PARENT] = itr_node;
						// cout << "sink " << idx_sink << " PARENT = " << dist_table[idx_sink][PARENT] << endl;
						update = true;
					}
				}
			}
			//
			else { // for sinks
				// last half indices belong to sinks
				int sink_idx = itr_node - N_SORC; 
				// 
				for (int itr = 0; itr < N_SORC; ++itr) {
					int weight = rGraph[itr][sink_idx].b_dist;
					int idx_sorc = itr;
					//
					if (weight == 0) {
						// cout << itr << " <- " << itr_node << " " << weight << endl;
						continue;
					}
					// update all connected sources
					if (dist_table[itr_node][DISTANCE] + weight < dist_table[idx_sorc][DISTANCE]) {
						dist_table[idx_sorc][DISTANCE] = dist_table[itr_node][DISTANCE] + weight;
						// dist_table[idx_sorc][PARENT] = idx_sorc;
						// cout << "[sink " << itr_node << " PARENT = " << dist_table[idx_sorc][PARENT] << endl;
						dist_table[idx_sorc][PARENT] = itr_node;
						// cout << "sorc " << idx_sorc << " PARENT = " << dist_table[idx_sorc][PARENT] << endl;
						update = true;
					}
				}
			}
			// if (update == false){
			// 	cout << "Early stopped due to no update on edges." << endl;
			//	break;
			//}
		}
		// if (update == false)
		//	break;
		//cout << endl;
	}
	show_distanceTable();
	// for (u, v) in E:
	// 	if v.distance > u.distance + weight(u, v):
	// 		print "A negative weight cycle exists"
	bool found = false;
	for (int i = 0; i < N_SORC; ++i) {
		for (int j = 0; j < N_SINK; ++j) {
			int weight, u, v, index;
			// convert source/sink index to dist_table's index
			u = rGraph[i][j].src;
			v = rGraph[i][j].dest + N_SORC;
			// forward
			weight = rGraph[i][j].f_dist;
			if (weight != 0) {
				if (dist_table[v][DISTANCE] > dist_table[u][DISTANCE] + weight) {
					// cout << "Negative cycle found on src#" << rGraph[i][j].src << " to snk#" << rGraph[i][j].dest << endl;
					cout << "Negative cycle found" << endl;
					for (index = v;;) {
						cout << index << " <- ";
						negative_cycle.push_back(index);
						index = dist_table[index][PARENT];
						if (index == v) {
							found = true;
							break;
						}
					}
					cout << index << endl;
				}
			} // end checking forward edge
			// backward
			weight = rGraph[i][j].b_dist;
			if (weight != 0) {
				// convert source/sink index to dist_table's index
				if (dist_table[u][DISTANCE] > dist_table[v][DISTANCE] + weight) {
					// cout << "Negative cycle found on snk#" << rGraph[i][j].dest << " to src#" << rGraph[i][j].src << endl;
					cout << "Negative cycle found" << endl;
					for (index = v;;) {
						cout << index << " <- ";
						negative_cycle.push_back(index);
						index = dist_table[index][PARENT];
						if (index == v) {
							found = true;
							break;
						}
					} // end printing
					cout << index << endl;
				}
			} // end checking backward edge
		} // end itr sink
		if (found)
			break;
	} // end itr sorc
	//
	if (found) {
		/*
		for (int i = 0; i < negative_cycle.size(); ++i) 
			cout << negative_cycle[i] << " -> ";
		cout << endl;
		*/
	}
	return negative_cycle;
}