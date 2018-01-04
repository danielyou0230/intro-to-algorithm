#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>

#define N 0
#define E 1
#define S 2
#define W 3

#define DIRECTIONS 12
#define CHAR_TO_DIR(c) 				\
(									\
	(c == 'N')?(0):					\
	(								\
		(c == 'E')?(1):				\
		(							\
			(c == 'S')?(2):			\
			(						\
				(c == 'W')?(3):(4)	\
			)						\
		)							\
	)								\
)

#define DIR_TO_CHAR(c) 				\
(									\
	(c == 0)?('N'):					\
	(								\
		(c == 1)?('E'):				\
		(							\
			(c == 2)?('S'):			\
			(						\
				(c == 3)?('W'):('X')\
			)						\
		)							\
	)								\
)

using namespace std;

int conflict_table[16][16] = 
	//   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
{	//  NN  NE  NS  NW  EN  EE  ES  EW  SN  SE  SS  SW  WN  WE  WS  WW
	  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // NN 1
	  { -1, -1, -1, -1,  1, -1,  0,  0,  0,  0, -1,  0,  0,  0,  0, -1 }, // NE 2
	  { -1, -1, -1, -1,  1, -1, -1,  0,  1,  1, -1,  0,  0,  0,  0, -1 }, // NS 3
	  { -1, -1, -1, -1,  1, -1,  0,  0,  1,  1, -1,  0,  1,  1,  1, -1 }, // NW 4
	  { -1,  1,  1,  1, -1, -1, -1, -1,  0,  1, -1,  0,  0,  1,  1, -1 }, // EN 5
	  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // EE
	  { -1,  0,  0,  0, -1, -1, -1, -1,  0,  1, -1,  0,  0,  0,  0, -1 }, // ES
	  { -1,  0,  0,  0, -1, -1, -1, -1,  0,  1, -1,  0,  0,  1,  1, -1 }, // EW
	  { -1,  0,  1,  1,  0, -1,  0,  0, -1, -1, -1, -1,  0,  0,  1, -1 }, // SN
	  { -1,  0,  1,  1,  1, -1,  1,  1, -1, -1, -1, -1,  0,  0,  1, -1 }, // SE
	  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // SS
	  { -1,  0,  0,  0,  0, -1,  0,  0, -1, -1, -1, -1,  0,  0,  1, -1 }, // SW
	  { -1,  0,  0,  1,  0, -1,  0,  0,  0,  0, -1,  0, -1, -1, -1, -1 }, // WN
	  { -1,  0,  0,  1,  1, -1,  0,  1,  0,  0, -1,  0, -1, -1, -1, -1 }, // WE
	  { -1,  0,  0,  1,  1, -1,  0,  1,  1,  1, -1,  1, -1, -1, -1, -1 }, // WS
	  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // WW
};

typedef struct{
	int op; 
	int time;
	int sched;
} car;

inline int convert_index(int come_from, int go_to) {
	if (go_to == '0')
		return -1;
	else
		return come_from * 4 + CHAR_TO_DIR(go_to);
}

inline void safe_ptr_add(int ptr, int sched_limit) {
	if (ptr != sched_limit - 1)
		++ptr;
}

int main(int argc, char const *argv[])
{
	ifstream f_in(argv[1]);
	ofstream f_out(argv[2]);
	string buffer;
	string tmp;
	int come_from = 0;
	char from[4] = {'N', 'E', 'S', 'W'};
	vector<car> v[4];
	bool finish[4] = {};
	/*/
	for (int i = 0; i < 16; ++i) {
		for (int j = 0; j < 16; ++j) 
			cout << setw(2) << conflict_table[i][j] << " ";
		cout << endl;
	}
	/*/
	while(getline(f_in, buffer)){
		stringstream ss(buffer);
		ss >> tmp;
		int arrival = 0;
		int i = 0;
		while(ss >> tmp){
			car tmp_car;
			tmp_car.op = convert_index(come_from, tmp[1]);
			tmp_car.time = arrival;
			// if (tmp[1] == '0') {
			// 	// cout << "0" << tmp[1] << " " << setw(2) << CHAR_TO_DIR(tmp[1]) << " ( X)  ";
			// 	cout << "0" << tmp[1] << " " << " ( X)  ";
			// }
			// else {
			// 	// cout << from[come_from] << tmp[1] << " " << setw(2) << CHAR_TO_DIR(tmp[1]) << " (";
			// 	cout << from[come_from] << tmp[1] << " " << " (";\
			// 	cout << setw(2) << convert_index(come_from, tmp[1]) << ")  ";
			// }
			// cout << "t = " << setw(2) << tmp_car.time << " | " << setw(2) << tmp_car.op << " | ";
			++arrival;
			v[come_from].push_back(tmp_car);
			// cout << setw(3) << tmp_car.op << " ";
		}
		++come_from;
		// cout << endl;
	}
	/*for (int k = 0; k < 4; ++k) {
		for (int i = 0; i < v[k].size(); ++i) {
			cout << from[k] << tmp[1] << " " << " (";\
			cout << setw(2) << v[k][i].op << ")  ";
			cout << "t = " << setw(2) << v[k][i].time << " | ";
		}
		cout << "\n";
	}*/
	f_in.close();
	int ptr[4] = {};
	int sched_limit = v[0].size();
	int current_time = 0;
	cout << "Scheduling... ";
	while (1) {
		int N_op = (ptr[N] == sched_limit) ? -3 : v[N][ptr[N]].op;
		int E_op = (ptr[E] == sched_limit) ? -3 : v[E][ptr[E]].op;
		int S_op = (ptr[S] == sched_limit) ? -3 : v[S][ptr[S]].op;
		int W_op = (ptr[W] == sched_limit) ? -3 : v[W][ptr[W]].op;
		/*
		cout << "Scheduling for time " << setw(2) << current_time << " | ";
		cout << setw(2) << N_op << "  " << setw(2) << E_op << " " << setw(2) << S_op << " " << setw(2) << W_op << endl;
		cout << ptr[N] << " " << ptr[E] << " " << ptr[S] << " " << ptr[W] << endl;
		*/
		if (N_op > 15 || E_op > 15 || S_op > 15 || W_op > 15) {
			cout << ptr[N] << " " << ptr[E] << " " << ptr[S] << " " << ptr[W] << endl;
			cout << "Failed" << endl;
			return 0;
		}
		//
		if (conflict_table[N_op][E_op] == 1 && N_op != -1 && E_op != -1 && N_op != -3 && E_op != -3) {
			// cout << N_op << " & " << E_op << ": "<< conflict_table[N_op][E_op] <<  endl;
			// cout << "[TWO] Pass car in N" << ptr[N] << " & E" << ptr[E] << endl;
			v[N][ptr[N]].sched = current_time;
			v[E][ptr[E]].sched = current_time;
			++ptr[N];
			++ptr[E];
		}
		else if (conflict_table[N_op][S_op] == 1 && N_op != -1 && S_op != -1 && N_op != -3 && S_op != -3) {
			// cout << N_op << " & " << S_op << ": "<< conflict_table[N_op][S_op] <<  endl;
			// cout << "[TWO] Pass car in N" << ptr[N] << " & S" << ptr[S] << endl;
			v[N][ptr[N]].sched = current_time;
			v[S][ptr[S]].sched = current_time;
			++ptr[N];
			++ptr[S];
		}
		else if (conflict_table[N_op][W_op] == 1 && N_op != -1 && W_op != -1 && N_op != -3 && W_op != -3) {
			// cout << N_op << " & " << W_op << ": "<< conflict_table[N_op][W_op] <<  endl;
			// cout << "[TWO] Pass car in N" << ptr[N] << " & W" << ptr[W] << endl;
			v[N][ptr[N]].sched = current_time;
			v[W][ptr[W]].sched = current_time;
			++ptr[N];
			++ptr[W];
		}
		else if (conflict_table[E_op][S_op] == 1 && E_op != -1 && S_op != -1 && E_op != -3 && S_op != -3) {
			// cout << E_op << " & " << S_op << ": "<< conflict_table[E_op][S_op] <<  endl;
			// cout << "[TWO] Pass car in E" << ptr[E] << " & S" << ptr[S] << endl;
			v[E][ptr[E]].sched = current_time;
			v[S][ptr[S]].sched = current_time;
			++ptr[E];
			++ptr[S];
		}
		else if (conflict_table[E_op][W_op] == 1 && E_op != -1 && W_op != -1 && E_op != -3 && W_op != -3) {
			// cout << E_op << " & " << W_op << ": "<< conflict_table[E_op][W_op] <<  endl;
			// cout << "[TWO] Pass car in E" << ptr[E] << "& W" << ptr[W] << endl;
			v[E][ptr[E]].sched = current_time;
			v[W][ptr[W]].sched = current_time;
			++ptr[E];
			++ptr[W];
		}
		else if (conflict_table[S_op][W_op] == 1 && S_op != -1 && W_op != -1 && S_op != -3 && W_op != -3) {
			// cout << S_op << " & " << W_op << ": "<< conflict_table[S_op][W_op] <<  endl;
			// cout << "[TWO] Pass car in S" << ptr[S] << " & W" << ptr[W] << endl;
			v[S][ptr[S]].sched = current_time;
			v[W][ptr[W]].sched = current_time;
			++ptr[S];
			++ptr[W];
		}
		else {
			// cout << "Only one car can pass at time " << setw(2) << current_time << endl;
			// if only one car can pass
			int minimum = sched_limit * 100;
			minimum = (minimum > v[N][ptr[N]].time && N_op != -1 && N_op != -3) ? N : minimum;
			minimum = (minimum > v[E][ptr[E]].time && E_op != -1 && E_op != -3) ? E : minimum;
			minimum = (minimum > v[S][ptr[S]].time && S_op != -1 && S_op != -3) ? S : minimum;
			minimum = (minimum > v[W][ptr[W]].time && W_op != -1 && W_op != -3) ? W : minimum;
			if (minimum == -1) {
				cout << "Error" << endl;
				return -1;
			}
			// if the there exist a car not 00
			if (minimum <= 3) {
				// cout << "Pass Car #" << setw(2) << minimum << endl;
				// v[minimum][ptr[minimum]].sched = current_time;
				if (ptr[minimum] < sched_limit) {
					v[minimum][ptr[minimum]].sched = current_time;
					++ptr[minimum];
				}
			}
		}
		if (N_op == -1 && finish[N] == false) {
			v[N][ptr[N]].sched = current_time;
			++ptr[N];
		}
		if (E_op == -1 && finish[E] == false) {
			v[E][ptr[E]].sched = current_time;
			++ptr[E];
		}
		if (S_op == -1 && finish[S] == false) {
			v[S][ptr[S]].sched = current_time;
			++ptr[S];
		}
		if (W_op == -1 && finish[W] == false) {
			v[W][ptr[W]].sched = current_time;
			++ptr[W];
		}
		//
		if (ptr[N] == sched_limit)
			finish[N] = true;
		if (ptr[E] == sched_limit)
			finish[E] = true;
		if (ptr[S] == sched_limit)
			finish[S] = true;
		if (ptr[W] == sched_limit)
			finish[W] = true;
		//
		++current_time;
		if (finish[N] && finish[E] && finish[S] && finish[W])
			break;
		if (current_time > sched_limit * 4)
			break;
	}
	//
	cout << "Total rounds: " << current_time << endl;
	int show_idx[4] = {};
	cout << "Result" << endl;
	for (int i = 0; i < 4; ++i) {
		cout << from[i] << ": ";
		f_out << from[i] << ": ";
		for (int j = 0; j < current_time; ++j) {
			if (v[i][show_idx[i]].op == -1) {
				cout << 0 << 0 << " ";
				f_out << 0 << 0 << " ";
				++show_idx[i];
			}
			else if (v[i][show_idx[i]].sched == j) {
				// cout << DIR_TO_CHAR(v[i][show_idx[i]].op / 4) << DIR_TO_CHAR(v[i][show_idx[i]].op % 4) << " ";
				cout << 1 << DIR_TO_CHAR(v[i][show_idx[i]].op % 4) << " ";
				f_out << 1 << DIR_TO_CHAR(v[i][show_idx[i]].op % 4) << " ";
				//cout << " (" << setw(2) << v[i][show_idx[i]].op << ") ";
				++show_idx[i];
			}
			else {
				cout << 0 << 0 << " ";
				f_out << 0 << 0 << " ";
			}
		}
		cout << endl;
		f_out << "\n";
	}
	f_out.close();
	return 0;
}
