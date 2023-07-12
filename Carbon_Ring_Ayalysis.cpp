#pragma once
#include<map>
#include<vector>
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector> 
#include <string> 
#include <fstream> 
#include <sstream>
using namespace std;

class Solve_Circle {
	int VERTEX_SIZE;
	int **ADJMatrix;
	int *visitedFlag;

	int seed=0;
	int innerStep = 0;
	int recallVertex;
	int isRecall=0;
	vector<int> loop_vector;
public:
	map<int, int> c_c;
	~Solve_Circle() {
		delete ADJMatrix;
		delete visitedFlag;
	}

	void initialVisitedFlagArray()
	{
		for (int i = 0; i < VERTEX_SIZE; i++)
			visitedFlag[i] = 0;
	}
	void printVisitedVertex(int vertexID)
	{
		printf("visited: %d \n", vertexID);
	}
	void setVisitedFlag(int vertexID, int value)
	{
		visitedFlag[vertexID] = value;
	}
	int firstAdjacentVertex(int vertexID)
	{
		for (int i = 0; i < VERTEX_SIZE; ++i) {
			if (ADJMatrix[vertexID][i] == 1)
				return i;
		}
		return -1;
	}
	int nextAdjacentVertex(int vertexID, int nextVertexID)
	{
		for (int i = nextVertexID + 1; i < VERTEX_SIZE; ++i) {
			if (ADJMatrix[vertexID][i] == 1)
				return i;
		}
		return -1;
	}

	void printVector(vector<int> vct)
	{
		if (!vct.empty()) {
			for (auto v : vct) cout << v << "->";
			cout << vct[0] << '\n';
		}
	}

	void DFS(int vertex)
	{
		
		setVisitedFlag(vertex, 1);
		int nextVertex;
		loop_vector.push_back(vertex);
		nextVertex = firstAdjacentVertex(vertex);
		innerStep++;
		while (1) {
			if (innerStep > 7) {
				isRecall = 1;
				break;
			}
			if (nextVertex != -1) {
				if (visitedFlag[nextVertex] == 1 &&
					nextVertex == seed &&
					innerStep == 2) {
					nextVertex = nextAdjacentVertex(vertex, nextVertex);
					continue;
				}
				else if (visitedFlag[nextVertex] == 1 &&
					nextVertex == seed &&
					innerStep != 2) {
					c_c[innerStep]++;
					//printVector(loop_vector);
					nextVertex = nextAdjacentVertex(vertex, nextVertex);
					continue;
				}
				else if (visitedFlag[nextVertex] == 0) {
					DFS(nextVertex);
				}
				if (isRecall == 1) {
					innerStep--;
					recallVertex = nextVertex;
					nextVertex = nextAdjacentVertex(vertex, nextVertex);
					loop_vector.pop_back();
					setVisitedFlag(recallVertex, 0);
					isRecall = 0;
					continue;
				}
				nextVertex = nextAdjacentVertex(vertex, nextVertex);
			}
			else if (nextVertex == -1) {
				isRecall = 1;
				break;
			}
		}
	}

	void DFSTraverse()
	{
		initialVisitedFlagArray();
		for (seed = 0; seed < VERTEX_SIZE; ++seed) {
			for (int i = 0; i < VERTEX_SIZE; ++i)
				visitedFlag[i] = 0;

			if (visitedFlag[seed] == 0) {
				//cout << "\n-------------------the loop start and end with "
					//<< seed << "-------------------\n";
				loop_vector.clear();
				innerStep = 0;
				isRecall = 0;
				DFS(seed);
			}
		}
	}

	int circle_count(map<int, int>& m_i_i,string &trjstep_,ofstream &outfile_, int idex = 0)
	{
		int n = 0;
		int circles[5] = { 0,0,0,0,0 };
		cout << "------------Final results-------------" << endl;
		cout << trjstep_ << endl;
		for (auto m : m_i_i) {
			cout  << m.first << "-atom circle: " << m.second / m.first / 2 << endl;
			outfile_ << trjstep_<<" "<< m.first << " " << m.second / m.first / 2 << endl;
			n += m.second / m.first;
		}

		if (idex)
			return n;
		else
			return n / 2;
	}
	Solve_Circle(int** adjmatrix, int atoms);
};
Solve_Circle::Solve_Circle(int ** adjmatrix,int atoms) {
	VERTEX_SIZE=atoms;
	ADJMatrix= adjmatrix;
	visitedFlag=new int [VERTEX_SIZE];
}

class Read_lammpstrj {
	ifstream infile;
	string filename,timestep,temp;
	
	
	double* px, *py, *pz;
	char* element;
	double cutoff, box_xl = 0, box_xh = 0, box_yl = 0, box_yh = 0, box_zl = 0, box_zh = 0;
public:
	int** ADJmatrix;
	int* type;
	int atom_number,carbon_number;
	void readtrj();
	void bond();
	void carbon_bond();
	Read_lammpstrj(string filename_, string timestep_);
	~Read_lammpstrj() {

	}
	
};
Read_lammpstrj::Read_lammpstrj(string filename_,string timestep_) {
	filename = filename_;
	//cout << "Please input the timestep:";
	//cin >> timestep;
	timestep = timestep_;
	cutoff = 1.8;
	atom_number = 0;
	carbon_number = 0;
}
void Read_lammpstrj::readtrj() {
	infile.open(filename, ios::in);
	stringstream mysstream;
	while (getline(infile, temp))            
	{
		if (temp.compare("ITEM: TIMESTEP") == 0) {
			getline(infile, temp);
			if (temp.compare(timestep) == 0) {
				for (int i = 0; i < 8; i++) {
					getline(infile, temp);
					if (i == 1) {
						mysstream << temp;
						mysstream >> atom_number;
						mysstream.clear();
						mysstream.str("");
					}
					else if (i == 3) {
						mysstream << temp;
						mysstream >> box_xl >> box_xh;
						mysstream.clear();
						mysstream.str("");
					}
					else if (i == 4) {
						mysstream << temp;
						mysstream >> box_yl >> box_yh;
						mysstream.clear();
						mysstream.str("");
					}
					else if (i == 5) {
						mysstream << temp;
						mysstream >> box_zl >> box_zh;
						mysstream.clear();
						mysstream.str("");
					}
				}
				break;

			}
		}
	}
	px = new double[atom_number];
	py = new double[atom_number];
	pz = new double[atom_number];
	element = new char[atom_number];
	type = new int[atom_number];
	int atomid, takeup;

	
	for (int j = 0; j < atom_number; j++) {
		mysstream << temp;
		mysstream >> atomid >> type[j] >> px[j] >> py[j] >> pz[j] >> takeup >> takeup >> takeup >> takeup >> takeup;
		//mysstream >> element[j] >> px[j] >> py[j] >> pz[j];
		mysstream.clear();
		mysstream.str("");
		getline(infile, temp);
		//cout<< atomid <<" "<< type << " " << px[j] << " " << py[j] << " " << pz[j] <<endl;
	}
	infile.close();
}
void Read_lammpstrj::bond() {
	int bond_num = 0;
	double distance = 0, dx = 0, dy = 0, dz = 0;
	ADJmatrix = new int* [atom_number];
	for (int i = 0; i < atom_number; i++) {
		ADJmatrix[i] = new int[atom_number];
	}
	for (int i = 0; i < atom_number; i++)
		for (int j = 0; j < atom_number; j++)
			ADJmatrix[i][j] = 0;

	for (int i = 0; i < atom_number; i++) {
		for (int j = 0; j < atom_number; j++) {
			dx = px[i] - px[j];
			dy = py[i] - py[j];
			dz = pz[i] - pz[j];
			dx = abs(dx - box_xh * round(dx / box_xh));
			dy = abs(dy - box_yh * round(dy / box_yh));
			dz = abs(dz - box_zh * round(dz / box_zh));
			distance = dx * dx + dy * dy + dz * dz;
			if (distance < cutoff * cutoff&&i!=j) {
				bond_num++;
				ADJmatrix[i][j] = 1;
				
			}
		}
	}
}

void Read_lammpstrj::carbon_bond() {
	int bond_num = 0;
	double distance = 0, dx = 0, dy = 0, dz = 0;
	for (int i = 0; i < atom_number; i++) {
		if (type[i] == 2) {
			carbon_number++;
		}
	}
	double* Cpx, *Cpy, *Cpz;
	Cpx = new double[carbon_number];
	Cpy = new double[carbon_number];
	Cpz = new double[carbon_number];
	carbon_number = 0;
	for (int i = 0; i < atom_number; i++) {
		if (type[i] == 2) {
			Cpx[carbon_number] = px[i];
			Cpy[carbon_number] = py[i];
			Cpz[carbon_number] = pz[i];
			carbon_number++;
		}
	}
	ADJmatrix = new int* [carbon_number];
	for (int i = 0; i < carbon_number; i++) {
		ADJmatrix[i] = new int[carbon_number];
	}
	for (int i = 0; i < carbon_number; i++)
		for (int j = 0; j < carbon_number; j++)
			ADJmatrix[i][j] = 0;

	for (int i = 0; i < carbon_number; i++) {
		for (int j = 0; j < carbon_number; j++) {
			dx = Cpx[i] - Cpx[j];
			dy = Cpy[i] - Cpy[j];
			dz = Cpz[i] - Cpz[j];
			dx = abs(dx - box_xh * round(dx / box_xh));
			dy = abs(dy - box_yh * round(dy / box_yh));
			dz = abs(dz - box_zh * round(dz / box_zh));
			distance = dx * dx + dy * dy + dz * dz;
			if (distance < cutoff * cutoff && i != j) {
				bond_num++;
				ADJmatrix[i][j] = 1;

			}
		}
	}
}

int main(int argc, char* argv[])
{
	string trjstep;
	ofstream outfile;
	outfile.open("./Circle_results.txt");
	outfile << "timestep" << " " << "Circle atoms" << " " << "Number" << endl;
	int startstep=6000, intervalstep =24000,sample=50;
	cout << "Input the start and interval of timesteps:";
	cin >> startstep >> intervalstep;
	cout << "Input the numble of sample points:";
	cin >> sample;
	for (int step = 0; step < sample; step++) {
		trjstep = to_string(startstep + step * intervalstep);
		//trjstep = "0";
		Read_lammpstrj mytrj(argv[1], trjstep);
		mytrj.readtrj();
		mytrj.carbon_bond();
		Solve_Circle mycircle(mytrj.ADJmatrix, mytrj.carbon_number);
		mycircle.DFSTraverse();
		cout << "\ncircle nums r: " << mycircle.circle_count(mycircle.c_c,trjstep,outfile) << '\n';
		
	}

	return EXIT_SUCCESS;
}