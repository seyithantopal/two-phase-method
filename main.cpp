#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

#define PARAMETER_SIZE 20
#define RESTRICTION_ROW 10
#define RESTRICTION_COL 21

enum class Sign {
	LESS,
	EQUAL,
	GREATER
};

int pivotColumn = -1;
int pivotRow = -1;
double pivot;
double biggest = 0;
double lowest = 99999;
bool isMax;

void findTheBiggestCj_Zj(vector<double>& cj_zj) {
	for(int i = 0; i < cj_zj.size() - 1; ++i) {
		if(isMax == true) {
			if(cj_zj[i] > biggest) {
				biggest = cj_zj[i];
				pivotColumn = i;
			}
		} else {
			if(cj_zj[i] < biggest) {
				biggest = cj_zj[i];
				pivotColumn = i;
			}
		}
	}
}

void findTheSmallestRatio(vector<double>& ratio) {
	for(int i = 0; i < ratio.size(); ++i) {
		if(ratio[i] > 0) {
			if(ratio[i] < lowest) {
				lowest = ratio[i];
				pivotRow = i;
			}
		}
	}
}

void calculateZj(double* table, vector<double>& zj, vector<int>& cbj, int tableRowSize, int tableColumnSize) {
	double sum;
	if(!zj.empty()) zj.clear();
	for(int i = 0; i < tableColumnSize; ++i) {
		sum = 0;
		for(int j = 0; j < tableRowSize; ++j) {
			sum += *((table + j * tableColumnSize)+i) * cbj[j];
		}
		zj.push_back(sum);
	}
}

void calculateCj_Zj(vector<double>& cj_zj, int* cj, vector<double>& zj, int tableColumnSize) {
	if(!cj_zj.empty()) cj_zj.clear();	
	for(int i = 0; i < tableColumnSize; ++i) {
		if(i == tableColumnSize - 1) {
			cj_zj.push_back(0 - zj[i]);
		} else {
			cj_zj.push_back(*(cj + i) - zj[i]);
		}
	}
}

void reset() {
	pivotColumn = -1;
	pivotRow = -1;
	pivot = 0;
	biggest = 0;
	lowest = 99999;
}


int main() {
	// Variables
	int aCount = 0;
	int sCount = 0;
	int both = 0;
	double value = 0;
	
	vector<double> A;
	vector<double> S;
	vector<double> Solution;
	A.push_back(0);
	S.push_back(0);

	// Settings
	isMax = true;
	array<int, PARAMETER_SIZE> Ze{
		5, 3, 5, 2, 1, 20, 20, 12, 8, 1, 1, 10, 1, 1, 1, 1, 1, 20, 1, 2
	};

	array<Sign, RESTRICTION_ROW> signs{
		Sign::LESS,
		Sign::LESS,
		Sign::LESS,
		Sign::GREATER,
		Sign::GREATER,
		Sign::GREATER,
		Sign::EQUAL,
		Sign::EQUAL,
		Sign::EQUAL,
		Sign::EQUAL
	};

	array<array<int, RESTRICTION_COL>, RESTRICTION_ROW> Re{
		0, 0, 0, 0, 0, 0, 0, 12, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 113,
		2, 0, 0, 0, 0, 0, 0, 0, 4, 5, 0, 0, 0, 0, 0, 0, 6, 2, 1, 1, 125,
		0, 1, 0, 2, 1, 1, 1, 1, 0, 1, 2, 1, 1, 1, 1, 2, 0, 2, 1, 5, 280,
		1, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0, 80,
		0, 1, 2, 0, 5, 1, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 125,
		0, 0, 5, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 2, 1, 2, 1, 0, 150,
		0, 0, 0, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 2, 1, 2, 0, 0, 0, 0, 100,
		0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 120,
		0, 3, 3, 2, 0, 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 0, 0, 0, 0, 4, 230,
		1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 2, 1, 1, 1, 0, 110
	};

	// Adding additinal and slack variables
	for(int i = 0; i < RESTRICTION_ROW; ++i) {
		auto aPos = A.cbegin() + i + 1;
		auto sPos = S.cbegin() + i + 1;

		if(signs.at(i) == Sign::GREATER) {
			A.insert(aPos, 1);
			S.insert(sPos, -1);
			aCount += 1;
			sCount += 1;
			both += 1;
		} else if(signs.at(i) == Sign::EQUAL) {
			A.insert(aPos, 1);
			S.insert(sPos, 0);
			aCount += 1;
		} else if(signs.at(i) == Sign::LESS) {
			A.insert(aPos, 0);
			S.insert(sPos, 1);
			sCount += 1;
		}

	}

	int CjSize = Ze.size() + aCount + sCount;
	int Cj[CjSize] = {0};

	int Cj2[CjSize] = {};
	for(int i = Ze.size() + sCount; i < CjSize; ++i) {
		if(isMax == true) Cj[i] = -1;
		else if(isMax == false) Cj[i] = 1;
	}
	for(int i = 0; i < Ze.size(); ++i) {
		Cj2[i] = Ze[i];
	}	

	// CBj
	vector<int> CBj;
	vector<int> CBj2;

	for(int i = 1; i < A.size(); ++i) {
		

		if((A.at(i) && S.at(i)) != 0) {
			if(isMax == true) CBj.push_back(-1);
			else if(isMax == false) CBj.push_back(1);
			CBj2.push_back(0);
		} else if(S.at(i) == 0 && A.at(i) != 0) {
			if(isMax == true) CBj.push_back(-1);
			else if(isMax == false) CBj.push_back(1);
			CBj2.push_back(0);
		} else if(A.at(i) == 0 && S.at(i) != 0) {
			CBj.push_back(0);
			CBj2.push_back(0);
		}
	}

	int tableRowSize = (aCount + sCount) - both;
	int tableColumnSize = CjSize + 1;
	double table[tableRowSize][tableColumnSize] = {};
	int tablePos = Ze.size(), a;
	
	for(int i = 0; i < tableRowSize; ++i) {
		for(int j = 0; j < Ze.size(); ++j) {
			table[i][j] = Re[i][j];
		}
	}

	for(int i = 0; i < tableRowSize; ++i) {
		a = 0;
		for(int j = 0; j < S.size(); ++j) {
			if(S[j] != 0) {
				if((i + 1) == j) {
					table[i][tablePos + a] = S[i + 1];
				} else {
					table[i][tablePos + a] = 0;
				}
				a += 1;
			}
		}
	}

	tablePos += a;
	for(int i = 0; i < tableRowSize; ++i) {
		a = 0;
		for(int j = 0; j < A.size(); ++j) {
			if(A[j] != 0) {
				if((i + 1) == j) {
					table[i][tablePos + a] = A[i + 1];
				} else {
					table[i][tablePos + a] = 0;
				}
				a += 1;
			}
		}
	}

	tablePos += a;
	for(int i = 0; i < tableRowSize; ++i) {
		a = 0;
		for(int j = 0; j < 1; ++j) {
			table[i][tablePos + a] = Re[i][RESTRICTION_COL - 1];
			a += 1;
		}
	}

	// Zj
	vector<double> Zj;
	calculateZj((double*)table, Zj, CBj, tableRowSize, tableColumnSize);

	// Cj - Zj
	vector<double> Cj_Zj;
	calculateCj_Zj(Cj_Zj, Cj, Zj, tableColumnSize);


	// Ratio
	vector<double> Ratio;

	// Phase-1
	double tableTemp[tableRowSize][tableColumnSize] = {};
	while(true) {
		// Find the biggest Cj-Zj
		findTheBiggestCj_Zj(Cj_Zj);
		if(biggest == 0 && pivotColumn == -1) break;
		
		// Find the Ratio column
		for(int i = 0; i < tableRowSize; ++i) {
			double division = table[i][tablePos] / table[i][pivotColumn];
			Ratio.push_back(division);
		}

		// Find the smallest Ratio
		findTheSmallestRatio(Ratio);
		if(lowest == 99999 && pivotRow == -1) break;
		
		// Find pivot element
		pivot = table[pivotRow][pivotColumn];

		// Divide pivot row by pivot
		for(int j = 0; j < tableColumnSize; ++j) {
			double division = table[pivotRow][j] / pivot;
			tableTemp[pivotRow][j] = division;
		}

		// Change artificial variables with variables
		CBj[pivotRow] = Cj[pivotColumn];
		CBj2[pivotRow] = Cj2[pivotColumn];

		// Fill in the table
		for(int i = 0; i < tableRowSize; ++i) {
			for(int j = 0; j < tableColumnSize; ++j) {
				if(i != pivotRow) {
					tableTemp[i][j] = table[i][j] - (table[pivotRow][j] * table[i][pivotColumn] / pivot);
				}
			} 
		}

		reset();
		Ratio.clear();

		calculateZj((double*)tableTemp, Zj, CBj, tableRowSize, tableColumnSize);

		// Cj - Zj
		calculateCj_Zj(Cj_Zj, Cj, Zj, tableColumnSize);
		findTheBiggestCj_Zj(Cj_Zj);
		if(biggest == 0 && pivotColumn == -1) break;
		copy(&tableTemp[0][0], &tableTemp[0][0] + tableRowSize * tableColumnSize, &table[0][0]);
	}
	// Phase-2
	while(true) {
		copy(&tableTemp[0][0], &tableTemp[0][0] + tableRowSize * tableColumnSize, &table[0][0]);
		
		calculateZj((double*)table, Zj, CBj2, tableRowSize, tableColumnSize);
		calculateCj_Zj(Cj_Zj, Cj2, Zj, tableColumnSize);
		for(int i = 0; i < tableRowSize; ++i) {
			for(int j = Ze.size() + sCount; j < tableColumnSize - 1; ++j) {
				table[i][j] = 0;
				Cj_Zj[j] = 0;
				Zj[j] = 0;
			}
		}

		findTheBiggestCj_Zj(Cj_Zj);
		if(biggest == 0 && pivotColumn == -1) break;

		// Find the Ratio column
		for(int i = 0; i < tableRowSize; ++i) {
			double division = table[i][tablePos] / tableTemp[i][pivotColumn];
			Ratio.push_back(division);
		}

		// Find the smallest Ratio
		findTheSmallestRatio(Ratio);
		if(lowest == 99999 && pivotRow == -1) break;
		
		// Find pivot element
		pivot = table[pivotRow][pivotColumn];

		// Divide pivot row by pivot
		for(int j = 0; j < tableColumnSize; ++j) {
			double division = table[pivotRow][j] / pivot;
			tableTemp[pivotRow][j] = division;
		}

		// Change artificial variables with variables
		CBj2[pivotRow] = Cj2[pivotColumn];

		// Fill in the table
		for(int i = 0; i < tableRowSize; ++i) {
			for(int j = 0; j < tableColumnSize; ++j) {
				if(i != pivotRow) {
					tableTemp[i][j] = table[i][j] - (table[pivotRow][j] * table[i][pivotColumn] / pivot);
				}
			} 
		}

		reset();
		Ratio.clear();

		calculateZj((double*)tableTemp, Zj, CBj2, tableRowSize, tableColumnSize);

		// Cj - Zj
		calculateCj_Zj(Cj_Zj, Cj2, Zj, tableColumnSize);
		findTheBiggestCj_Zj(Cj_Zj);
		if(biggest == 0 && pivotColumn == -1) break;
		copy(&tableTemp[0][0], &tableTemp[0][0] + tableRowSize * tableColumnSize, &table[0][0]);
	}

	for(int i = 0; i < Ze.size(); ++i) {
		if(Cj_Zj[i] == 0) {
			Solution.push_back(tableTemp[i][tablePos]);
		} else {
			Solution.push_back(0);
		}
	}
	value = Zj[tableColumnSize - 1];
	

	cout << "tableTemp: " << endl;;
	for(int i = 0; i < tableRowSize; ++i) {
		for(int j = 0; j < tableColumnSize; ++j) {
			cout << fixed << setprecision(2) << tableTemp[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;

	cout << "Solution: \n";
	for(int i = 0; i < Solution.size(); ++i) {
		cout << "X" << i+1 << " = " << Solution[i] << endl;
	}
	cout << (isMax == true ? "Max Value: " : "Min Value: ") << fixed << value << endl;
}