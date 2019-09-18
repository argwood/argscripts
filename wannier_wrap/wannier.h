#ifndef WANNIER_H
#define WANNIER_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;


int count_wannier(vector<Snapshot> &traj, struct Params params, int snap, int i, vector<int> &assigned);

void wannier(vector<Snapshot> &traj, struct Params params);

#endif
