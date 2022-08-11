#ifndef _LKH_INTERFACE_H
#define _LKH_INTERFACE_H

extern "C"
{
#include "LKH.h"
#include "Genetic.h"
}

#include <Eigen/Eigen>
#include <bits/stdc++.h>

class LKH_TSP_Solver
{
private:
    std::string dir;
    std::string file_name;
    std::string trace_level;
    std::string runs;
    std::string gains23;
    std::string max_trials;
    std::string seed;
    std::string move_type;
    int dimension;
    bool is_symmetric;
    int scale;
    std::vector<int> indices;

protected:
    void createTSPFile(Eigen::MatrixXd &cost_mat);
    void createPARFile();
    void readTSPResultFile();

public:
    LKH_TSP_Solver(std::string dir, std::string file_name);
    void setCostMatrix(Eigen::MatrixXd &mat);
    int solve();
    std::vector<int> getTour();
};

#endif