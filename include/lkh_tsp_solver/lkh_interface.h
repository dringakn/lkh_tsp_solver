#ifndef _LKH_INTERFACE_H
#define _LKH_INTERFACE_H

extern "C"
{
#include "LKH.h"
#include "Genetic.h"
}

/**
 * @brief The interface function for the TSP problem. "~/personal_ws/src/lkh_tsp_solver/resource/single.par"
 * 
 * @param input_file The name of the input file (file_name.par) to read the problem specifications.
 * @return int 0 if success.
 */
int solveTSPLKH(const char *input_file);

#endif