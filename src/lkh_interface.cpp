#include <lkh_tsp_solver/lkh_interface.h>

LKH_TSP_Solver::LKH_TSP_Solver(std::string dir, std::string file_name)
{
  this->dir = dir;
  this->file_name = file_name;
  this->trace_level = "1";
  this->runs = "10";
  this->move_type = "2";
  this->seed = "1";
  this->gains23 = "NO";
  this->is_symmetric = false;
  this->scale = 1000;
}

void LKH_TSP_Solver::createTSPFile(Eigen::MatrixXd &cost_mat)
{
  // Write params and cost matrix to problem file
  // std::cout << (this->dir + this->file_name + ".tsp").c_str() << std::endl;
  std::string fn("");
  fn.append(this->dir).append(this->file_name).append(".tsp");
  std::ofstream prob_file(fn.c_str());

  // Problem specification part, follow the format of TSPLIB (Matrix Format):
  std::string prob_spec = "NAME : " + this->file_name + "\n" +
                          "TYPE : " + ((this->is_symmetric) ? "TSP" : "ATSP") + "\n" +
                          "COMMENT : Dr. -Ing. Ahmad Kamal Nasir\n" +
                          "DIMENSION : " + std::to_string(this->dimension) + "\n" +
                          "EDGE_WEIGHT_TYPE : EXPLICIT\n" +
                          "EDGE_WEIGHT_FORMAT : FULL_MATRIX\n" +
                          "EDGE_WEIGHT_SECTION\n";

  prob_file << prob_spec;
  // Problem data part
  if (this->is_symmetric)
  {
    // Use symmetric TSP
    for (int i = 1; i < this->dimension; ++i)
    {
      for (int j = 0; j < i; ++j)
      {
        int int_cost = cost_mat(i, j) * this->scale;
        prob_file << int_cost << " ";
      }
      prob_file << "\n";
    }
  }
  else
  {
    // Use Asymmetric TSP
    for (int i = 0; i < this->dimension; ++i)
    {
      for (int j = 0; j < this->dimension; ++j)
      {
        int int_cost = cost_mat(i, j) * this->scale;
        prob_file << int_cost << " ";
      }
      prob_file << "\n";
    }
  }

  prob_file << "EOF";
  prob_file.close();
}

void LKH_TSP_Solver::createPARFile()
{
  std::string fn("");
  fn.append(this->dir).append(this->file_name).append(".par");
  std::ofstream par_file(fn.c_str());
  fn.clear();
  par_file << "PROBLEM_FILE = " << fn.append(this->dir).append(this->file_name).append(".tsp") << "\n";
  fn.clear();
  par_file << "OUTPUT_TOUR_FILE = " << fn.append(this->dir).append(this->file_name).append(".txt") << "\n";
  par_file << "TRACE_LEVEL = " << this->trace_level << "\n";
  par_file << "RUNS = " << this->runs << "\n";
  par_file << "SEED = " << this->seed << "\n";
  par_file << "MOVE_TYPE = " << this->move_type << "\n";
  par_file << "MAX_TRIALS = " << this->dimension << "\n";
  par_file << "GAIN23 = " << this->gains23 << "\n";
  par_file.close();
}

void LKH_TSP_Solver::readTSPResultFile()
{
  // Read optimal tour from the tour section of result file
  std::string fn("");
  fn.append(this->dir).append(this->file_name).append(".txt");
  std::ifstream res_file(fn.c_str());
  std::string res;
  while (getline(res_file, res))
  {
    // Go to tour section
    if (res.compare("TOUR_SECTION") == 0)
      break;
  }

  this->indices.clear();
  if (this->is_symmetric)
  {
    // Read path for Symmetric TSP formulation
    getline(res_file, res); // Skip current pose
    getline(res_file, res);
    int id = stoi(res);
    bool rev = (id == this->dimension); // The next node is virutal depot?

    while (id != -1)
    {
      this->indices.push_back(id - 2);
      getline(res_file, res);
      id = stoi(res);
    }
    if (rev)
      reverse(this->indices.begin(), this->indices.end());
    this->indices.pop_back(); // Remove the depot
  }
  else
  {
    // Read path for ATSP formulation
    while (getline(res_file, res))
    {
      // Read indices of frontiers in optimal tour
      int id = stoi(res);
      // if (id == 1) // Ignore the current state
      //   continue;
      if (id == -1)
        break;
      this->indices.push_back(id-1); // Subtract 1 to create 0 based index
      // this->indices.push_back(id - 2); // Idx of solver-2 == Idx of frontier
    }
  }

  res_file.close();
}

std::vector<int> LKH_TSP_Solver::getTour()
{
  return this->indices;
}

int LKH_TSP_Solver::solve()
{
  GainType Cost, OldOptimum;
  double Time, LastTime = GetTime();

  /* Read the specification of the problem */
  std::string fn("");
  fn.append(this->dir).append(this->file_name).append(".par");
  ParameterFileName = const_cast<char *>(fn.c_str());
  std::cout << "Parameter File: " << ParameterFileName << std::endl;
  ReadParameters();
  MaxMatrixDimension = 20000;
  MergeWithTour = Recombination == IPT ? MergeWithTourIPT : MergeWithTourGPX2;
  ReadProblem();

  if (SubproblemSize > 0)
  {
    if (DelaunayPartitioning)
      SolveDelaunaySubproblems();
    else if (KarpPartitioning)
      SolveKarpSubproblems();
    else if (KCenterPartitioning)
      SolveKCenterSubproblems();
    else if (KMeansPartitioning)
      SolveKMeansSubproblems();
    else if (RohePartitioning)
      SolveRoheSubproblems();
    else if (MoorePartitioning || SierpinskiPartitioning)
      SolveSFCSubproblems();
    else
      SolveTourSegmentSubproblems();
    return EXIT_SUCCESS;
  }
  AllocateStructures();
  CreateCandidateSet();
  InitializeStatistics();

  if (Norm != 0)
    BestCost = PLUS_INFINITY;
  else
  {
    /* The ascent has solved the problem! */
    Optimum = BestCost = (GainType)LowerBound;
    UpdateStatistics(Optimum, GetTime() - LastTime);
    RecordBetterTour();
    RecordBestTour();
    WriteTour(OutputTourFileName, BestTour, BestCost);
    WriteTour(TourFileName, BestTour, BestCost);
    Runs = 0;
  }

  /* Find a specified number (Runs) of local optima */
  for (Run = 1; Run <= Runs; Run++)
  {
    LastTime = GetTime();
    Cost = FindTour(); /* using the Lin-Kernighan heuristic */
    if (MaxPopulationSize > 1)
    {
      /* Genetic algorithm */
      int i;
      for (i = 0; i < PopulationSize; i++)
      {
        GainType OldCost = Cost;
        Cost = MergeTourWithIndividual(i);
        if (TraceLevel >= 1 && Cost < OldCost)
        {
          printff("  Merged with %d: Cost = " GainFormat, i + 1, Cost);
          if (Optimum != MINUS_INFINITY && Optimum != 0)
            printff(", Gap = %0.4f%%", 100.0 * (Cost - Optimum) / Optimum);
          printff("\n");
        }
      }
      if (!HasFitness(Cost))
      {
        if (PopulationSize < MaxPopulationSize)
        {
          AddToPopulation(Cost);
          if (TraceLevel >= 1)
            PrintPopulation();
        }
        else if (Cost < Fitness[PopulationSize - 1])
        {
          i = ReplacementIndividual(Cost);
          ReplaceIndividualWithTour(i, Cost);
          if (TraceLevel >= 1)
            PrintPopulation();
        }
      }
    }
    else if (Run > 1)
      Cost = MergeTourWithBestTour();
    if (Cost < BestCost)
    {
      BestCost = Cost;
      RecordBetterTour();
      RecordBestTour();
      WriteTour(OutputTourFileName, BestTour, BestCost);
      WriteTour(TourFileName, BestTour, BestCost);
    }
    OldOptimum = Optimum;
    if (Cost < Optimum)
    {
      if (FirstNode->InputSuc)
      {
        Node *N = FirstNode;
        while ((N = N->InputSuc = N->Suc) != FirstNode)
          ;
      }
      Optimum = Cost;
      printff("*** New optimum = " GainFormat " ***\n\n", Optimum);
    }
    Time = fabs(GetTime() - LastTime);
    UpdateStatistics(Cost, Time);
    if (TraceLevel >= 1 && Cost != PLUS_INFINITY)
    {
      // printff("Run %d: Cost = " GainFormat, Run, Cost);
      // if (Optimum != MINUS_INFINITY && Optimum != 0)
      //   printff(", Gap = %0.4f%%", 100.0 * (Cost - Optimum) / Optimum);
      // printff(", Time = %0.2f sec. %s\n\n", Time,
      //         Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
    }
    if (StopAtOptimum && Cost == OldOptimum && MaxPopulationSize >= 1)
    {
      Runs = Run;
      break;
    }
    if (PopulationSize >= 2 && (PopulationSize == MaxPopulationSize || Run >= 2 * MaxPopulationSize) &&
        Run < Runs)
    {
      Node *N;
      int Parent1, Parent2;
      Parent1 = LinearSelection(PopulationSize, 1.25);
      do
        Parent2 = LinearSelection(PopulationSize, 1.25);
      while (Parent2 == Parent1);
      ApplyCrossover(Parent1, Parent2);
      N = FirstNode;
      do
      {
        if (ProblemType != HCP && ProblemType != HPP)
        {
          int d = C(N, N->Suc);
          AddCandidate(N, N->Suc, d, INT_MAX);
          AddCandidate(N->Suc, N, d, INT_MAX);
        }
        N = N->InitialSuc = N->Suc;
      } while (N != FirstNode);
    }
    SRandom(++Seed);
  }
  PrintStatistics();
  this->readTSPResultFile();
  return EXIT_SUCCESS;
}

void LKH_TSP_Solver::setCostMatrix(Eigen::MatrixXd &cost_mat)
{
  this->dimension = cost_mat.rows();
  this->createTSPFile(cost_mat); // Ucomment after testing
  // this->createPARFile();         // Uncomment after testing
}
