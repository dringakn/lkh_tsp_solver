# Sample LKH TSP solver parameter file format

Specifies the level of detail of the output given during the solution process.
The value 0 signifies a minimum amount of output. The higher the value is the
more information is given. Default: 1.
`TRACE_LEVEL = 100`

Specifies the name of the problem file
`PROBLEM_FILE = /home/ahmad/personal_ws/src/lkh_tsp_solver/resource/pr2392.tsp`

The input tour file
`INPUT_TOUR_FILE = /home/ahmad/personal_ws/src/lkh_tsp_solver/resource/pr2392.txt`

Specifies the name of a file to which the best tour is to be written.
`TOUR_FILE = /home/ahmad/personal_ws/src/lkh_tsp_solver/resource/pr2392.txt`

Specifies the name of a file to which penalties (π-values determined by the ascent)
is to be written. If the file already exits, the penalties are read from the file, and the
ascent is skipped.
`PI_FILE =` 


The total number of runs.
`RUNS = 10`

The maximum number of trials in each run. Default: number of nodes (DIMENSION, given in the problem file)
`MAX_TRIALS = 2392`

Known optimal tour length. A run will be terminated as soon as a tour length less
than or equal to optimum is achieved. Default: DBL_MAX.
`OPTIMUM = 378032`

The maximum number of candidate edges to be associated with each node.
The integer may be followed by the keyword SYMMETRIC, signifying that the
candidate set is to be complemented such that every candidate edge is associated
with both its two end nodes. Default: 5.
`MAX_CANDIDATES = 5 SYMMETRIC`

The number of candidate edges to be associated with each node during the ascent.
The candidate set is complemented such that every candidate edge is associated
with both its two end nodes.Default: 50.
`ASCENT_CANDIDATES = 50`

The maximum α-value allowed for any candidate edge is set to EXCESS times the
absolute value of the lower bound of a solution tour (determined by the ascent). Default: 1.0/DIMENSION.
`EXCESS = 0`

The length of the first period in the ascent. Default: DIMENSION/2 (but at least 100).
`INITIAL_PERIOD = 100`

The initial step size used in the ascent. Default: 1.
`INITIAL_STEP_SIZE = 1`

The internal precision in the representation of transformed distances:
dij = PRECISION*cij + πi + πj, where dij, cij, πi and πj are all integral.Default: 100 (which corresponds to 2 decimal places).
`PRECISION = 100`

Specifies the initial seed for random number generation.Default: 1.
`SEED = 1`

Specifies whether the π-values should be determined by subgradient optimization. Default: YES.
`SUBGRADIENT = YES`

Specifies the move type to be used in local search. The value can be 2, 3, 4 or 5 and signifies
whether 2-opt, 3-opt, 4-opt or 5-opt move is to be used. Default: 5.
`MOVE_TYPE = 5`

`PATCHING_C = 3`

`PATCHING_A = 2`

```
BACKBONE_TRIALS = <integer>
BACKTRACKING = {YES | NO }
CANDIDATE_SET_TYPE = {ALPHA | DELAUNAY [PURE ] | NEAREST-NEIGHBOR | QUADRANT }
EXTRA_CANDIDATES = <integer> [SYMMETRIC ]
EXTRA_CANDIDATE_SET_TYPE = {NEAREST-NEIGHBOR | QUADRANT }
GAIN_CRITERION = {YES | NO }
INITIAL_TOUR_ALGORITHM = {BORUVKA | GREEDY | NEAREST-NEIGHBOR | QUICK-BORUVKA | SIERPINSKI | WALK }
INITIAL_TOUR_FRACTION = <real in [0;1]>
KICKS = <integer>
KICK_TYPE = <integer>
MAX_BREADTH = <integer>
MERGE_TOUR_FILE = <string>
NONSEQUENTIAL_MOVE_TYPE = <integer>
PATCHING_A = <integer> [RESTRICTED | EXTENDED ]
PATCHING_C = <integer> [RESTRICTED | EXTENDED ]
SUBPROBLEM_SIZE = <integer> [DELAUNAY | KARP | K-MEANS | ROHE | SIERPINSKI ] [COMPRESSED ] [BORDERS ]
SUBPROBLEM_TOUR_FILE = <string>}
SUBSEQUENT_MOVE_TYPE = <integer>
SUBSEQUENT_PATCHING = {YES | NO }
```