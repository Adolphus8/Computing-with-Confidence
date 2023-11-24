## Import the necessary packages:
using ProbabilityBoundsAnalysis, UncLogic, PyPlot, Tables, CSV

################################################################################################################################################
## Consider the Parallel-Series bridge system with 5 components:

# Component 1: No. of working states, k1 = 23; No. of tests, n1 = 24;
# Component 2: No. of working states, k2 = 23; No. of tests, n2 = 24;
# Component 3: No. of working states, k3 = 14; No. of tests, n3 = 17;
# Component 4: No. of working states, k4 = 14; No. of tests, n4 = 17;
# Component 5: No. of working states, k5 = 12; No. of tests, n5 = 12;

C1 = KN(23, 24); # Component 1
C2 = KN(23, 24); # Component 2
C3 = KN(14, 17); # Component 3
C4 = KN(14, 17); # Component 4
C5 = KN(12, 12); # Component 5

# Case 1: Independence between components, indepedence between path-sets:
rho_com = 0; rho_path = 0; a = 0.025;
path1 = and(C1, C3, rho_com);
path2 = and(C2, C4, rho_com);
path3 = and(C1, and(C4, C5, rho_com), rho_com);
path4 = and(C2, and(C3, C5, rho_com), rho_com);

R_case1 = or(path4, or(path3, or(path1, path2, rho_path), rho_path), rho_path);
cut1a = cut(R_case1, a); cut1b = cut(R_case1, 1-a);

# Case 2: Independence between components, uncertain dependence between path-sets:
rho_com = 0; rho_path = interval(-1,1); a = 0.025;
path1 = and(C1, C3, rho_com);
path2 = and(C2, C4, rho_com);
path3 = and(C1, and(C4, C5, rho_com), rho_com);
path4 = and(C2, and(C3, C5, rho_com), rho_com);

R_case2 = or(path4, or(path3, or(path1, path2, rho_path), rho_path), rho_path);
cut2a = cut(R_case2, a); cut2b = cut(R_case1, 1-a);

# Case 3: Uncertain dependence between components, independence between path-sets:
rho_com = interval(-1,1); rho_path = 0; a = 0.025;
path1 = and(C1, C3, rho_com);
path2 = and(C2, C4, rho_com);
path3 = and(C1, and(C4, C5, rho_com), rho_com);
path4 = and(C2, and(C3, C5, rho_com), rho_com);

R_case3 = or(path4, or(path3, or(path1, path2, rho_path), rho_path), rho_path);
cut3a = cut(R_case3, a); cut3b = cut(R_case1, 1-a);

# Case 4: Uncertain dependence between components, uncertain dependence between path-sets:
rho_com = interval(-1,1); rho_path = interval(-1,1); a = 0.025;
path1 = and(C1, C3, rho_com);
path2 = and(C2, C4, rho_com);
path3 = and(C1, and(C4, C5, rho_com), rho_com);
path4 = and(C2, and(C3, C5, rho_com), rho_com);

R_case4 = or(path4, or(path3, or(path1, path2, rho_path), rho_path), rho_path);
cut4a = cut(R_case4, a); cut4b = cut(R_case1, 1-a);