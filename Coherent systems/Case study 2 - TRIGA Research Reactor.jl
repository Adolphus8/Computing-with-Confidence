## Import the necessary packages:
using ProbabilityBoundsAnalysis, UncLogic, PyPlot, Tables, CSV

################################################################################################################################################
## Consider the Nuclear reactor TRR-1/M1 with 6 components:

# Valve-3: No. of working states, k1 = 8; No. of tests, n1 = 42 yr;
# Valve-3: Operator failure, p = 3.49E-05 yr-1
# Valve-4: No. of working states, k2 = 8; No. of tests, n2 = 42 yr;
# Valve-3: Operator failure, p = 3.49E-05 yr-1
# Inlet pipe: No. of working states, k3 = 1; No. of tests, n3 = 42 yr;

C1 = KN(1, 42);                               # Component 1 - Inlet pipe (Mechanical failure) [yr-1]
C2 = KN(8, 42);                               # Component 2 - Valve 3 (Mechanical failure) [yr-1]
C3 = makepbox(interval(3.49e-05, 3.49e-05));  # Component 3 - Valve 3 (Human operator failure) [yr-1]
C4 = KN(8, 42);                               # Component 4 - Valve 4 (Mechanical failure) [yr-1]
C5 = makepbox(interval(3.49e-05, 3.49e-05));  # Component 5 - Valve 4 (Human operator failure) [yr-1]

# Case 1: Independence between in the AND / OR relations:
rho = 0; a = 0.025;

E4 = or(C4, C5, rho); cutE4a_indep = cut(E4, a); cutE4b_indep  = cut(E4, 1-a); 
E3 = or(C2, C3, rho); cutE3a_indep = cut(E3, a); cutE3b_indep  = cut(E3, 1-a);
E2 = or(E3, E4, rho); cutE2a_indep = cut(E2, a); cutE2b_indep  = cut(E2, 1-a); 
E1 = and(C1, E2, rho); cutE1a_indep = cut(E1, a); cutE1b_indep  = cut(E1, 1-a); 

R_cbox_indep = zeros(200,2); cbox = E1;
R_cbox_indep[:,1] = cbox.d; R_cbox_indep[:,2] = cbox.u; table_R = Tables.table(R_cbox_indep);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E1_cbox_indep_TRIGA.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_indep = zeros(200,2); cbox = E2;
R_cbox_indep[:,1] = cbox.d; R_cbox_indep[:,2] = cbox.u; table_R = Tables.table(R_cbox_indep);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E2_cbox_indep_TRIGA.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_indep = zeros(200,2); cbox = E3;
R_cbox_indep[:,1] = cbox.d; R_cbox_indep[:,2] = cbox.u; table_R = Tables.table(R_cbox_indep);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E3_cbox_indep_TRIGA.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_indep = zeros(200,2); cbox = E4;
R_cbox_indep[:,1] = cbox.d; R_cbox_indep[:,2] = cbox.u; table_R = Tables.table(R_cbox_indep);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E4_cbox_indep_TRIGA.csv", table_R) # Save Pbox output cbox as CSV


# Case 2: Uncertain dependence between in the AND / OR relations:
rho = interval(-1, 1); a = 0.025;

E4 = or(C4, C5, rho); cutE4a_frechet = cut(E4, a); cutE4b_frechet  = cut(E4, 1-a); 
E3 = or(C2, C3, rho); cutE3a_frechet = cut(E3, a); cutE3b_frechet = cut(E3, 1-a);  
E2 = or(E3, E4, rho); cutE2a_frechet = cut(E2, a); cutE2b_frechet  = cut(E2, 1-a); 
E1 = and(C1, E2, rho); cutE1a_frechet = cut(E1, a); cutE1b_frechet  = cut(E1, 1-a); 

R_cbox_frechet = zeros(200,2); cbox = E1;
R_cbox_frechet[:,1] = cbox.d; R_cbox_frechet[:,2] = cbox.u; table_R = Tables.table(R_cbox_frechet);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E1_cbox_frechet_TRIGA.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_frechet = zeros(200,2); cbox = E2;
R_cbox_frechet[:,1] = cbox.d; R_cbox_frechet[:,2] = cbox.u; table_R = Tables.table(R_cbox_frechet);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E2_cbox_frechet_TRIGA.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_frechet = zeros(200,2); cbox = E3;
R_cbox_frechet[:,1] = cbox.d; R_cbox_frechet[:,2] = cbox.u; table_R = Tables.table(R_cbox_frechet);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E3_cbox_frechet_TRIGA.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_frechet = zeros(200,2); cbox = E4;
R_cbox_frechet[:,1] = cbox.d; R_cbox_frechet[:,2] = cbox.u; table_R = Tables.table(R_cbox_frechet);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E4_cbox_frechet_TRIGA.csv", table_R) # Save Pbox output cbox as CSV

################################################################################################################################################
## Consider the Pressurised tank system with 6 components using the mean reliability estimates:

C1 = 1/42;      # Component 1 - Inlet pipe (Mechanical failure) [yr-1]
C2 = 8/42;      # Component 2 - Valve 3 (Mechanical failure) [yr-1]
C3 = 3.49E-05;  # Component 3 - Valve 3 (Human operator failure) [yr-1]
C4 = 8/42;      # Component 4 - Valve 4 (Mechanical failure) [yr-1]
C5 = 3.49E-05;  # Component 5 - Valve 4 (Human operator failure) [yr-1]

# Case 1: Independence between in the AND / OR relations:
rho = 0; 

E4 = or(C4, C5, rho);  
E3 = or(C2, C3, rho);  
E2 = or(E3, E4, rho);  
E1 = and(C1, E2, rho); 

# Case 2: Uncertain dependence between in the AND / OR relations:
rho = interval(-1, 1); 

E4 = or(C4, C5, rho);  
E3 = or(C2, C3, rho);  
E2 = or(E3, E4, rho);  
E1 = and(C1, E2, rho); 