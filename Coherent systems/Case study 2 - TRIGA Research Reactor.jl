## Import the necessary packages:
using ProbabilityBoundsAnalysis, UncLogic, PyPlot, Tables, CSV

################################################################################################################################################
## Consider the Nuclear reactor TRR-1/M1 with 6 components:

# Valve-3: No. of working states, k1 = 0; No. of tests, n1 = 2000;
# Valve-3: Operator failure, p = 3.49E-05
# Valve-4: No. of working states, k2 = 0; No. of tests, n2 = 2000;
# Valve-3: Operator failure, p = 3.49E-05
# Inlet pipe: No. of working states, k5 = 7; No. of tests, n5 = 10000;
# Component 6: No. of working states, k6 = 0; No. of tests, n6 = 380;

C1 = KN(0, 2000);  # Component 1
C2 = KN(3, 500);   # Component 2
C3 = KN(1, 150);   # Component 3
C4 = KN(0, 460);   # Component 4
C5 = KN(7, 10000); # Component 5
C6 = KN(0, 380);   # Component 6

# Case 1: Independence between in the AND / OR relations:
rho = 0; a = 0.025;

E5 = or(C5, C6, rho); cutE5a_indep = cut(E5, a); cutE5b_indep  = cut(E5, 1-a);
E4 = or(C4, E5, rho); cutE4a_indep = cut(E4, a); cutE4b_indep  = cut(E4, 1-a);
E3 = and(C3, E4, rho); cutE3a_indep = cut(E3, a); cutE3b_indep  = cut(E3, 1-a);
E2 = or(C2, E3, rho); cutE2a_indep = cut(E2, a); cutE2b_indep  = cut(E2, 1-a);
E1 = or(C1, E2, rho); cutE1a_indep = cut(E1, a); cutE1b_indep  = cut(E1, 1-a);

R_cbox_indep = zeros(200,2); cbox = E1;
R_cbox_indep[:,1] = cbox.d; R_cbox_indep[:,2] = cbox.u; table_R = Tables.table(R_cbox_indep);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E1_cbox_indep.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_indep = zeros(200,2); cbox = E2;
R_cbox_indep[:,1] = cbox.d; R_cbox_indep[:,2] = cbox.u; table_R = Tables.table(R_cbox_indep);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E2_cbox_indep.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_indep = zeros(200,2); cbox = E3;
R_cbox_indep[:,1] = cbox.d; R_cbox_indep[:,2] = cbox.u; table_R = Tables.table(R_cbox_indep);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E3_cbox_indep.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_indep = zeros(200,2); cbox = E4;
R_cbox_indep[:,1] = cbox.d; R_cbox_indep[:,2] = cbox.u; table_R = Tables.table(R_cbox_indep);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E4_cbox_indep.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_indep = zeros(200,2); cbox = E5;
R_cbox_indep[:,1] = cbox.d; R_cbox_indep[:,2] = cbox.u; table_R = Tables.table(R_cbox_indep);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E5_cbox_indep.csv", table_R) # Save Pbox output cbox as CSV


# Case 2: Uncertain dependence between in the AND / OR relations:
rho = interval(-1, 1); a = 0.025;

E5 = or(C5, C6, rho); cutE5a_frechet = cut(E5, a); cutE5b_frechet  = cut(E5, 1-a);
E4 = or(C4, E5, rho); cutE4a_frechet = cut(E4, a); cutE4b_frechet  = cut(E4, 1-a);
E3 = and(C3, E4, rho); cutE3a_frechet = cut(E3, a); cutE3b_frechet  = cut(E3, 1-a);
E2 = or(C2, E3, rho); cutE2a_frechet = cut(E2, a); cutE2b_frechet  = cut(E2, 1-a);
E1 = or(C1, E2, rho); cutE1a_frechet = cut(E1, a); cutE1b_frechet  = cut(E1, 1-a);

R_cbox_frechet = zeros(200,2); cbox = E1;
R_cbox_frechet[:,1] = cbox.d; R_cbox_frechet[:,2] = cbox.u; table_R = Tables.table(R_cbox_frechet);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E1_cbox_frechet.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_frechet = zeros(200,2); cbox = E2;
R_cbox_frechet[:,1] = cbox.d; R_cbox_frechet[:,2] = cbox.u; table_R = Tables.table(R_cbox_frechet);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E2_cbox_frechet.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_frechet = zeros(200,2); cbox = E3;
R_cbox_frechet[:,1] = cbox.d; R_cbox_frechet[:,2] = cbox.u; table_R = Tables.table(R_cbox_frechet);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E3_cbox_frechet.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_frechet = zeros(200,2); cbox = E4;
R_cbox_frechet[:,1] = cbox.d; R_cbox_frechet[:,2] = cbox.u; table_R = Tables.table(R_cbox_frechet);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E4_cbox_frechet.csv", table_R) # Save Pbox output cbox as CSV

R_cbox_frechet = zeros(200,2); cbox = E5;
R_cbox_frechet[:,1] = cbox.d; R_cbox_frechet[:,2] = cbox.u; table_R = Tables.table(R_cbox_frechet);
CSV.write("C:\\Users\\snrltsa\\Desktop\\E5_cbox_frechet.csv", table_R) # Save Pbox output cbox as CSV

################################################################################################################################################
## Consider the Pressurised tank system with 6 components using the mean reliability estimates:

C1 = 0;  
C2 = 3/500;
C3 = 1/150;
C4 = 0;
C5 = 7/10000;
C6 = 0;

# Case 1: Independence between in the AND / OR relations:
rho = 0; 

E5 = or(C5, C6, rho); 
E4 = or(C4, E5, rho); 
E3 = and(C3, E4, rho); 
E2 = or(C2, E3, rho); 
E1 = or(C1, E2, rho); 

# Case 2: Uncertain dependence between in the AND / OR relations:
rho = interval(-1, 1); 

E5 = or(C5, C6, rho); 
E4 = or(C4, E5, rho); 
E3 = and(C3, E4, rho); 
E2 = or(C2, E3, rho); 
E1 = or(C1, E2, rho); 