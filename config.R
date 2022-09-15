MASS_RANGE = c(50,2000)

# fraction of total retention time, or of ranks of retention time
# used to determine coelution of ions ad hoc
RETENTION_TIME_TOLERANCE_FRAC = 0.02

pos_default <- c('M[1+]', 'M+H[1+]', 'M+2H[2+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]','M+Na[1+]', 'M+H+Na[2+]', 'M+HCOONa[1+]')

primary_ions = c('M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]')

dict_weight_adduct = c( 'M[1+]'=5, 
                        'M+H[1+]'=5,
                        'M+2H[2+]'=3,
                        'M+3H[3+]'=1,
                        'M(C13)+H[1+]'=2,
                        'M(C13)+2H[2+]'=1,
                        'M(C13)+3H[3+]'=1,
                        'M(S34)+H[1+]'=1,
                        'M(Cl37)+H[1+]'=1,
                        'M+Na[1+]'=3, 
                        'M+H+Na[2+]'=2,
                        'M+K[1+]'=2, 
                        'M+H2O+H[1+]'=1, 
                        'M-H2O+H[1+]'=1, 
                        'M-H4O2+H[1+]'=1,
                        'M-NH3+H[1+]'=1,
                        'M-CO+H[1+]'=1,
                        'M-CO2+H[1+]'=1,
                        'M-HCOOH+H[1+]'=1,
                        'M+HCOONa[1+]'=1,
                        'M-HCOONa+H[1+]'=1,
                        'M+NaCl[1+]'=1, 
                        'M-C3H4O2+H[1+]'=1,
                        'M+HCOOK[1+]'=1,
                        'M-HCOOK+H[1+]'=1,
                        # negative
                        'M-H[-]'=5,
                        'M-2H[2-]'=3,
                        'M(C13)-H[-]'=2,
                        'M(S34)-H[-]'=1,
                        'M(Cl37)-H[-]'=1,
                        'M+Na-2H[-]'=2,
                        'M+K-2H[-]'=1,
                        'M-H2O-H[-]'=1,
                        'M+Cl[-]'=1,
                        'M+Cl37[-]'=1,
                        'M+Br[-]'=1,
                        'M+Br81[-]'=1,
                        'M+ACN-H[-]'=1,
                        'M+HCOO[-]'=1,
                        'M+CH3COO[-]'=1,
                        'M-H+O[-]'=1)

SEARCH_STEPS = 4
MODULE_SIZE_LIMIT = 100

NUM_PERM = 100

SIGNIFICANCE_CUTOFF = 0.05



