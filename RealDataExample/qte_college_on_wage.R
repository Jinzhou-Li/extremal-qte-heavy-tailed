# The code in this file is used to estimate (extremal) QTEs of college education on wage data.
library(haven)
library(boot)
library(qte)

source("Functions/LoadFunctions.R")

##### Read data
df.raw <- read_dta("RealDataExample/hhv-neo_v3.dta")
df.raw <- df.raw[!is.na(df.raw$wage) & (df.raw$collegedeg | df.raw$somecollege | df.raw$hsgrad),]       # select only people which graduaded high-school
df.nan <- df.raw[,c("wage", "black", "hisp", "south", "west", "northeast", "curban", "broken", 
                                    "age80", "mhgc_mi", "fhgc_mi", "faminc79_th", "numsibs", 
                                    "sasvab1", "sasvab2", "sasvab3", "sasvab4", "sasvab5", "sasvab6" ,
                                    "sgr9_lang_gpa", "sgr9_math_gpa", "sgr9_scosci_gpa", "sgr9_sci_gpa"
                                    )]  
df.nan <- cbind(college = (df.raw$collegedeg | df.raw$somecollege), df.nan)                             # college is 1 if participant went to college for some time
df <- df.nan[complete.cases(df.nan),]                                                                   # remove missing values

###########################################################################
#####                         Results from CV                        ######
###########################################################################

## x.formula is the result of the cross-validtion procedure in crossval.R
## As crossval.R takes some time to run, we also provide the result

# Uncomment the following line if you want to execute CV and comment out the definition of base.formula
# source("crossval.R")
# x.formula <- paste("~", strsplit(base.formula, split="~")[[1]][2])

x.formula <- "~black+hisp+south+west+northeast+curban+broken+age80+
  mhgc_mi+fhgc_mi+faminc79_th+numsibs+sasvab1+sasvab2+sasvab3+
  sasvab4+sasvab5+sasvab6+sgr9_lang_gpa+sgr9_math_gpa+sgr9_scosci_gpa+
  sgr9_sci_gpa+I(mhgc_mi^2)+I(fhgc_mi^2)+I(faminc79_th^2)+I(numsibs^2)+
  I(sasvab1^2)+I(sasvab5^2)+black:west+black:age80+black:mhgc_mi+black:sasvab6+
  hisp:curban+south:broken+south:sasvab1+south:sgr9_scosci_gpa+west:age80+
  west:mhgc_mi+west:fhgc_mi+west:sasvab5+west:sgr9_scosci_gpa+curban:sasvab2+
  age80:faminc79_th+age80:sasvab1+age80:sasvab5+sasvab1:sgr9_lang_gpa+sasvab2:sasvab6+
  sasvab2:sgr9_lang_gpa+sasvab4:sgr9_sci_gpa"

Y = df$wage
D = df$college

###########################################################################
#####                         Extremal QTE                           ######
###########################################################################
set.seed(123)

qte_level <- 0.999
ks <- 35

fit <- glm(as.formula(paste("college", x.formula)), data=df, family="binomial")
prop.pred <- as.numeric(fitted(fit))

# zhang's b out of n bootstrap
set.seed(2021)
qte_firpo_zhang_result <- qte_firpo_zhang(Y, X=NULL, D, pn=1-qte_level, CI_level=0.95, N_bootstrap = 1000, prop_scores=prop.pred, replacement=TRUE)
qte_firpo_zhang_result

# our proposed extremal QTE 
qte_hill_result <- qte_extrapolation_hill(Y=Y, X=NULL, D=D, pn=1-qte_level, ks, CI_level=0.95, prop_scores=prop.pred)
qte_hill_result
