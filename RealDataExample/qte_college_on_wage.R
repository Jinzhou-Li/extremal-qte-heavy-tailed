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

### Generate categorica variables 'race' and 'region of residence' based on 
# ('black', 'hisp') and ('south', 'west', 'northeast'), respectively.
race <- rep("OtherRace", nrow(df))
race[df[,3]==1] <- "black"
race[df[,4]==1] <- "hisp"
df <- data.frame(df, race=as.factor(race))

region <- rep("OtherRegion", nrow(df))
region[df[,5]==1] <- "south"
region[df[,6]==1] <- "west"
region[df[,7]==1] <- "northeast"
df <- data.frame(df, region=as.factor(region))

df <- df[,-c(3:7)]

### Estimate extreme QTE
Y = df$wage
D = df$college

qte_level <- 0.999
ks <- 85

### Fit a full logistic regression model for estimating propensity scores
fit <- glm(as.formula(paste0("college~",paste(names(df)[-c(1,2)], collapse = '+'))),
           data=df, family="binomial")
prop.pred <- as.numeric(fitted(fit))

# zhang's b out of n bootstrap
set.seed(2021)
qte_firpo_zhang_result <- qte_firpo_zhang(Y, X=NULL, D, pn=1-qte_level, CI_level=0.9, N_bootstrap = 1000, prop_scores=prop.pred, replacement=TRUE)
qte_firpo_zhang_result

# our proposed extremal QTE 
qte_hill_result <- qte_extrapolation_hill(Y=Y, X=NULL, D=D, pn=1-qte_level, ks, CI_level=0.9, prop_scores=prop.pred)
qte_hill_result
