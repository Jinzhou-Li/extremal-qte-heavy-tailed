# This file contains the crossvalidation procedure for selecting the propensity score model.

library(boot)

#### The following section for reading the data has to be uncommented if crossval.R is used on its own and not called by qte_college_on_wage.R
# library(haven)
# 
# df.raw <- read_dta("hhv-neo_v3.dta")
# df.raw <- df.raw[!is.na(df.raw$wage) & (df.raw$collegedeg | df.raw$somecollege | df.raw$hsgrad),]       # select only people which graduaded high-school
# df.nan <- df.raw[,c("wage", "black", "hisp", "south", "west", "northeast", "curban", "broken",
#                                     "age80", "mhgc_mi", "fhgc_mi", "faminc79_th", "numsibs",
#                                     "sasvab1", "sasvab2", "sasvab3", "sasvab4", "sasvab5", "sasvab6" ,
#                                     "sgr9_lang_gpa", "sgr9_math_gpa", "sgr9_scosci_gpa", "sgr9_sci_gpa"
#                                     )]
# df.nan <- cbind(college = (df.raw$collegedeg | df.raw$somecollege), df.nan)                             # college is 1 if participant went to college for some time
# df <- df.nan[complete.cases(df.nan),]                                                                   # remove missing values

####################
controls <- c("black", "hisp", "south", "west", "northeast", "curban", "broken",
		    "age80", "mhgc_mi", "fhgc_mi", "faminc79_th", "numsibs", 
                    "sasvab1", "sasvab2", "sasvab3", "sasvab4", "sasvab5", "sasvab6" ,
                    "sgr9_lang_gpa", "sgr9_math_gpa", "sgr9_scosci_gpa", "sgr9_sci_gpa")
races <- c("black", "hisp")
regions <- c("south", "west", "northeast")

## Log loss 
cost <- function(obs, pred) -mean(obs*log(pred) + (1-obs)*log(1-pred))
set.seed(1)

## CV score for linear model
base.formula <- paste("college~", paste(names(df[,-(1:2)]), collapse="+"))
formula <- as.formula(base.formula)
fit <- glm(formula, data=df, family="binomial")
cv.scores <- cv.glm(df, fit, cost=cost)

## Successively add square terms
cv.min <- cv.scores$delta[1]
n_cont <- length(controls)
print("Square Terms (Continuous only)")
for (i in 8:n_cont){
  print(i)
  name <- controls[i]
  formula.new <- paste(base.formula, "+I(", name, "^2)", sep="")
  formula <- as.formula(formula.new)
  fit <- glm(formula, data=df, family="binomial")
  cv.scores <- cv.glm(df, fit, cost=cost)
  if (cv.scores$delta[1] < cv.min){
    cv.min <- cv.scores$delta[1]
    print(cv.min)
    base.formula <- formula.new 
    print(paste("added ", name, "^2", sep=""))
  }
}

## Successively add interaction terms
print("Interaction Terms")
for (i in 1:(n_cont-1)){
  print(paste("i", i))
  for (j in (i+1):n_cont){
    print(j)
    name1 <- controls[i]
    name2 <- controls[j]
    ## Skip if interaction is not possible
    if (name1 %in% races && name2 %in% races)
	    next
    if (name1 %in% regions && name2 %in% regions)
	    next
    formula.new <- paste(base.formula, "+", name1, ":", name2, sep="")
    formula <- as.formula(formula.new)
    fit <- glm(formula, data=df, family="binomial")
    cv.scores <- cv.glm(df, fit, cost=cost)
      if (cv.scores$delta[1] < cv.min){
        cv.min <- cv.scores$delta[1]
        print(cv.min)
    	  base.formula <- formula.new
        newname <- paste(name1, ":", name2, sep="")
    	print(paste("added", newname))
      }
    }
}

## Final model 
fit <- glm(base.formula, data=df, family="binomial")
summary(fit)
print(base.formula)
print(cv.min)
