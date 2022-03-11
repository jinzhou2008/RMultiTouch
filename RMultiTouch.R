 
####################################################################################################################
# RMultiTouch: R Tool For Multi Touch Attribution                                                                  #
#                                                                                                                  #                            
# Objective: Allocate channel attribution                                                                          #
#                                                                                                                  #
# Methods include:                                                                                                 #
# First Touch                                                                                                      #
# Last Touch                                                                                                       #
# Linear (Even)                                                                                                    #
# Time Decay                                                                                                       #
# Position Based                                                                                                   #
# Probability                                                                                                      #
# Logistic                                                                                                         #
# Coxph                                                                                                            #
#                                                                                                                  #
# Output:                                                                                                          #
# Overall channel attribution                                                                                      #
# Channel attribution for each order                                                                               #
####################################################################################################################


######################## RMultiTouch READ IN PARAMETERS ######################

RMultiTouch<-function(
indir                  = "C:/jimzh_work/RMultiTouch",
outdir                 = "C:/jimzh_work/RMultiTouch",
Datafile               = "peak_new_60days.csv",
Channel                = c("AF", "PA","EM", "NS", "PS", "OD", "SS", "OG"),
Event                  = "BUY_IND",
User_Id                = "VISITOR_DIM_ID",
Time                   = "DAYS",
Half_Life              = 7 ,
Position               = c(0.4, 0.2, 0.4) 
) 


{
 
#start.time  <-  Sys.time()

options(scipen=999) #remove scientific notation in printing 
offsets <- 0.001   # beta + abs(min(beta)) + offsets if negative



require.package<-function(pckg)
{
package.installed<-try(require(pckg, character.only =TRUE))
if (!package.installed) {
cat(paste("Installing", pckg,  "from CRAN\n", sep=" "))
install.packages(pckg,  repos = "http://cran.r-project.org")
# rJava does not work JAVA_HOME is there and not correct
#1. set the PATH environmental variable in my computer
#2. set JAVA_HOME as null
if (Sys.getenv("JAVA_HOME")!="")  Sys.setenv(JAVA_HOME="")
require(pckg, character.only =TRUE)
}#if
}#require.package

require.package("survival")
#require.package("imputation")
require.package("psych")
#require.package("XLConnect")

Channel       <- toupper(Channel)
Event         <- toupper(Event)
User_Id       <- toupper(User_Id)
Time          <- toupper(Time) 



data <- read.csv(Datafile,header=T, stringsAsFactors=F)
dim(data) #4442   26
names(data)<- toupper(names(data))
dim(data)
head(data)
describe(data, skew = FALSE)



#Imputing missing variables
#x<-meanImpute(x)$x

#Remove variables with std=0 
if (any(sapply(data[, Channel],function(x) sd(x)==0))) {
cat(Channel[sapply(data[, Channel],function(x) sd(x)==0)], "has zero variance, please check;", "\n\n")
}
Channel<-Channel[sapply(data[, Channel],function(x) sd(x)!=0)] 



############################### ROUTINE METHODS BEGINS ##################################

#FirstTouch
#LastTouch
#Linear(Even)
#TimeDeacy
#PositionBased

attrbtn_routine <-function(data){

# users with at least one order and order date
data_order <-  data[data[, Event]==1, c(User_Id, Time) ]  #all orders
data_order$Time_At_Lasttouch <- data_order[, Time]
data_order[, Time] <- NULL

#only keep those users with at least one order and days <= order date
data <- merge(data, data_order, by.x=User_Id , by.y=User_Id , all.x=F, all.y=T)
data <- data[ data[, Time] <= data$Time_At_Lasttouch, ]
data <- data[order(data[, User_Id], data$Time_At_Lasttouch, data[, Time]),]
data <- data[!duplicated(data[, c(User_Id, Time)]),]


#first touch day
temp <- aggregate(as.formula(paste0(Time, "~", User_Id, "+", "Time_At_Lasttouch")), data=data, FUN=function(x) min(x,na.rm=T))
temp$Time_At_Firsttouch <- temp[, Time]
temp[, Time] <- NULL


data <- merge(data, temp, by.x=c(User_Id, "Time_At_Lasttouch") , by.y=c(User_Id, "Time_At_Lasttouch"), all.x=T, all.y=F)
data$Time_Until_Lasttouch  <- data$Time_At_Lasttouch  -  data[, Time]
data$Time_From_Firsttouch <- data[, Time] - data$Time_At_Firsttouch 

#order data by user_id, time
data <- data[order(data[, User_Id], data[, Time]),]

####### define function to calulate attribution for each routine method ######
single_method <- function(data=data, Method='_FirstTouch') {
if (Method =='_FirstTouch') {
#weight of that day is 1 if first touch occurs on that day and 0 else
data$Credit <-  ifelse(data$Time_From_Firsttouch==0, 1, 0) 
} else if (Method =='_LastTouch') {
#weight of that day is 1 if last touch occurs on that day and 0 else
data$Credit <-  ifelse(data$Time_Until_Lasttouch==0, 1, 0)
} else if (Method =='_Linear') {
#Linear (Even) method  weight of that day decays as days from last touch goes
data$Credit <- 1
} else if (Method =='_TimeDecay') {
# assign more credit to the interactions which are closest in time to the conversion, half-life t_0.5=10*ln(2)=10*0.693.
#when Half_Life=7, rate=7/0693=10.1
rate <- Half_Life/log(2) 
data$Credit <- exp(-data$Time_Until_Lasttouch/rate)
} else if (Method =='_PositionBased') {
# Assign more credit to the interactions which are closer to first and last touch and less credit to interactions in the middle of the path
data$Credit <-  ifelse(data$Time_From_Firsttouch==0,                                Position[1], 0) 
data$Credit <-  ifelse(data$Time_From_Firsttouch!=0 & data$Time_Until_Lasttouch!=0, Position[2], data$Credit)
data$Credit <-  ifelse(data$Time_Until_Lasttouch==0,                                Position[3], data$Credit)
}

#each day has a decay factor
data[, paste0(Channel, Method)] <- data[, Channel]*data$Credit 
fmla <- paste0("cbind(", paste(Channel, Method, sep='', collapse=','), ")~", User_Id, "+Time_At_Lasttouch")
attrbtn  <- aggregate(as.formula(fmla), data=data, FUN=function(x) sum(x,na.rm=T))
attrbtn <- attrbtn[order(attrbtn[, User_Id], attrbtn$Time_At_Lasttouch),]
head(attrbtn)

attrbtn <- t(apply(as.matrix(attrbtn[, -(1:2)]), 1, FUN=function(x) x/sum(x, na.rm=T)))
attrbtn[is.na(attrbtn)] <- 0 
#attrbtn_individual <- merge(data[data[, Event]==1, !names(data) %in% tempname], attrbtn, by=c(User_Id, "Time_At_Lasttouch"), all.x=T, all.y=F)

tempname <- c("Time_At_Firsttouch",  "Time_Until_Lasttouch", "Time_From_Firsttouch",  "Credit", paste0(Channel, Method))
temp <- data[data[, Event]==1, !names(data) %in% tempname]
temp <- temp[order(temp[, User_Id], temp$Time_At_Lasttouch),]
attrbtn_individual <- cbind(temp, attrbtn)


attrbtn_overall <-  apply(as.matrix(attrbtn), 2, FUN=function(x) sum(x, na.rm=T))
attrbtn_overall <-  attrbtn_overall/sum(attrbtn_overall)

result <- list(attrbtn_overall, attrbtn_individual)

return(result)
} #single_method 

result1 <- single_method(data, '_FirstTouch') 
result2 <- single_method(data, '_LastTouch') 
result3 <- single_method(data, '_Linear') 
result4 <- single_method(data, '_TimeDecay') 
result5 <- single_method(data, '_PositionBased') 

attrbtn_overall <- data.frame(Channel=Channel, FirstTouch=result1[[1]], LastTouch=result2[[1]], Linear=result3[[1]], TimeDecay=result4[[1]], PositionBased=result5[[1]])
row.names(attrbtn_overall) <- Channel

result <- list(attrbtn_overall, 
               attrbtn_individual_1=result1[[2]],
               attrbtn_individual_2=result2[[2]],
               attrbtn_individual_3=result3[[2]],
               attrbtn_individual_4=result4[[2]],
               attrbtn_individual_5=result5[[2]]
               )

return(result)
}
############################### ROUTINE METHODS ENDS #######################################



############################ PROBABILISTIC METHOD BEGINS ##################################
attrbtn_prob <- function(data){

x <- data[, c(Channel), drop=FALSE]  # all channel touches
y <- data[, Event]  # all y s 
n <- length(Channel)

# c(y|x1)= 1/2*p(y|x1) - 1/(2(N-1))*( sum of p(y|xi) for all i's  - sum of p(y|x1,xi) for all i's  ), N= # of channels, inc. 1 itself
# sum of contribution 1/(N-1))* sum of p(y|xi,xj) for all i, j, i^=j  
# p(y|xi)= sum(y*xi)/sum(xi)
# p(y|xi, xj)= sum(y*(xi OR xj)/sum(xi OR xj)


# y <- sample(c(0,1), replace=TRUE, size=5, prob=c(0.5,0.5))
# n <- 3
# x <- matrix(c(1,1,1,1,0,1,1,1,1), nrow=3)

nx <- nrow(x)
p1way <- colSums(y*x)/colSums(x)  #one-channel prob p(y|xi)
p1way[is.infinite(p1way)] <-0

x2 <- apply(x, 2, function(z) x|z)   #two-channel touch xi OR xj
yx2 <- y*x2
sumij  <- apply(x2,  2, function(z) {colSums(matrix(z, nrow=nx))} )
sumyij <- apply(yx2, 2, function(z) {colSums(matrix(z, nrow=nx))} )
p2way <-  sumyij / sumij          #two-channel prob p(y|xi, xj)
p2way[is.infinite(p2way)] <-0

contrbtn <- 1/2*p1way - 1/(2*(n-1))* (sum(p1way) - colSums(p2way)) 
contrbtn  # for channels with smaller marginal prob, contribution may be negative

if (any(contrbtn<=0)) {contrbtn  <- contrbtn + abs(min(contrbtn)) + offsets}  #offsets negative contribution
attrbtn <- abs(contrbtn)/sum(abs(contrbtn))

attrbtn_overall <- data.frame(Channel=Channel, p1way=p1way, p2way=p2way, contrbtn=contrbtn, attrbtn=attrbtn)



x<-data[, Channel, drop=FALSE]
attrbtn <- matrix(rep(attrbtn, nrow(x)), nrow=nrow(x), byrow=T)
attrbtn <- attrbtn* x
attrbtn <- t(apply(as.matrix(attrbtn), 1, FUN=function(x) x/sum(x)))
colnames(attrbtn) <- paste0(Channel, "_Attrbtn")

attrbtn_individual <- cbind(data, attrbtn)
#only show attrbtn of buyers
attrbtn_individual <- attrbtn_individual[attrbtn_individual[,Event]==1,]
attrbtn[is.na(attrbtn)] <- 0 
describe(attrbtn_individual, skew = FALSE)

result <- list(attrbtn_overall, attrbtn_individual)

return(result)

}
############################ PROBABILISTIC METHOD ENDS ##################################




############################ LINEAR REGRESSION BEGINS ##################################
########################## IF REVENUE is availble ###############
attrbtn_linear <-function(data){
fmla <- paste0(Revenue, "~", paste0(Channel, collapse="+"))
fit <- glm(fmla, data=data,  family =gaussian)
summary(fit)
names(fit)
names(summary(fit))
names(anova(fit))
y <- fit$model[,1]
x <- fit$model[-1]
b <- fit$coef[-1]
yhat <- predict(fit)  
ux <- sapply(x, mean)
uy <- mean(y)
sx <- sapply(x, sd)
sy <- sd(y)
rxy <- cor(x, y);   # dimnames(rxy)<-NULL #if y <- fit$model[1]
R2 <- cor(y, yhat)^2  #R2 <- 1-(sum((y-yhat)^2)/sum((y-uy)^2))   #R2 <- summary(fit)$r.squared  #if fit <- lm(fmla, data=data)
beta <- b*sx/sy
if (any(beta<=0)) {beta <- beta + abs(min(beta, na.rm=T)) + offsets}  #offsets negative beta
#attrbtn <- beta*rxy       # b*sx*rxy/sy
attrbtn <- beta*ux/sx    # b*ux/sy
attrbtn <- abs(attrbtn)/sum(abs(attrbtn), na.rm=T)

attrbtn_overall <- data.frame(Channel=Channel, b=b, ux=ux, sx=sx, sy=sy, beta=beta, R2=R2, rxy=rxy, attrbtn=attrbtn)


#yhat_newdata <- predict.lm(fit, newdata=newdata, type="terms", terms=Channel)
#head(yhat_newdata)
x <- data[, Channel, drop=FALSE]
sx <- sapply(x, sd) 
sx <- matrix(rep(sx, nrow(x)), nrow=nrow(x), byrow=T)
beta <- matrix(rep(beta, nrow(x)), nrow=nrow(x), byrow=T)
attrbtn <- beta * x/sx
attrbtn <- abs(attrbtn)
attrbtn <- t(apply(as.matrix(attrbtn), 1, FUN=function(x) x/sum(x, na.rm=T)))
colnames(attrbtn) <- paste0(Channel, "_Attrbtn")

attrbtn_individual <- cbind(data, attrbtn)
#only show attrbtn of buyers
attrbtn_individual <- attrbtn_individual[attrbtn_individual[,Event]==1,]
attrbtn[is.na(attrbtn)] <- 0 
describe(attrbtn_individual, skew = FALSE)

result <- list(attrbtn_overall, attrbtn_individual)

return(result)
}

############################ LINEAR REGRESSION ENDS ##################################


############################ LOGISTIC REGRESSION BEGINS ##################################
attrbtn_logistic<-function(data){
fmla <- paste0(Event, "~", paste0(Channel, collapse="+"))
fit <- glm(fmla, data=data,  family = binomial(link = "logit"))
summary(fit)
names(fit)
names(summary(fit))
names(anova(fit))
y <- fit$model[,1]
x <- fit$model[-1]
b <- fit$coef[-1]
yhat <- predict(fit, type="response") #log odds ratio= log(p)/(1-log(p)) #predict(fit, type=c("link", "response", "terms")) 
ux <- sapply(x, mean)
uy <- mean(y)
sx <- sapply(x, sd)
sy <- sd(y)
rxy <- cor(x, y)   # dimnames(rxy)<-NULL #if y <- fit$model[1]
R2 <- cor(y, yhat)^2  #R2 <- 1-(sum((y-yhat)^2)/sum((y-uy)^2))   #R2 <- summary(fit)$r.squared  #if fit <- lm(fmla, data=data)
beta <- b*sx/sy
if (any(beta<=0)) {beta <- beta + abs(min(beta, na.rm=T)) + offsets}  #offsets negative beta
#attrbtn <- beta*rxy       # b*sx*rxy/sy
attrbtn <- beta*ux/sx    # b*ux/sy
attrbtn <- abs(attrbtn)/sum(abs(attrbtn), na.rm=T)



attrbtn_overall <- data.frame(Channel=Channel, b=b, ux=ux, sx=sx, sy=sy, beta=beta, R2=R2, rxy=rxy, attrbtn=attrbtn)


#yhat_newdata <- predict.lm(fit, newdata=newdata, type="terms", terms=Channel)
#head(yhat_newdata)
x <- data[, Channel, drop=FALSE]
sx <- sapply(x, sd) 
sx <- matrix(rep(sx, nrow(x)), nrow=nrow(x), byrow=T)
beta <- matrix(rep(beta, nrow(x)), nrow=nrow(x), byrow=T)
attrbtn<- beta * x/sx
attrbtn <- abs(attrbtn)
attrbtn <- t(apply(as.matrix(attrbtn), 1, FUN=function(x) x/sum(x, na.rm=T)))
colnames(attrbtn) <- paste0(Channel, "_Attrbtn")

attrbtn_individual <- cbind(data, attrbtn)
#only show attrbtn of buyers
attrbtn_individual <- attrbtn_individual[attrbtn_individual[,Event]==1,]
describe(attrbtn_individual, skew = FALSE)

result <- list(attrbtn_overall, attrbtn_individual)

return(result)

}

############################ LOGISTIC REGRESSION ENDS ##################################




############################ COX REGRESSION BEGINS ##################################
attrbtn_coxph<-function(data){
fmla <- paste("Survtime ~", paste(Channel, collapse="+"))
Survtime <- Surv(data[,Time], data[,Event])
#Survtime <- Surv(data[,Start], data[,Stop], data[,Event])
#coxph(Surv(start, stop, event) ~ x + sex, robust=FALSE,  model=FALSE, x=FALSE, y=TRUE)
fit <- coxph(as.formula(fmla), data, ties=c("efron"), model=TRUE, x=TRUE, y=TRUE)
summary(fit)

y <- fit$model[,1][,1]   #Only take time from (time, status)
x <- fit$model[-1]
b <- fit$coef   #no intercept
yhat <- predict(fit, type="risk") #hazard # predict(fit, type=c("lp", "risk", "expected", "terms")
ux <- sapply(x, mean)
uy <- mean(y)
sx <- sapply(x, sd)
sy <- sd(y)
rxy <- cor(x, y);   # dimnames(rxy)<-NULL #if y <- fit$model[1]
R2 <- cor(y, yhat)^2  #R2 <- 1-(sum((y-yhat)^2)/sum((y-uy)^2))   #R2 <- summary(fit)$r.squared  #if fit <- lm(fmla, data=data)
beta <- b*sx/sy
if (any(beta<=0)) {beta <- beta + abs(min(beta, na.rm=T)) + offsets}  #offsets negative beta
#attrbtn <- beta*rxy       # b*sx*rxy/sy
attrbtn <- beta*ux/sx    # b*ux/sy
attrbtn <- abs(attrbtn)/sum(abs(attrbtn), na.rm=T)

attrbtn_overall <- data.frame(Channel=Channel, b=b, ux=ux, sx=sx, sy=sy, beta=beta, R2=R2, rxy=rxy, attrbtn=attrbtn)


#yhat_newdata <- predict.lm(fit, newdata=newdata, type="terms", terms=Channel)
#head(yhat_newdata)
x <- data[, Channel, drop=FALSE]
sx <- sapply(x, sd) 
sx <- matrix(rep(sx, nrow(x)), nrow=nrow(x), byrow=T)
beta <- matrix(rep(beta, nrow(x)), nrow=nrow(x), byrow=T)
attrbtn<- beta * x/sx
attrbtn <- abs(attrbtn)
attrbtn <- t(apply(as.matrix(attrbtn), 1, FUN=function(x) x/sum(x, na.rm=T)))
colnames(attrbtn) <- paste0(Channel, "_Attrbtn")

attrbtn_individual <- cbind(data, attrbtn)
#only show attrbtn of buyers
attrbtn_individual <- attrbtn_individual[attrbtn_individual[,Event]==1,]

result <- list(attrbtn_overall, attrbtn_individual)

return(result)

}

############################ COX REGRESSION ENDS ##################################






############################ SCORE DATA BEGINS ##################################
#Routine
if (!is.na(Time)) {
routine_attr   <- attrbtn_routine(data)
}#if (!is.na(Time))

#Simple
prob_attr <- attrbtn_prob(data)

#Linear
if (FALSE) {
linear_attr   <- attrbtn_linear(data)
}

#Logistic
logistic_attr <- attrbtn_logistic(data)

#Cox
coxph_attr   <- attrbtn_coxph(data)

overall_attr <- data.frame(routine_attr[[1]], Probability=prob_attr[[1]]$attrbtn, Logistic=logistic_attr[[1]]$attrbtn, Coxph=coxph_attr[[1]]$attrbtn)
############################ SCORE DATA ENDS ##################################



################################ Output ##################################
if (FALSE) {
if (file.exists("Output.xlsx")) file.remove("Output.xlsx")

wb <- loadWorkbook("Output.xlsx", create = TRUE)

sheet <- "Overall_Attrbution"
createSheet(wb, name = sheet)
writeWorksheet(wb, overall_attr,       sheet=sheet , startRow = 1 , startCol = 1, header = TRUE )
setColumnWidth(wb, sheet=sheet, column = 1:10, width = rep(4000,10))  #1/256

sheet <- "Indiv_FirstTouch"
createSheet(wb, name = sheet)
writeWorksheet(wb, routine_attr[[2]],  sheet=sheet , startRow = 1 , startCol = 1, header = TRUE )

sheet <- "Indiv_LastTouch"
createSheet(wb, name = sheet)
writeWorksheet(wb, routine_attr[[3]],  sheet=sheet , startRow = 1 , startCol = 1, header = TRUE )

sheet <- "Indiv_Linear"
createSheet(wb, name = sheet)
writeWorksheet(wb, routine_attr[[4]],  sheet=sheet , startRow = 1 , startCol = 1, header = TRUE )

sheet <- "Indiv_TimeDecay"
createSheet(wb, name = sheet)
writeWorksheet(wb, routine_attr[[5]],  sheet=sheet , startRow = 1 , startCol = 1, header = TRUE )

sheet <- "Indiv_PositionBased"
createSheet(wb, name = sheet)
writeWorksheet(wb, routine_attr[[6]],  sheet=sheet , startRow = 1 , startCol = 1, header = TRUE )

sheet <- "Indiv_Probability"
createSheet(wb, name = sheet)
writeWorksheet(wb, prob_attr[[2]],  sheet=sheet , startRow = 1 , startCol = 1, header = TRUE )


sheet <- "Indiv_Logistic"
createSheet(wb, name = sheet)
writeWorksheet(wb, logistic_attr[[2]],  sheet=sheet , startRow = 1 , startCol = 1, header = TRUE )

sheet <- "Indiv_Coxph"
createSheet(wb, name = sheet)
writeWorksheet(wb, coxph_attr[[2]],  sheet=sheet , startRow = 1 , startCol = 1, header = TRUE )


saveWorkbook (wb)
}

write.csv(overall_attr,         paste0(outdir, "/Overall_Attrbution.csv"),   row.names=F)
write.csv(routine_attr[[2]],    paste0(outdir, "/Indiv_FirstTouch.csv"),     row.names=F)
write.csv(routine_attr[[3]],    paste0(outdir, "/Indiv_LastTouch.csv"),      row.names=F)
write.csv(routine_attr[[4]],    paste0(outdir, "/Indiv_Linear.csv"),         row.names=F)
write.csv(routine_attr[[5]],    paste0(outdir, "/Indiv_TimeDecay.csv"),      row.names=F)
write.csv(routine_attr[[6]],    paste0(outdir, "/Indiv_PositionBased.csv"),  row.names=F)
write.csv(prob_attr[[2]],       paste0(outdir, "/Indiv_Probability.csv"),    row.names=F)
write.csv(logistic_attr[[2]],   paste0(outdir, "/Indiv_Logistic.csv"),       row.names=F)
write.csv(coxph_attr[[2]],      paste0(outdir, "/Indiv_Coxph.csv"),          row.names=F)


################################ Output ##################################


#cat("Total Time Used: ", format(Sys.time()-start.time), "\n\n")


}#RMultiTouch


