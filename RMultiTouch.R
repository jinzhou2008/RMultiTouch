####################################################################################################################
# RMultiTouch: R Tool For Multi Touch Attribution                                                                  #
#                                                                                                                  #                            
# Objective: Allocate channel attribution                                                                          #
#                                                                                                                  # 
# RMultiTouch(wd, searchlist, user, n, since, until)                                                               #
# Parameters:                                                                                                      #
# wd:  work directory                                                                                              #
####################################################################################################################
 


######################## RMultiTouch READ IN PARAMETERS ######################

RMultiTouch<-function(
wd                     = "C:/jimzh_work/RMultiTouch",
Datafile               = "peak_new_lastday.csv",
Channel                = c("AF_AGG", "PA_AGG","EM_AGG", "NS_AGG", "PS_AGG", "OD_AGG", "SS_AGG"),
Dependent              = "BUY_IND",
Dependent_Binary       = TRUE,
Time                   = "DAYS",
Start                  = "DAYS1", 
Stop                   = "DAYS", 
Output_File            = "RMultiTouch_Output.xlsx"
)

{

setwd(wd)
start.time <- Sys.time()


cat("Total Time Used: ", format(Sys.time()-start.time), "\n\n")
}

