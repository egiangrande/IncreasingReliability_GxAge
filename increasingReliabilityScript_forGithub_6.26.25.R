#Fitting continuous time dynamic models of changes in twin correlations across cognitive development
#using ctsem. These models allow for estimation of G x Age interaction, and for comparison
#of reliability-corrected estimates to uncorrected estimates. 

#-See the Giangrande et al. (2025, preprint) manuscript and supplement for model details. 

#----------------------------------------------------------------------------------------
library(ctsem)
library(data.table)
library(ggplot2)
library(abind)

### Read-in cleaned data and prepare for analysis:
ltsData <- #
names(ltsData)[names(ltsData)=="momId"] <- "id"
names(ltsData)[names(ltsData)=="ageYrs"] <- "time"

#Note: Zygosity is coded as mz = -1, dz = 1

#Rescale the cognitive existing scores by subtracting test mean and then dividing by test SD,
#separately for each test.

#Bayley Scales of Infant Development, Second Edition, mean 100, SD 15
#WPPSI & WISC: mean 100, SD 15
#Stanford-Binet (SB) 1972 norms: mean 100, SD 16
#McCarthy General cognitive index: mean 100, SD 16

ltsData$adjValue1 <- ifelse(ltsData$test=="bayley"|ltsData$test=="wppsi"|ltsData$test=="wisc", (ltsData$value.1-100)/15, 
                            ifelse(ltsData$test=="sb"|ltsData$test=="mccarthy", (ltsData$value.1-100)/16,
                                   NA))
ltsData$adjValue2 <- ifelse(ltsData$test=="bayley"|ltsData$test=="wppsi"|ltsData$test=="wisc", (ltsData$value.2-100)/15, 
                            ifelse(ltsData$test=="sb"|ltsData$test=="mccarthy", (ltsData$value.2-100)/16,
                                   NA))

#### Adding binary test variables
ltsData$test1 <- ifelse(ltsData$test=="bayley",1,0)
ltsData$test2 <- ifelse(ltsData$test=="sb",1,0)
ltsData$test3 <- ifelse(ltsData$test=="mccarthy",1,0)
ltsData$test4 <- ifelse(ltsData$test=="wppsi",1,0)
ltsData$test5 <- ifelse(ltsData$test=="wisc",1,0)

#### Add phantom time variable to every twin pair to account for different starting ages.
#To do so, create phantom 3mo. variable for people with missing data at age 3 mos. 
have3mo <-subset(ltsData, time==0.25) 
haveIds <- have3mo$id
phant <- subset(ltsData, !(id %in% haveIds) &  time!=0.25)
phant <- phant[!duplicated(phant$id),]
phant[, c("variable","value.1","value.2","adjValue1","adjValue2",
          "test","test1","test2","test3","test4","test5"
)] <- NA
phant$time <- .25
ltsData <- rbind(ltsData, phant)
ltsData <- ltsData[order(ltsData$id, ltsData$time), ]

ltsData <- data.table(ltsData)
ltsData$id <- factor(ltsData$id)

###Divide time by 10 to help with numerical stability 
ltsData$time <- ltsData$time/10

#------------------------------------------------------------------------------------------------------
#Model specification
latentNames=c('level1',"level2",'slope1','slope2',"agetrend")
  
basecor <- '2/(1 + exp(-basecor)) - 1'
diffcor <- "2/(1 + exp(-(basecor + diffcor + diffcorByAge*agetrend + 
                diffcorByAgeSq*(agetrend)^2))) - 1"
diffsd <- "log1p(exp(diffsd + diffsdbyage*agetrend))+1e-5"
drift <- "-log1p(exp(-(drift+driftbyage*agetrend)))-1e-6"

  T0VAR = matrix(c(
    "t0levelsd",0,0,0,0,
    basecor,"t0levelsd", 0,0,0,
    0,0,'baseslopesd',0,0,
    0,0,diffcor,'baseslopesd',0,
    0,0,0,0,0),byrow=TRUE,5,5)
  
  DRIFT=matrix(c(
    0,0,1,0,0,
    0,0,0,1,0,
    'levelonslope',0,drift,0,0,
    0,'levelonslope',0,drift,0,
    0,0,0,0,"trenddrift"),byrow=TRUE,5,5)
  
  DIFFUSION=matrix(c(
    0,0,0,0,0,
    0,0,0,0,0,
    0,0,diffsd,0,0,
    0,0,diffcor,diffsd,0,
    0,0,0,0,0),byrow=TRUE,5,5)
  
  CINT=c(0,0,'cint1||FALSE',"cint1||FALSE",1) #Disable individual variation 
  
  T0MEANS=c('baselevel||FALSE',"baselevel||FALSE",'baseslope||FALSE','baseslope||FALSE',0) #Disable individual variation 
  
  LAMBDA = matrix(c(
    1,0,0,0,0,
    0,1,0,0,0),byrow=TRUE,2,5)

#---------------------------------------------------------------------------------------------
#Adjust reliability correction. Do so by adjusting how measurement error is specified in 
# the measurement model, and then updating the PARS block accordingly. 
  
#Below is the specification for each model.
#Run this script separately for each model after commenting out the appropriate blocks.
  
  #1) Full model (correcting for reliability changes across age and test battery)
  ressd <- "log1p_exp(ressdByAge*agetrend + ressdByTest1*tdpreds[rowi,1] +
          ressdByTest2*tdpreds[rowi,2] + ressdByTest3*tdpreds[rowi,3] +
            ressdByTest4*tdpreds[rowi,4] + ressdByTest5*tdpreds[rowi,5])"
  
  ressdCor <- '2/(1 + exp(-ressdCor + ressdCorByAge*agetrend)) - 1'  
  
  PARS=c('basecor ||||zyg',
           'diffcor|||| zyg',
           "diffcorByAge|||| zyg","diffcorByAgeSq|||| zyg",
           "diffsd","diffsdbyage",
           'drift','driftbyage',
           "ressdByAge", "ressdByTest1", "ressdByTest2",
           "ressdByTest3", "ressdByTest4", "ressdByTest5",
           "ressdCor",
           "ressdCorByAge")
  
  #2) Test-corrected model (correcting for reliability changes across test, but not age)
    # ressd <- "log1p_exp(ressdByTest1*tdpreds[rowi,1] +
    #             ressdByTest2*tdpreds[rowi,2] + ressdByTest3*tdpreds[rowi,3] +
    #                    ressdByTest4*tdpreds[rowi,4] + ressdByTest5*tdpreds[rowi,5])"
    # 
    # ressdCor <- '2/(1 + exp(-ressdCor)) - 1'
    # 
    # PARS=c(PARS,
    #        "diffcorByAge|||| zyg","diffcorByAgeSq|||| zyg",
    #        "diffsd","diffsdbyage",
    #        'drift','driftbyage',
    #        "ressdByTest1", "ressdByTest2",
    #        "ressdByTest3", "ressdByTest4", "ressdByTest5",
    #        "ressdCor")
  
  #3) Age-corrected model (correcting for reliability changes across age, but not test)
      # ressd <- "log1p_exp(ressdByAge*agetrend)" 
      # 
      # ressdCor <- '2/(1 + exp(-ressdCor + ressdCorByAge*agetrend)) - 1'
      # 
      # PARS=c(PARS,
      #        "diffcorByAge|||| zyg","diffcorByAgeSq|||| zyg",
      #        "diffsd","diffsdbyage",
      #        'drift','driftbyage',
      #        "ressdByAge", 
      #        "ressdCor",
      #        "ressdCorByAge")
  
  #4) Uncorrected model 
    # ressd <- "ressd"
    # 
    # ressdCor <- '2/(1 + exp(-ressdCor)) - 1'
    # 
    # PARS=c(PARS,
    #      "diffcorByAge|||| zyg","diffcorByAgeSq|||| zyg",
    #      "diffsd","diffsdbyage",
    #      'drift','driftbyage',
    #      "ressdCor")

  #------------------------------------------------------  
  
  m1 <- ctModel(type = 'stanct',tipredDefault = FALSE,
                manifestNames = c('adjValue1',"adjValue2"), 
                latentNames=latentNames,
                TIpredNames = c('zyg'),
                TDpredNames = c("test1","test2","test3","test4","test5"),
                T0VAR = T0VAR,
                DRIFT=DRIFT,
                DIFFUSION=DIFFUSION,
                CINT=CINT,
                T0MEANS=T0MEANS,
                TDPREDEFFECT = 0,
                LAMBDA = LAMBDA,
                MANIFESTMEANS=matrix(c('td1*tdpreds[rowi,1]+td2*tdpreds[rowi,2]+td3*tdpreds[rowi,3]+td4*tdpreds[rowi,4]+td5*tdpreds[rowi,5]',
                                       'td1*tdpreds[rowi,1]+td2*tdpreds[rowi,2]+td3*tdpreds[rowi,3]+td4*tdpreds[rowi,4]+td5*tdpreds[rowi,5]'),2,1),
                MANIFESTVAR=c(ressd,0,
                             "ressdCor",ressd), 
                PARS = c(PARS,
                         "td1", "td2", "td3", "td4", "td5")
  )
  
  m1$pars$sdscale[!grepl('drift',m1$pars$param)] <- 10

  #--------------------------------------------------------------------
  # Fit the model
  cores = 15 #adjust depending on compute setup/capacity
  
  f <- ctStanFit(datalong = ltsData, ctstanmodel = m1,cores=cores,verbose=1,
                 nopriors = F, 
                 plot = 10, 
                 saveComplexPars = T) 
  options(max.print=1000000)
  
  summary(f)
  
  save(f, file = paste0("./fitCtsemModel_f.RData"))
  
  z <- ltsData$id[c(match(-1,ltsData$zyg),match(1,ltsData$zyg))]
  
  f$stanfit$transformedparsfull$subj_T0cov[1,z[1],,] 
  cov2cor(f$stanfit$transformedparsfull$subj_T0cov[1,z[1],,])
  cov2cor(f$stanfit$transformedparsfull$subj_DIFFUSIONcov[1,z[1],,])
  
  f$stanfit$transformedparsfull$subj_T0cov[1,z[2],,] 
  cov2cor(f$stanfit$transformedparsfull$subj_T0cov[1,z[2],,])
  cov2cor(f$stanfit$transformedparsfull$subj_DIFFUSIONcov[1,z[2],,])
  
  #--------------------------------------------------------------------------------------------------------
  #Plot parameter change over time, separately for each zygosity:
  
  rstan::expose_stan_functions(ctsem:::stanmodels$ctsm)
  e <- ctExtract(f,subjectMatrices = T,cores = cores,nsamples = 500) #Can set nsamples to low value for fast approximate version
  nlpars <- e$calcs
  dimnames(nlpars) <- list(
    Iteration=1:dim(nlpars)[1],
    Row=1:dim(nlpars)[2],
    Parameter= c(f$ctstanmodel$modelmats$calcs$calcNames))
  
  nlT0VAR <- e$pop_T0VAR[1,,] 
  nlDIFFUSION <- e$pop_DIFFUSION[1,,]
  nlcors <- array(NA,c(dim(nlpars)[1:2],3))
  for(i in 1:dim(nlcors)[1]){
    print(i)
    for(j in 1:dim(nlcors)[2]){
      nlT0VAR[2,1] <- nlpars[i,j,"T0VAR[2, 1]"]
      nlT0VAR[4,3] <- nlpars[i,j,"T0VAR[4, 3]"]
      nlDIFFUSION[4,3] <- nlpars[i,j,"DIFFUSION[4, 3]"]
      
      nlcors[i,j,] <- c(
        cov2cor(sdcovsqrt2cov(mat = nlT0VAR,choleskymats = 0L))[cbind(c(2,4),c(1,3))]
        ,cov2cor(sdcovsqrt2cov(mat = nlDIFFUSION,choleskymats = 0L))[cbind(c(4),c(3))]
      )
    }
  }
  
  nlpars <- abind(nlpars,nlcors)
  dimnames(nlpars) <- list(
    Iteration=1:dim(nlpars)[1],
    Row=1:dim(nlpars)[2],
    Parameter= c(f$ctstanmodel$modelmats$calcs$calcNames,'t0cor21','t0cor43','diffusion43'))
  
  nlextra <- data.table(Age=f$standata$time, Subject=f$standata$subject,Row=1:f$standata$ndatapoints)
  nltipreds <- data.table(Zygosity=factor(f$standata$tipredsdata[,1]),Subject=1:f$standata$nsubjects)
  nlpars <- as.data.table(nlpars)
  nlpars[,Row:=as.numeric(Row)]
  nlpars <- merge(nlpars, nlextra,by = 'Row')
  nlpars <- merge(nlpars, nltipreds,by = 'Subject')

  print(ggplot(nlpars,aes(y=value,x=Age,fill=Zygosity, colour=Zygosity))+
    stat_summary(fun.data=mean_sdl,fun.args=list(mult=2,na.rm=TRUE), geom="ribbon", alpha=0.2)+
    stat_summary(fun=mean,geom='point',size=1)+
    theme_bw()+
    facet_wrap(vars(Parameter), scales = "free_y"))
  
  #------------------------------------------------------------------------------------------
  #Plot estimated latent twin correlations over time by zygosity
  
  k1=ctKalman(f,timestep=.1,removeObs = TRUE,subjects = z[1]) #compute expectations over time, independent of observations but dependent on covariates
  k0=ctKalman(f,timestep=.1,removeObs = TRUE,subjects = z[2])
  k1=data.table(k1);k0=data.table(k0)
  
  plot(  k1[Element == 'etapriorcov' & Col =='level1' & Row=='level2',Time],
         k1[Element == 'etapriorcov' & Col =='level1' & Row=='level2',value] / (
           sqrt(k1[Element == 'etapriorcov' & Col =='level1' & Row=='level1',value]) * 
             sqrt(k1[Element == 'etapriorcov' & Col =='level2' & Row=='level2',value])),
         ylab='Corr',ylim=c(0,1),xlab='Time',type='l',col='black',main='Latent Correlations')
  
  points(  k0[Element == 'etapriorcov' & Col =='level1' & Row=='level2',Time],
           k0[Element == 'etapriorcov' & Col =='level1' & Row=='level2',value] / (
             sqrt(k0[Element == 'etapriorcov' & Col =='level1' & Row=='level1',value]) * 
               sqrt(k0[Element == 'etapriorcov' & Col =='level2' & Row=='level2',value])),
           ylab='Corr',ylim=c(0,1),xlab='Time',type='l',col='red')

  legend('bottomleft',legend = c('Monozygotic','Dizygotic'),
         lty=c(1,1),col=c(1,2),pch=c(NA,NA),bty='n')
  
 #------------------------------------------------------------------------------------------
  #Plot model-expected raw twin correlations vs. observed raw twin correlations by zygosity.
  #Expected correlations are plotted as solid lines, observed are dashed lines. 
  
  #Calculate twin corrs for observed cognitive data by zygosity, across ages
  ltsData[, "obscors":= list(cor(adjValue1,adjValue2, use = "pairwise.complete.obs")), by=list(time,zyg)] 
  
  #To facilitate plotting, extract information from participants who have complete data (i.e., 
  #non-missing observations at each time point).
  completeObs <- subset(ltsData, id == idnumberX | id == idnnumberY)
  completeObs <- subset(completeObs, select = c(id, zyg, time, obscors))
  
  plot(  k1[Element == 'ypriorcov' & Col =='adjValue1' & Row=='adjValue2',Time],
         k1[Element == 'ypriorcov' & Col =='adjValue1' & Row=='adjValue2',value] / (
           sqrt(k1[Element == 'ypriorcov' & Col =='adjValue1' & Row=='adjValue1',value]) * 
             sqrt(k1[Element == 'ypriorcov' & Col =='adjValue2' & Row=='adjValue2',value])),
         ylab='Correlation',ylim=c(0,1),xlab='Age (Years)', col="black", type='l')
  
  points(  k0[Element == 'ypriorcov' & Col =='adjValue1' & Row=='adjValue2',Time],
           k0[Element == 'ypriorcov' & Col =='adjValue1' & Row=='adjValue2',value] / (
             sqrt(k0[Element == 'ypriorcov' & Col =='adjValue1' & Row=='adjValue1',value]) * 
               sqrt(k0[Element == 'ypriorcov' & Col =='adjValue2' & Row=='adjValue2',value])),
           ylab='Corr',ylim=c(0,1),xlab='Time',type='l',col='red')
  
  points(completeObs[id== id[which(zyg==1)[1]],time],completeObs[id==id[which(zyg==1)[1]],obscors],ylim=c(0,1),type='l',xlab='Time',ylab='Corr.',lty=2,col=2)
  points(completeObs[id== id[which(zyg==-1)[1]],time],completeObs[id==id[which(zyg==-1)[1]],obscors],ylim=c(0,1),type='l',lty=2)
  
  legend('bottomleft',legend = c('Monozygotic','Dizygotic'),
         lty=c(1,1),col=c(1,2),pch=c(NA,NA),bty='n')
  
