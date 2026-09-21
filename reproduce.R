suppressPackageStartupMessages(library(survival))
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
stopifnot(length(script_arg)==1)
root <- dirname(sub("^--file=", "", script_arg))
w <- file.path(root, "inputs")
mode <- if(length(args)>0) args[1] else "smoke"
stopifnot(mode %in% c("smoke", "full"))
out_dir <- if(length(args)>1) args[2] else file.path(root, paste0("run_",mode))
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
cat("Mode:",mode," Output:",out_dir,"\n")
defs <- c('OX3','OX4','OX5','OX5_noNRF2','OX6')
extra <- c('ferroptosis_c','apoptosis_c','necroptosis_c','parthanatos_c','mitophagy_c','oxstress_c','S_prolif_c','S_mito_c')
fits <- list(); bgrows <- list(); phrows <- list()
fitone <- function(d,sc,cv,time='OS_time',ev='OS_event',ph=FALSE) {
  d$stage4 <- droplevels(factor(d$stage4,levels=c('I','II','III','IV','Unknown')))
  cv <- unique(cv); warnings <- character()
  f <- as.formula(paste0('Surv(',time,',',ev,') ~ ',paste(c(sc,cv),collapse='+')))
  fit <- tryCatch(withCallingHandlers(coxph(f,d,ties='efron',x=TRUE),warning=function(e){warnings <<- c(warnings,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(fit,'error'))return(data.frame(HR=NA,lo95=NA,hi95=NA,p=NA,n=nrow(d),events=sum(d[[ev]]),parameters=NA,iter=NA,score_residual_norm=NA,warning=conditionMessage(fit),ph_global=NA))
  s<-summary(fit);i<-which(rownames(s$coefficients)==sc)
  pv<-if(ph) tryCatch(cox.zph(fit,transform='rank',terms=TRUE)$table['GLOBAL','p'],error=function(e) NA_real_) else NA_real_
  data.frame(HR=s$conf.int[i,1],lo95=s$conf.int[i,3],hi95=s$conf.int[i,4],p=s$coefficients[i,5],n=fit$n,events=fit$nevent,parameters=sum(!is.na(coef(fit))),iter=fit$iter,score_residual_norm=sqrt(sum(colSums(as.matrix(residuals(fit,type='score')))^2)),warning=paste(warnings,collapse=' | '),ph_global=pv)
}
for(name in c('TCGA','GSE39582','GSE17538')){
  co<-if(name=='TCGA')'TCGA-COADREAD' else name
  d<-read.csv(file.path(w,paste0(name,'_model_inputs.csv')),check.names=FALSE)
  for(sc in defs)for(m in c('M0','M1','M2')){
    cv<-switch(m,M0=character(),M1=c('age','male','stage4'),M2=c('age','male','stage4',extra))
    rr<-fitone(d,sc,cv,ph=(sc=='OX4'));rr$cohort<-co;rr$score<-sc;rr$endpoint<-'OS';rr$variant<-'frozen';rr$model<-m;fits[[length(fits)+1]]<-rr
  }
  if(mode=='full'){
  ep<-read.csv(file.path(w,paste0(name,'_endpoint_inputs.csv')),check.names=FALSE)
  vv<-list(V0_old=ep,V1_drop_timeLe0=ep[ep$ep_time>0,],V2_drop_timeLe0_noIV=ep[ep$ep_time>0 & ep$stage4!='IV',],V2b_primary=ep[ep$ep_time>0 & !ep$stage4 %in% c('IV','Unknown'),],V3_src_stageII_III=ep[ep$stage4 %in% c('II','III'),],V4_src_stageII_III_timeGt0=ep[ep$stage4 %in% c('II','III') & ep$ep_time>0,])
  if(name!='TCGA'){
    dd<-vv$V4_src_stageII_III_timeGt0;dd$ep_event[dd$ep_time>60]<-0;dd$ep_time<-pmin(dd$ep_time,60);vv[['V5_maxfidelity_cens60']]<-dd
  }
  for(v in names(vv))for(sc in defs)for(m in c('M0','M1','M2')){
    cv<-switch(m,M0=character(),M1=c('age','male','stage4'),M2=c('age','male','stage4',extra))
    rr<-fitone(vv[[v]],sc,cv,'ep_time','ep_event',ph=(sc=='OX4' && m=='M1' && v=='V5_maxfidelity_cens60'));rr$cohort<-co;rr$score<-sc;rr$endpoint<-'recurrence_diagnostic';rr$variant<-v;rr$model<-m;fits[[length(fits)+1]]<-rr
  }
  b<-read.csv(file.path(w,paste0(name,'_random_scores.csv.gz')),check.names=FALSE)
  reg<-read.csv(file.path(w,paste0(name,'_random_registry.csv')),check.names=FALSE)
  stopifnot(identical(d$sample_id,b$sample_id))
  bd<-d[,c('OS_time','OS_event','age','male','stage4')]
  for(i in seq_len(nrow(reg))){
    bd$BG<-b[[reg$audit_score[i]]]
    rr<-fitone(bd,'BG',c('age','male','stage4'));rr$cohort<-co;rr$design<-reg$design[i];rr$set_id<-reg$set_id[i];rr$upstream_logHR<-reg$logHR_M1[i];rr$logHR_error<-abs(log(rr$HR)-rr$upstream_logHR);bgrows[[length(bgrows)+1]]<-rr
    if(i%%500==0)cat(co,i,'random M1 fits completed\n')
  }
  }
  write.csv(do.call(rbind,fits),file.path(out_dir,'R_observed_and_endpoint_fits.csv'),row.names=FALSE)
  if(length(bgrows)>0)write.csv(do.call(rbind,bgrows),file.path(out_dir,'R_all6000_matched_M1.csv'),row.names=FALSE)
}
capture.output(sessionInfo(),file=file.path(out_dir,'R_session_v04.txt'))
cat('Completed',length(fits),'observed/diagnostic and',length(bgrows),'random fits\n')
