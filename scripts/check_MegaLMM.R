#!/usr/bin/env Rscript

library('MegaLMM')

time="WD_0712"
n=1000
run=3

run_id=sprintf('MegaLMM/MegaLMM_%s_test_%s_%s',time,n,run)

MegaLMM_state=readRDS(sprintf('MegaLMM/%s/MegaLMM_state_base.rds',run_id))
MegaLMM_state$current_state=readRDS(sprintf('MegaLMM/%s/current_state.rds',run_id))
MegaLMM_state$Posterior=reload_Posterior(MegaLMM_state,'F')

pdf('test.pdf')
Image(MegaLMM_state$current_state$F)
Image(get_posterior_mean(MegaLMM_state,F))
dev.off()
