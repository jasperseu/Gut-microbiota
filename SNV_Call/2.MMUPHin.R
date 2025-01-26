library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)

abd=read.csv('abundance_taxa.csv',sep=',',row.names='run')  
meta=read.csv('metadata_new.csv',sep=',',row.names='run')   

meta$study_id=as.factor(meta$study_id)
meta$response=as.factor(meta$response)

fit_adjust_batch <- adjust_batch(feature_abd = abd,
                                 batch = "study_id",   # batch_size的修改条件
                                 data = meta,
                                 covariates='response'
                                 control = list(verbose = FALSE))

#CRC_abd_adj <- fit_adjust_batch$feature_abd_adj
#D_after <- vegdist(t(CRC_abd_adj))
#set.seed(1)
#fit_adonis_after <- adonis2(D_after ~ response, data = meta)
#fit_adonis_after

fit_lm_meta <- lm_meta(feature_abd = CRC_abd_adj,
                       batch = "study_id",
                       exposure = "response",
                       data = meta,
                       covariates=c('age','BMI','sex')
                       control = list(verbose = FALSE))

meta_fits <- fit_lm_meta$meta_fits

meta_fits %>% 
  filter(pval < 0.05) %>% 
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>% 
  ggplot(aes(y = coef, x = feature)) +
  geom_bar(stat = "identity") +
  coord_flip()
