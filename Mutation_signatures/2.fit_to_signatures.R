
library(MutationalPatterns)
signatures = read.csv('/media/desk16/szz077/team/metagene/post_treatment/code/96muntation/result/matrix/W_new.csv')
mut_mat=read.csv('/media/desk16/szz077/team/metagene/post_treatment/code/96muntation/result/matrix/sample_H.csv')

signatures=as.matrix(signatures)
mut_mat=as.matrix(mut_mat)
fit_res <- fit_to_signatures(mut_mat, signatures)

p <- plot_contribution(fit_res$contribution, coord_flip = FALSE, mode = "absolute")
p
data=fit_res$contribution
write.csv(data,file='/media/desk16/szz077/team/metagene/post_treatment/code/96muntation/cosine_similarity/contribution.csv')
