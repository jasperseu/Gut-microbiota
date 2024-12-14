library(vegan)
library(ggplot2)
library(ggrepel)

#MEL_SNV
group<-read.csv('/media/desk16/szz077/team/metagene/metadate_new.csv',sep = ',',row.names = 'run')
#path = '/media/desk16/szz077/team/metagene/result/PERMANOVA/MEL/clear_data_no_empty_row'
path='/media/desk16/szz077/team/metagene/result/PERMANOVA_INDEL/MEL/clear_data_no_empty_row'
files=list.files(path,full.names = TRUE)
for (file in files){
  print(file)
  df=read.csv(file,sep=',',row.names = 'run')
  df=t(df)
  distance<-vegdist(df,method='bray')
  group1=group[match(rownames(df),rownames(group)),]
  adonis_result_study <- adonis2(distance~response, group1, permutations = 999)
  print(adonis_result_study)
}

#MEL_Bacterial_abundance
df<- read.csv('data/combine_abd_MEL.csv',sep = ',')
group<-read.csv('metadate_new.csv',sep = ',')
group=group[which(group$cancer=="STAD"),]
df1<-t(df)
distance<-vegdist(df1,method='bray')
pcoa<- cmdscale(distance,k=(nrow(df1)-1),eig=TRUE)
plot_data<-data.frame({pcoa$point})[1:2]

names(plot_data)[1:2]<-c('PCoA1','PCoA2') 
eig=pcoa$eig
group1<-group['response']
data<-plot_data[match(rownames(group),rownames(plot_data)),]
data<-data.frame(group,plot_data)
head(data)
tail(data)
adonis_result_study <- adonis2(distance~response, data, permutations = 999)
print(adonis_result_study)



#NSCLC_SNV
group<-read.csv('/media/desk16/szz077/team/metagene/metadate_new.csv',sep = ',',row.names = 'run')
path='/media/desk16/szz077/team/metagene/result/PERMANOVA_INDEL/NSCLC/clear_data_no_empty_row'
files=list.files(path,full.names = TRUE)
for (file in files){
  print(file)
  df=read.csv(file,sep=',',row.names = 'run')
  df=t(df)
  distance<-vegdist(df,method='bray')
  group1=group[match(rownames(df),rownames(group)),]
  adonis_result_study <- adonis2(distance~response, group1, permutations = 999)
  print(adonis_result_study)
}

#NSCLC_Bacterial_abundance
df<- read.csv('data/combine_abd_NSCLC.csv',sep = ',')
group<-read.csv('metadate_new.csv',sep = ',')
group=group[which(group$cancer=="STAD"),]
df1<-t(df)
distance<-vegdist(df1,method='bray')
pcoa<- cmdscale(distance,k=(nrow(df1)-1),eig=TRUE)
plot_data<-data.frame({pcoa$point})[1:2]

names(plot_data)[1:2]<-c('PCoA1','PCoA2') 
eig=pcoa$eig
group1<-group['response']
data<-plot_data[match(rownames(group),rownames(plot_data)),]
data<-data.frame(group,plot_data)
head(data)
tail(data)
adonis_result_study <- adonis2(distance~response, data, permutations = 999)
print(adonis_result_study)


#RCC_SNV
group<-read.csv('/media/desk16/szz077/team/metagene/metadate_new.csv',sep = ',',row.names = 'run')
path='/media/desk16/szz077/team/metagene/result/PERMANOVA_INDEL/RCC/clear_data_no_empty_row'
files=list.files(path,full.names = TRUE)
for (file in files){
  print(file)
  df=read.csv(file,sep=',',row.names = 'run')
  df=t(df)
  distance<-vegdist(df,method='bray')
  group1=group[match(rownames(df),rownames(group)),]
  adonis_result_study <- adonis2(distance~response, group1, permutations = 999)
  print(adonis_result_study)
}

#RCC_Bacterial_abundance
df<- read.csv('data/combine_abd_RCC.csv',sep = ',')
group<-read.csv('metadate_new.csv',sep = ',')
group=group[which(group$cancer=="STAD"),]
df1<-t(df)
distance<-vegdist(df1,method='bray')
pcoa<- cmdscale(distance,k=(nrow(df1)-1),eig=TRUE)
plot_data<-data.frame({pcoa$point})[1:2]

names(plot_data)[1:2]<-c('PCoA1','PCoA2') 
eig=pcoa$eig
group1<-group['response']
data<-plot_data[match(rownames(group),rownames(plot_data)),]
data<-data.frame(group,plot_data)
head(data)
tail(data)
adonis_result_study <- adonis2(distance~response, data, permutations = 999)
print(adonis_result_study)


#GIC_SNV
group<-read.csv('/media/desk16/szz077/team/metagene/metadate_new.csv',sep = ',',row.names = 'run')
path='/media/desk16/szz077/team/metagene/result/PERMANOVA_INDEL/GIC/clear_data_no_empty_row'
files=list.files(path,full.names = TRUE)
for (file in files){
  print(file)
  df=read.csv(file,sep=',',row.names = 'run')
  df=t(df)
  distance<-vegdist(df,method='bray')
  group1=group[match(rownames(df),rownames(group)),]
  adonis_result_study <- adonis2(distance~response, group1, permutations = 999)
  print(adonis_result_study)
}

#GIC_Bacterial_abundance
df<- read.csv('data/combine_abd_GIC.csv',sep = ',')
group<-read.csv('metadate_new.csv',sep = ',')
group=group[which(group$cancer=="STAD"),]
df1<-t(df)
distance<-vegdist(df1,method='bray')
pcoa<- cmdscale(distance,k=(nrow(df1)-1),eig=TRUE)
plot_data<-data.frame({pcoa$point})[1:2]

names(plot_data)[1:2]<-c('PCoA1','PCoA2') 
eig=pcoa$eig
group1<-group['response']
data<-plot_data[match(rownames(group),rownames(plot_data)),]
data<-data.frame(group,plot_data)
head(data)
tail(data)
adonis_result_study <- adonis2(distance~response, data, permutations = 999)
print(adonis_result_study)