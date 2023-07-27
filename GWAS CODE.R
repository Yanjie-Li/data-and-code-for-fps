setwd("G:/qtl-nir/WHEAT4 DATA SET")  
# aver_pheno_snpspec <- readr:: read_rds('aver_pheno_snpspec.rds')
# aver_pheno_snpspec[is.na(aver_pheno_snpspec)] = 0
aver_pheno_snpspec <- readr:: read_rds('all.rds')
aver_pheno_snpspec <-  as.data.frame(aver_pheno_snpspec)


aver_pheno_snpspec <- aver_pheno_snpspec%>% mutate(Fam = factor(Fam),
                                                   month = factor(month))%>%dplyr:: select(-c(chl,N,P,S))
aver_pheno_snpspec <- aver_pheno_snpspec%>% mutate(sample = Fam)%>%
  dplyr::  select(sample,everything())%>%
  separate(Fam,c("site","block",'fam'),sep = "([_])") %>% dplyr:: select(month,pre.chl,pre.n,everything())
# aver_pheno_snpspec <- aver_pheno_snpspec%>% drop_na()
 
aver_pheno_snpspec <- aver_pheno_snpspec%>%drop_na(site)
aver_pheno_snpspec <-aver_pheno_snpspec%>%filter(complete.cases(.))
nir2 <- sapply(aver_pheno_snpspec[,-c(1:12)], function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))

aver_pheno_snpspec <-cbind.data.frame(aver_pheno_snpspec[, c(1:12)], nir2) %>%
  mutate(block = as.numeric(block),
         site = as.numeric(site),
         fam = as.numeric(fam))%>%
  group_by(month,sample) %>%
  summarise_at(vars(-treeID), funs(mean(., na.rm=TRUE)))

write_rds(aver_pheno_snpspec,'aver_pheno_snpspec_xy.rds')
aver_pheno_snpspec <- readr:: read_rds('G:/qtl-nir/WHEAT4 DATA SET/aver_pheno_snpspec_xy.rds')




library(statgenGWAS)
Pheno <- aver_pheno_snpspec[,c(1:11)] %>% dplyr::select(sample,month,site, 
                                                 block, fam,max_height, crownArea,pre.chl,pre.n)%>%
  rename(Chl=pre.chl,
         N =pre.n
  )
Pheno <- as.data.frame(Pheno)
nir <- as.data.frame(aver_pheno_snpspec[,-c(3:4,5:11)]) %>% dplyr::select(sample,everything())

summary(factor(aver_pheno_snpspec$sample))
library(data.table) # v1.9.5+
library(tidyr)
Markers <- pivot_wider(nir, names_from = "month", values_from = c(3:28), names_sep = ".")


# Markers <- reshape(nir, idvar = "sample", timevar = "month", direction = "wide")
names(Markers)[1] <-  "genotype"
Markers <- as.data.frame(Markers)
# Markers<- Markers%>% drop_na()
rownames(Markers) <-Markers$genotype
map <-  data.frame(t(Markers))

map$che <-sapply(str_split(rownames(map), '[.]'), function(x){
  y=x[2]
})
factor(map$che)
 
map$che <-  factor(map$che,c( 'Jan','Mar','Apr', 'May', 'June', 'July','Aug', 'Sep','Oct','Nov','Dec'))
str(Markers)
map$snpnames <- rownames(map)
map <- map%>% dplyr::select(snpnames,che)
map <-map[-1,]

map$position <-  sapply(str_split(rownames(map), '[.]'), function(x){
  y=x[1]
})

names(map)[2:3] <- c("chr", "pos")
map$pos <- as.numeric(as.factor(map$pos))
map$chr <- as.numeric(as.factor(map$chr))
gDataDrops <- createGData(geno = Markers, map = map)
factor(Pheno$site)

colnames(Pheno)[colnames(Pheno) == "sample"] <- "genotype"
# Pheno$exper <- '1'
## Select relevant columns and convert data to a list.
dropsPhenoList <- split(x = Pheno[c("genotype", "max_height",'crownArea',
                                    'Chl','N')], 
f = Pheno[["site"]])

gDataDrops <- createGData(gData = gDataDrops, pheno = dropsPhenoList)



# summary(gDataDrops, trials = "1")
gDataDropsDedup <- codeMarkers(gDataDrops, impute =T, verbose = TRUE,removeDuplicates = T) 

gDataDropsDedup
# ## Copy gData object.
# gDataDropsMiss <- gDataDrops
# ## Add random missing values to 1% of the values in the marker matrix.
# set.seed(1)
# nVal <- nrow(gDataDropsMiss$markers) * ncol(gDataDropsMiss$markers)
# gDataDropsMiss$markers[sample(x = 1:nVal, size = nVal / 100)] <- NA
# 
# gDataDropsImputed <- codeMarkers(gData = gDataDropsMiss, 
#                                  nMissGeno = 0.01, 
#                                  removeDuplicates = F,
#                                  nMiss = 0.01, 
#                                  impute = TRUE, 
#                                  imputeType = "random", 
#                                  verbose = TRUE)
library(tictoc)
tictoc::tic()
GWASDrops <- runSingleTraitGwas(gData = gDataDropsDedup,
                                GLSMethod = "multi",
                                # trials ='1',
                                traits =c("max_height",'crownArea',
                                          'Chl','N'
                                ))
tictoc::toc()
readr::write_rds(GWASDrops,'2022-02-15_GWASDrops.rds')

print(GWASDrops$GWAResult$`1` , row.names = FALSE)
print(GWASDrops$signSnp$`1`, row.names = FALSE)


rownames(GWASDrops$kinship$`4`)
GWASDrops <- readRDS('g:/qtl-nir/WHEAT4 DATA SET/2022-02-15_GWASDrops.rds')

fgr <-GWASDrops$GWAResult$Sep
plot(GWASDrops , plotType = "manhattan",  trial=1,
     trait = c( "max_height"),yThr = 1.5 )
plot(GWASDrops, plotType = "manhattan",  trial=1,
     trait = c( "crownArea"),yThr = 1.5 )
plot(GWASDrops, plotType = "manhattan",  trial=1,
     trait = c( "Chl"),yThr = 1.5 )
plot(GWASDrops, plotType = "manhattan",  trial=1,
     trait = c( "N"),yThr = 1.5 )
plot(GWASDrops, plotType = "manhattan", 
     trait = c( "P"),yThr = 1.5 )
plot(GWASDrops, plotType = "manhattan", 
     trait = c( "S"),yThr = 1.5 )

kinsp <- GWASDrops$kinship$`10`
heatmap(kinsp)


plot(GWASDrops,  trial=2,plotType = "qtl")
library(qqman)
library(tidyverse)
snpsOfInterest
 
datamatou <- GWASDrops$GWAResult$`2`

datamatou<- as.data.frame(datamatou)
datamatou_re <- reshape(datamatou, idvar = "snp", timevar = "trait", direction = "wide")

datamatou_re <- datamatou_re %>%dplyr::select(snp,chr.max_height,pos.max_height,pValue.max_height ,pValue.crownArea,
                                        pValue.Chl, pValue.N
)%>% as_data_frame()
names(datamatou_re)[1:3] <- c( 'SNP', 'CHR', 'BP')


datamatou_re <- datamatou_re %>% 
  mutate(across(c('pValue.max_height', 'pValue.crownArea',
                  'pValue.Chl', 'pValue.N'),
                ~replace_na(.x, max(rnorm(n = length(datamatou_re[is.na(datamatou_re)] ),
                                          mean = 0.6, sd = 0.1)))))
 


datamatou_re <- datamatou_re %>%
  # mutate(month=as.character(month))%>%
  # Replacing values
  dplyr:: mutate(CHR = replace(CHR, CHR == 1,'Jan'),
                 CHR = replace(CHR, CHR == 2,'Mar'),
                 CHR = replace(CHR, CHR == 3,'Apr'),
                 CHR = replace(CHR, CHR == 4,'May'),
                 CHR = replace(CHR, CHR == 5,'June'),
                 CHR = replace(CHR, CHR == 6,'July'),
                 CHR = replace(CHR, CHR == 7, 'Aug'),
                 CHR = replace(CHR, CHR == 8,'Sep'),
                 CHR = replace(CHR, CHR == 9,'Oct'),
                 CHR = replace(CHR, CHR == 10,'Nov'),
                 CHR = replace(CHR, CHR == 11,'Dec')
  )

datamatou_re <- as.data.frame(datamatou_re)
names(datamatou_re)[4:5] <- c('Height','CA')
# datamatou_re <- datamatou_re%>% drop_na()
SNPs <- list(
  datamatou_re$SNP[datamatou_re$pValue.height<0.001]
  # datamatou_re$SNP[datamatou_re$pValue.ca<0.01],
  # datamatou_re$SNP[datamatou_re$pValue.Chl<0.01],
  # datamatou_re$SNP[datamatou_re$pValue.N<0.01],
  # datamatou_re$SNP[datamatou_re$pValue.P<0.01],
  # datamatou_re$SNP[datamatou_re$pValue.S<0.01]
)
library("CMplot")

CMplot(datamatou_re[,c(1:3,4)],  
       plot.type="m",multracks=TRUE,threshold=c(0.005,0.0001),threshold.lty=c(1,2), band=1,
       threshold.lwd=c(5,5), threshold.col=c("black","blue"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green","blue"),
       signal.cex=2, 
       mar= c(5,15,5,5),
       trait.legend.ncol=2,
       dpi=30,file.output=F,verbose=TRUE,width=22,height=12,file="pdf",memo="",
       highlight.col=c("red"),highlight.cex=1,highlight.pch=16, highlight.text=SNPs,      
       highlight.text.col=c("red"),
       highlight.text.cex=1,
       # cir.legend=F,
       cir.legend.cex=1,
       cex.axis=1,
       highlight=SNPs
)


SNPs <- list(
  datamatou_re$SNP[datamatou_re$Height<0.00001],
  datamatou_re$SNP[datamatou_re$CA<0.00001]
  # datamatou_re$SNP[datamatou_re$pValue.Chl<0.01],
  # datamatou_re$SNP[datamatou_re$pValue.N<0.001] 
  # datamatou_re$SNP[datamatou_re$pValue.P<0.001],
  # datamatou_re$SNP[datamatou_re$pValue.S<0.001]
)

CMplot(datamatou_re[,c(1:3,4:5)],  
       plot.type="m",threshold=c(0.01,0.001),threshold.lty=c(1,2), band=1,
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green","blue"),
       signal.cex=1, 
       mar= c(5,15,5,5),
       trait.legend.ncol=2,
       dpi=300,file.output=T,verbose=TRUE,width=22,height=12,file="pdf",memo="",
       highlight.col=c("red","blue","green"),highlight.cex=1,highlight.pch=16, 
       highlight.text=SNPs,      
       highlight.text.col=c("red","blue","green"),
       highlight.text.cex=1,
       cir.legend.cex=1,
       outward=TRUE,
       cex.axis=1,
       highlight=SNPs
)

c1 <- CMplot(datamatou_re[,c(1:3,4:5)],multracks=F,  plot.type="m", col=c("grey30","grey60"), 
       LOG10=TRUE,  threshold=c(0.0001,0.00001),highlight.text=SNPs,  
       threshold.lty=c(1,2), threshold.lwd=c(1,1), 
       threshold.col=c("black","grey"), highlight=SNPs,
       amplify=TRUE,chr.den.col=NULL,highlight.text.col=c("red","blue"),
       highlight.text.cex=1,
       signal.col=c("red"),signal.cex=c(1,1),
       signal.pch=c(19,19),file="tiff",memo="",
       dpi=300,file.output=F,verbose=TRUE)

heidwd <-datamatou_re[,c(1:3,4 )]
heidwd <-heidwd %>% drop_na()

heidwd$BP <- as.numeric(heidwd$BP)
# heidwd$CHR <- as.numeric(heidwd$CHR)
heidwd$CHR <- factor(heidwd$CHR,
                            unique(heidwd$CHR))

data_cum <- heidwd %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  dplyr::select(CHR, bp_add)

gwas_data <- heidwd %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))


library(ggrepel)
library(ggtext)
p4 <-  ggplot(gwas_data, aes(x = bp_cum, y =  -log10(Height), 
                             color = as_factor(CHR), size =  -log10(Height))) +
  # geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
  # scale_y_continuous(expand = c(0,0), limits = c(-1, 1)) +
  scale_color_manual(values = rep(c('grey','grey2'),11)) +
  scale_size_continuous(range = c(2,6)) +
  labs(x = NULL, 
       y="-log<sub>10</sub>(p)")+
  # y = "-log<sub>10</sub>(p)") + 
  theme_minimal(base_size = 20) +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(1,1,1,1,'cm'),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 0, size = 15, vjust = 0.5)
  )+  ggnewscale::new_scale_color()+
  geom_text_repel(data = subset(gwas_data, -log10(Height) > 5), col='red',
                  aes(label = SNP,size=4), 
                  box.padding = unit(0.45, "lines"))+
  geom_point(data = subset(gwas_data, -log10(Height) > 5), aes(col='red',size=7) ) +
  geom_point(data = subset(gwas_data, -log10(Height) > 4 &-log10(Height) < 4.9),
  aes(col='pink',size=3) ) +
  annotate("rect", xmin = c(1,26,52,78,104,130,156,184,210,236,262),
           
           xmax = c(26,52,78,104,130,156,184,210,236,262,288), ymin = -Inf, ymax =Inf,
           fill =  rep(c('grey2','white'),6)[1:11],
           alpha = .2)+geom_hline(yintercept=4, linetype="dashed", 
                                  color = "blue", size=1)+
  geom_hline(yintercept=5, linetype="dashed", 
     color = "grey2", size=1)

 

p4 
heidwd <-datamatou_re[,c(1:3,5 )]
heidwd <-heidwd %>% drop_na()

heidwd$BP <- as.numeric(heidwd$BP)
# heidwd$CHR <- as.numeric(heidwd$CHR)
heidwd$CHR <- factor(heidwd$CHR,
                     unique(heidwd$CHR))

data_cum <- heidwd %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  dplyr::select(CHR, bp_add)

gwas_data <- heidwd %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))


library(ggrepel)
p5 <-  ggplot(gwas_data, aes(x = bp_cum, y =  -log10(CA), 
                             color = as_factor(CHR), size =  -log10(CA))) +
  # geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
  # scale_y_continuous(expand = c(0,0), limits = c(-1, 1)) +
  scale_color_manual(values = rep(c('grey','grey2'),11)) +
  scale_size_continuous(range = c(2,6)) +
  labs(x = NULL, y=NULL)+
       # y="-log<sub>10</sub>(p)")+
  # y = "-log<sub>10</sub>(p)") + 
  theme_minimal(base_size = 20) +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(1,1,1,1,'cm'),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 0, size = 15, vjust = 0.5)
  )+  ggnewscale::new_scale_color()+
  geom_text_repel(data = subset(gwas_data, -log10(CA) > 5), col='red',
                  aes(label = SNP,size=4), 
                  box.padding = unit(0.45, "lines"))+
  geom_point(data = subset(gwas_data, -log10(CA) > 5), aes(col='red',size=7) ) +
  geom_point(data = subset(gwas_data, -log10(CA) > 4 &-log10(CA) < 4.9),
             aes(col='pink',size=3) ) +
  annotate("rect", xmin = c(1,26,52,78,104,130,156,184,210,236,262),
           
           xmax = c(26,52,78,104,130,156,184,210,236,262,288), ymin = -Inf, ymax =Inf,
           fill =  rep(c('grey2','white'),6)[1:11],
           alpha = .2)+geom_hline(yintercept=4, linetype="dashed", 
                                  color = "blue", size=1)+
  geom_hline(yintercept=5, linetype="dashed", 
             color = "grey2", size=1)
p5
library(patchwork) 
p4+ggtitle(label = 'Height')|p5+ggtitle(label = 'CA')   



CMplot(datamatou_re[,c(1:3,4:7)],plot.type="q",box=FALSE,file="jpg",memo="",dpi=300,
       conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
       file.output=F,verbose=TRUE,width=5,height=5)


#Note: if you are not supposed to change the color of signal, 
#          please set signal.col=NULL and highlight.col=NULL.

CMplot(datamatou_re[,c(1:3,4,5)], type="p", LOG10=FALSE,outward=TRUE,
       col=matrix(c("#4DAF4A",NA,NA,"dodgerblue4",
                    
                    "deepskyblue",NA,"dodgerblue1", "olivedrab3", "darkgoldenrod1"), nrow=3, byrow=TRUE),
       plot.type="c",multracks=TRUE,
       threshold.lty=c(1,2), band=0.1,
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green","blue"),
       signal.cex=1, 
       mar= c(5,15,5,5), 
       threshold=c(0.05,0.01),r=1.2,cir.chr.h=1.5,
       
       cir.band=1,
       dpi=30,file.output=T,verbose=TRUE,width=22,height=12,file="pdf",memo="",
       highlight.col=c("red","blue","green"),highlight.cex=1,highlight.pch=16, highlight.text=SNPs,      
       highlight.text.col=c("red","blue","green"),
       highlight.text.cex=1,
       cir.legend.cex=1,
       cex.axis=1,
       highlight=SNPs
)



trait_he <- datamatou%>% filter(trait == 'S')
trait_he <- trait_he%>% select(snp,chr,pos,pValue)
names(trait_he) <- c( 'SNP', 'CHR', 'BP', 'P')
manhattan(trait_he, chr="CHR", bp="BP", snp="SNP", p="P" )
manhattan(trait_he, chr="CHR", bp="BP", snp="SNP", p="P", annotatePval = 0.05)

qq(trait_he$P)
str(trait_he)
trait_he <- as.data.frame(trait_he)
library("CMplot")  
SNPs <- trait_he[trait_he[,4] <0.05,1]

#install.packages("CMplot")
CMplot(trait_he, plot.type="c", r=1.6, cir.legend=TRUE,file.output=F,highlight.text=genes,
       outward=TRUE, cir.legend.col="black", cir.chr.h=.1 ,chr.den.col="orange", file="jpg",
       memo="", dpi=300)

CMplot(trait_he[,c(1:3,4)], plot.type="m",LOG10=TRUE,highlight=SNPs,
       
       threshold.lwd=c(1,1), threshold.col=c("black","grey"),bin.size=1e6,
       chr.den.col=NULL,
       highlight.col=c("red","blue","green"),highlight.cex=1,highlight.pch=c(15:17), highlight.text=SNPs,      
       highlight.text.col=c("red","blue","green"),
       
       threshold=0.05,threshold.lty=2,
       signal.col=c("red","green","blue"), 
       amplify=T,file="jpg",memo="",dpi=30,file.output=F,verbose=TRUE,width=14,height=6)
# Note:
# 'highlight', 'highlight.text', 'highlight.text.xadj', 'highlight.text.yadj' could be vector or list, if it is a vector, 






datamatou <- GWASDrops$GWAResult$Nov

datamatou<- as.data.frame(datamatou)
datamatou_re <- reshape(datamatou, idvar = "snp", timevar = "trait", direction = "wide")

datamatou_re <- datamatou_re %>% select(snp,chr.height,pos.height,LOD.height,LOD.ca,LOD.Chl)%>% as_data_frame()
names(datamatou_re)[1:3] <- c( 'SNP', 'CHR', 'BP')
datamatou_re[is.na(datamatou_re)] = 0
datamatou_1 <- datamatou_re%>% select(CHR,BP,LOD.height)
plot(datamatou_1$CHR,datamatou_1$LOD.height)
out.hk <- scanone(datamatou_1, method="hk")
