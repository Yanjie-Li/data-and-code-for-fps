month <- c('Jan','Mar','Apr', 'May', 'June', 'July','Aug', 'Sep','Oct','Nov','Dec')
setwd("G:/qtl-nir/WHEAT4 DATA SET")  
fdfdr <- readr::read_rds('2022-02-15_GWASDrops.rds')
kinsp <- fdfdr$kinship 
names(kinsp) <- c('Jan','Mar','Apr', 'May', 'June', 'July','Aug', 'Sep','Oct','Nov','Dec')

std_mean <- function(x) sd(x)/sqrt(length(x))

aver_pheno_snpspec <- readr:: read_rds('G:/qtl-nir/WHEAT4 DATA SET/aver_pheno_snpspec_xy.rds')
aver_pheno_snpspec <- as.data.frame(aver_pheno_snpspec)
aver_pheno_snpspec$month

cor_PS <- function(month,kinsp,time=1:500){
  e <- new.env()
  e$h_cr_nir_GB <- list()
  e$h_cr_nir_GK <- list()
  e$h_cr_nir_GB_p <- list() 
  e$h_cr_nir_GK_p <- list() 
  
  e$ca_cr_nir_GB <- list()
  e$ca_cr_nir_GB_p <- list()
  e$dt <- list()
  e$cr_all <- list()
  
  lapply(time, function(z){ 
    lapply(month, function(jj){
      expr <- tryCatch({
        library(BGGE)
        library(prospectr)
        library(BGLR)
        library(tidyverse)
        library(tidyverse)
        library(sf)
        library(data.table)
        message( paste0(jj))
        dec_data <- aver_pheno_snpspec %>% filter(month == jj[1]) %>% drop_na(23:37)
        nir <-  dec_data[,c(2,12:37)]
        names(nir)
        # nir <-  as.matrix(dec_data[,c(12:37)])
        # nir <- nir %>%
        #   group_by(sample) %>%
        #   dplyr::summarise(
        #     result_Red = mean(Red,na.rm=T),
        #     result_Blue = mean(Blue,na.rm=T),
        #     result_Green = mean(Green,na.rm=T),
        #     result_RedEdge = mean(Rededge,na.rm=T),
        #     ndvi = mean(NDVI,na.rm=T),
        #     osavi = mean(OSAVI,na.rm=T),
        #     gndvi = mean(GNDVI,na.rm=T),
        #     savi = mean(SAVI,na.rm=T),
        #     msavi = mean(MSAVI,na.rm=T),
        #     gci = mean(GCI,na.rm=T),
        #     RECI = mean(RECI,na.rm=T),
        #     LCI = mean(LCI,na.rm=T),
        #     RGBVI = mean(RGBVI,na.rm=T),
        #     NDRE = mean(NDRE,na.rm=T),
        #     MACI = mean(MACI,na.rm=T),
        #     ARI = mean(ARI,na.rm=T),
        #     MARI = mean(MARI,na.rm=T),
        #     RVI = mean(RVI,na.rm=T),
        #     EVI = mean(EVI,na.rm=T),
        #     TVI = mean(TVI,na.rm=T),
        #     CVI = mean(CVI,na.rm=T),
        #     CI_G = mean(CI_G,na.rm=T),
        #     GLI = mean(GLI,na.rm=T),
        #     CI_RE = mean(CI_RE,na.rm=T),
        #     TGI = mean(TGI,na.rm=T),
        #     result_NIR = mean(NIR,na.rm=T)
        #   )
        rownames(nir) <- nir$sample
        d <- rownames(nir)
        nir <- nir[,-1]
        nir <- as.matrix(nir)
        rownames(nir) <- d
        GB1<-tcrossprod(nir)/ncol(nir)
        dist1<-as.matrix(dist(nir)^2)
        GK1 <-exp(-dist1/median(dist1))   # Gaussian Kernel  
        fd <- kinsp[names(kinsp) == jj[1]]
        # message(paste0(names(fd),'_',rownames(nir)))
        lapply(fd, function(x1){
          # x2 <-list(unlist(x1))
          # lapply(x2, function(ls){
          #   lapply(ls, function(ls2){
          tryCatch({
            library(tictoc)
            tic("start")
            un1 <- sort(union(rownames(x1), colnames(x1)))
            m1 <- matrix(0, dimnames = list(un1, un1), ncol=length(un1), nrow=length(un1))
            i1 <- match(rownames(GB1), rownames(x1), nomatch = 0)
            j1 <- match(colnames(GB1), colnames(x1), nomatch = 0)
            m1 <- x1[i1,j1]
            m1 <- as.matrix(m1)
          
            #  NIR1
            K0<-list(list(Kernel=GB1,Type="D")) #NIR1
            K_gk <-list(list(Kernel=GK1,Type="D")) #NIR1
            y <- dec_data[,10]#height
            names(y) <- dec_data$sample
            yna<-y
            # sel <- kenStone(GB1, k = 100, pc = .99)
            # set.seed(122);
            e$dt[[z]] <- sort(sample(length(y),length(y)*.8 ))
            test <- y[-(e$dt[[z]])]
            yna[names(test)] <-NA
            fit <- BGGE(y = yna, K = K0, ne =1, ite=3000,burn=200,thin = 1)
            e$h_cr_nir_GB[[z]]  <- round(cor(y[names(test)],fit$yHat[-(e$dt[[z]])] ),2)
            e$h_cr_nir_GB[[z]]   <- as.data.frame(e$h_cr_nir_GB[[z]]  )
            rownames(e$h_cr_nir_GB[[z]] ) <- colnames(e$h_cr_nir_GB[[z]])
            colnames(e$h_cr_nir_GB[[z]] ) <- 'value'
            ##GK
            fit_gk <- BGGE(y = yna, K = K_gk, ne =1, ite=3000,burn=200,thin = 1)
            e$h_cr_nir_GK[[z]]  <- round(cor(y[names(test)],fit_gk$yHat[-(e$dt[[z]])] ),2)
            e$h_cr_nir_GK[[z]]   <- as.data.frame(e$h_cr_nir_GK[[z]]  )
            rownames(e$h_cr_nir_GK[[z]] ) <- colnames(e$h_cr_nir_GK[[z]])
            colnames(e$h_cr_nir_GK[[z]] ) <- 'value'
            
            
            
            # #   pedigree +NIR1
            k1 <- list(list(Kernel=m1,Type="D"),list(Kernel=GB1,Type="D")) #pedigree +NIR1
            fit1 <- BGGE(y = yna, K = k1, ne =1, ite=3000,burn=200,thin = 1)
            e$h_cr_nir_GB_p[[z]] <-round(cor(y[names(test)],fit1$yHat[-(e$dt[[z]])] ),2)
            e$h_cr_nir_GB_p[[z]]  <- as.data.frame(e$h_cr_nir_GB_p[[z]])
            rownames(e$h_cr_nir_GB_p[[z]]) <- colnames(e$h_cr_nir_GB_p[[z]])
            colnames(e$h_cr_nir_GB_p[[z]]) <- 'value'
            ##GK
            
            k1_GK <- list(list(Kernel=m1,Type="D"),list(Kernel=GK1,Type="D")) #pedigree +NIR1
            fit1_GK <- BGGE(y = yna, K = k1_GK, ne =1, ite=3000,burn=200,thin = 1)
            e$h_cr_nir_GK_p[[z]] <-round(cor(y[names(test)],fit1_GK$yHat[-(e$dt[[z]])] ),2)
            e$h_cr_nir_GK_p[[z]]  <- as.data.frame(e$h_cr_nir_GK_p[[z]])
            rownames(e$h_cr_nir_GK_p[[z]]) <- colnames(e$h_cr_nir_GK_p[[z]])
            colnames(e$h_cr_nir_GK_p[[z]]) <- 'value'
            
            #   NIR1
            
            y2 <- dec_data$crownArea#ca
            names(y2) <-dec_data$sample
            yna2 <-y2
            # sel <- kenStone(GB1, k = 100, pc = .99)
            # set.seed(122);dt2 <- sort(sample(length(y2),length(y2)*.8 ))
            # test2<-y2[-dt2]
            yna2[names(test)]<-NA
            fit2 <- BGGE(y = yna2, K = K0, ne =1, ite=3000,burn=200,thin = 1)
            e$ca_cr_nir_GB[[z]]  <- round( cor(y2[names(test)],fit2$yHat[-(e$dt[[z]])] ),2)
            e$ca_cr_nir_GB[[z]]  <- as.data.frame(e$ca_cr_nir_GB[[z]] )
            rownames(e$ca_cr_nir_GB[[z]]) <- colnames(e$ca_cr_nir_GB[[z]])
            colnames(e$ca_cr_nir_GB[[z]]) <- 'value'
            #GK
            
            fit2_GK <- BGGE(y = yna2, K = K_gk, ne =1, ite=3000,burn=200,thin = 1)
            e$ca_cr_nir_GK[[z]]  <- round( cor(y2[names(test)],fit2_GK$yHat[-(e$dt[[z]])] ),2)
            e$ca_cr_nir_GK[[z]]  <- as.data.frame(e$ca_cr_nir_GK[[z]] )
            rownames(e$ca_cr_nir_GK[[z]]) <- colnames(e$ca_cr_nir_GK[[z]])
            colnames(e$ca_cr_nir_GK[[z]]) <- 'value'
            
            
            fit3 <- BGGE(y = yna2, K = k1, ne =1, ite=3000,burn=200,thin = 2)
            e$ca_cr_nir_GB_p[[z]] <- round( cor(y2[names(test)],fit3$yHat[-(e$dt[[z]])]),2)
            e$ca_cr_nir_GB_p[[z]]  <- as.data.frame(e$ca_cr_nir_GB_p[[z]] )
            rownames(e$ca_cr_nir_GB_p[[z]]) <- colnames(e$ca_cr_nir_GB_p[[z]])
            colnames(e$ca_cr_nir_GB_p[[z]]) <- 'value'
            
            fit3_GK <- BGGE(y = yna2, K = k1_GK, ne =1, ite=3000,burn=200,thin = 2)
            e$ca_cr_nir_GK_p[[z]] <- round( cor(y2[names(test)],fit3_GK$yHat[-(e$dt[[z]])]),2)
            e$ca_cr_nir_GK_p[[z]]  <- as.data.frame(e$ca_cr_nir_GK_p[[z]] )
            rownames(e$ca_cr_nir_GK_p[[z]]) <- colnames(e$ca_cr_nir_GK_p[[z]])
            colnames(e$ca_cr_nir_GK_p[[z]]) <- 'value'
            # c
            cr_all <- rbind(e$h_cr_nir_GB[[z]],
                            e$h_cr_nir_GB_p[[z]],
                            e$ca_cr_nir_GB[[z]] ,
                            e$ca_cr_nir_GB_p[[z]],
                            e$h_cr_nir_GK[[z]],
                            e$h_cr_nir_GK_p[[z]],
                            e$ca_cr_nir_GK_p[[z]],
                            e$ca_cr_nir_GK[[z]]
                            
            )
            cr_all <- data.frame(cr_all)
            cr_all$month <- jj[1]
            message( paste0(names(fd),'_','loop_',z))
            # cr_all1 <- cr_all     %>% invoke(rbind,.)
            ###calibration data 
            # cr_all1 <- as.data.frame(do.call(rbind,  cr_all, quote = T))#sperate r
            tictoc::toc()
           return(cr_all)
            
          },
          error = function(e) {
            message('Caught an error!')
            cat("ERROR :", conditionMessage(e), "\n")
            print(e)})
          
        })
        
      },error = function(e) {
        message('Caught an error!')
        cat("ERROR :", conditionMessage(e), "\n")
        print(e)}, 
      print(paste("processing Loop",jj, sep = "_"))
      )
    })
    
  })
}

cor_allmonth <- cor_PS(month,kinsp = kinsp,time =1:100)
write_rds(cor_allmonth,'cor_allmonth.rds')
cor_allmonth <- readRDS ('G:/qtl-nir/WHEAT4 DATA SET/cor_allmonth.rds')
library(dplyr)
library(data.table)
library(tidyverse)
cor_data <- cor_allmonth %>%  unlist(recursive = F)%>%  unlist(recursive = F) %>% invoke(rbind,.)
cor_data$type <- rownames(cor_data)
cor_data <- cor_data%>% tidyr:: separate(type, c('month1','trait','cor','datatype','model', 'p' ),'[_]|[.]' )
cor_data<- cor_data %>%  
  mutate(trait = replace(trait, trait == 'e$h', 'Height')) %>%
  mutate(trait = replace(trait, trait == 'e$ca', 'CA')) %>%
  mutate(datatype = replace(datatype, datatype %like% '[[z]]', 'nir')) %>%
  mutate(p = replace(p, p %like% '[[z]]', 'p')) %>% select(-3)

cor_data<- cor_data %>%  
  separate(model, c('model','dd' ),'[[]|[0-9]' ) %>% select(-dd)


cor_data <- cor_data%>% mutate(typed=paste0(datatype,'_', p) ) 
cor_dataDD <- cor_data %>% 
  group_by(month,trait,typed,model) %>%
  summarise_at(vars(starts_with("value")),
               funs(mean = mean(.),
                    min= min(.),
                    max= max(.),
                    sd = sd(.)))



cor_dataDD$month <-  factor(cor_dataDD$month,c( 'Jan','Mar','Apr', 'May', 'June', 'July','Aug', 'Sep','Oct','Nov','Dec'))



dat <- data.table(cor_dataDD)
dat[,y_min := mean*0.5, by = cor_dataDD]
dat[,y_max:= mean*1.5, by = cor_dataDD]

cor_dataDD <- cor_dataDD %>%
  arrange(case_when(
    month== "Jan" ~ -mean,
    month== "Mar" ~ -mean, 
    month== "Apr" ~ -mean ,
    month== "May" ~ -mean, 
    month== "June" ~ -mean ,
    month== "July" ~ -mean ,
    month== "Aug" ~ -mean, 
    month== "Sep" ~ -mean, 
    month== "Oct" ~ -mean, 
    month== "Nov" ~ -mean, 
    month== "Dec" ~ -mean 
    
    # endsWith(month, "Mar") ~ "Spr",
    # endsWith(month, "Apr") ~ "Spr",
    # endsWith(month, "May") ~ "Spr",
    # endsWith(month, "Jun") ~ "Sum",
    # endsWith(month, "July") ~ "Sum",
    # endsWith(month, "Aug") ~ "Sum",
    # endsWith(month, "Sep") ~ "Aut",
    # endsWith(month, "Oct") ~ "Aut",
    # endsWith(month, "Nov") ~ "Aut",
    # endsWith(month, "Dec") ~ "Win" 
  ))
library(tidytext)
library(ggrepel)

cor_dataDD <- cor_dataDD %>%   mutate(typed = as.factor(typed),
                                      names = reorder_within(month, mean, typed)) 


ggplot(cor_dataDD, aes( x = month,
                        y = mean,
                        ymin = mean - sd,
                        ymax = mean + sd,
                        fill = typed
)) +
  facet_grid(trait ~ model, scales = 'free') +
  geom_bar(position = position_dodge(1), stat = "identity") +
  geom_errorbar(
    size = 0.1,
    width = 0.2,
    position = position_dodge(1),
    colour = "black"
  ) + 
  
  geom_point(position = position_dodge(1), aes(y = mean, colour = 'red')) +
  # geom_text(aes(x=month,y=round(value,2),angle=90,
  #               label=round(value,2),hjust=0,vjust=0.2),
  #           size=6,
  #           position = position_dodge(1))+
  theme_classic(base_size = 20) +
  labs(y = 'Prediction ability', x = 'Month') +
  theme(
    axis.text = element_text(colour = 'black'),
    # plot.margin = margin(2,1,1,1,'cm'),
    legend.text = element_text(size = 18),
    legend.margin = margin(0, 1, 1, 1, 'cm'),
    legend.position = 'bottom'
  ) +
  # scale_fill_viridis_d()
  scale_fill_manual(
    name = "",
    values =viridis::viridis(3),
    labels = expression("MV",
                        'MV+P')
  ) +
  geom_blank(data = dat , aes(y = y_min)) +
  geom_blank(data = dat , aes(y = y_max)) +
  guides(colour = 'none')
# +geom_text(aes(label= round(mean,2)))





ggplot(cor_data, aes( x = model ,
                      y = value,
                      # ymin = mean - sd,
                      # ymax = mean + sd,
                      fill = typed
)) +
  facet_grid(trait ~ model, scales = 'free') +
  geom_boxplot()
# geom_bar(position = position_dodge(1), stat = "identity") +
geom_errorbar(
  size = 0.1,
  width = 0.2,
  position = position_dodge(1),
  colour = "black"
) + 
  
  geom_point(position = position_dodge(1), aes(y = mean, colour = 'red')) +
  # geom_text(aes(x=month,y=round(value,2),angle=90,
  #               label=round(value,2),hjust=0,vjust=0.2),
  #           size=6,
  #           position = position_dodge(1))+
  theme_classic(base_size = 28) +
  labs(y = 'Prediction ability', x = 'Month') +
  theme(
    axis.text = element_text(colour = 'black'),
    # plot.margin = margin(2,1,1,1,'cm'),
    legend.text = element_text(size = 18),
    legend.margin = margin(0, 1, 1, 1, 'cm'),
    legend.position = 'bottom'
  ) +
  # scale_fill_viridis_d()
  scale_fill_manual(
    name = "",
    values =viridis::viridis(3),
    labels = expression("MS",
                        'MS+P')
  ) +
  geom_blank(data = dat , aes(y = y_min)) +
  geom_blank(data = dat , aes(y = y_max)) +
  guides(colour = 'none')




