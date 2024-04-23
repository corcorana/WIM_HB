library(tidyverse)
library(mgcv)
library(mgcViz)
library(itsadug)
library(emmeans)
library(ggeffects)
library(RColorBrewer)
library(cowplot)

# global settings
options(contrasts = c("contr.treatment", "contr.poly"))
emm_options(rg.limit = 1e6)

theme_set(theme_classic())
fsize = 8
font_theme <- theme(
  axis.title.x = element_text(size = fsize+2, face="bold"),
  axis.title.y = element_text(size = fsize+2, face="bold"),
  axis.text.x = element_text(size = fsize),
  axis.text.y = element_text(size = fsize),
  strip.text.x = element_text(size = fsize, face="bold"),
  strip.text.y = element_text(size = fsize, face="bold"),
  legend.position = "none",  
  panel.border = element_rect(fill=NA),
  strip.background = element_rect(fill = "gray80"),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA)
)


# import trial data
dat_trl <- read.csv( "WIM_HB_trl_data.txt" )

# replace NaN with NA
dat_trl[dat_trl=="NaN"] <- NA


# recode variables
dat_trl <- dat_trl %>% 
  mutate(subj_id = as.factor(subj_id),
         stim_type = factor(stim_type, levels = c("D", "F"), labels = c("Digit", "Face")),
         trial_type = factor(trial_type, levels = c("G", "N"), labels = c("Go", "NoGo")),
         m_state = factor(ifelse(state>3, 3, state), levels = c("1", "2", "3"), labels = c("ON", "MW", "MB")),
         mstat = ifelse(state>3, 2, state-1),
         pOn = c(0, diff(dat_trl$probe_num) ),
         lagRespType = c(NA, resp_type[1:nrow(.)-1]), 
         lagRespType = ifelse(pOn==T, NA, lagRespType),
         pOff = c(diff(dat_trl$probe_num), 0 ),
         leadRespType = c(resp_type[2:nrow(.)], NA),
         leadRespType = ifelse(pOff==T, NA, leadRespType),
         resp_type = ifelse(resp_type=='C' & lagRespType=='M' & leadRespType=="M", NA, resp_type), 
         resp_type = factor(resp_type, levels = c("H", "M", "F", "C"), labels = c("Hit", "Miss", "FalseAlarm", "CorrectReject"))
  )


# import interpolated IBI and RS data
dat_int <- read.csv( "WIM_HB_IBI_RT_interp.txt" ) %>%
  mutate(subj_id = as.factor(subj_id),
         ln_IBcv_30.s = ( log(IBcv_30)-mean(log(IBcv_30)) ) / sd(log(IBcv_30)),
         ln_RScv_30.s = ( log(RScv_30)-mean(log(RScv_30)) ) / sd(log(RScv_30)),
         ln_IBcv_10.s = ( log(IBcv_10)-mean(log(IBcv_10)) ) / sd(log(IBcv_10)),
         ln_RScv_10.s = ( log(RScv_10)-mean(log(RScv_10)) ) / sd(log(RScv_10))
  )

# replace NaN with NA
dat_int[dat_int=="NaN"] <- NA

# exclude RS estimates derived from insufficient number of trials
dat_int <- filter(dat_trl, time_to_probe > -30) %>%
  group_by(subj_id, probe_num) %>%
  summarise( nRS_30 = sum(!is.na(rspeed)) ) %>%
  ungroup() %>%
  mutate(vRS_30 = nRS_30 >= 15) %>%
  left_join(dat_int, by = c("subj_id", "probe_num") )

dat_int <- filter(dat_trl, time_to_probe > -10) %>%
  group_by(subj_id, probe_num) %>%
  summarise( nRS_10 = sum(!is.na(rspeed)) ) %>%
  ungroup() %>%
  mutate(vRS_10 = nRS_10 >= 5) %>%
  left_join(dat_int, by = c("subj_id", "probe_num") )


# combine probe-locked trial and ibi data
epoch <- c(-30, 0)

prb_dat <- dat_trl %>%
  filter(time_to_probe > epoch[1] & time_to_probe < epoch[2]) %>%
  group_by(subj_id, block_num, probe_num, stim_type, m_state, mstat, vigil) %>%
  summarise(H = sum(resp_type == "Hit", .5, na.rm = T), 
            M = sum(resp_type == "Miss", .5, na.rm = T),
            FA = sum(resp_type == "FalseAlarm", .5, na.rm = T),
            CR = sum(resp_type == "CorrectReject", .5, na.rm = T),
            HR = H/(H+M),
            FAR = FA/(FA+CR)
  ) %>%
  ungroup() %>%
  mutate(dp = qnorm(HR) - qnorm(FAR), cn = rowMeans(cbind(qnorm(HR), qnorm(FAR))), 
         vigil.o = as.ordered(vigil), vigil.s = scale(vigil)
  ) %>%
  left_join(dat_int, by = c("subj_id", "probe_num") )

# unadjusted H / FA rates 
prb_dat <- dat_trl %>%
  filter(time_to_probe > epoch[1] & time_to_probe < epoch[2]) %>%
  group_by(subj_id, probe_num) %>%
  summarise(H_u = sum(resp_type == "Hit", na.rm = T), 
            M_u = sum(resp_type == "Miss", na.rm = T),
            FA_u = sum(resp_type == "FalseAlarm", na.rm = T),
            CR_u = sum(resp_type == "CorrectReject", na.rm = T),
            HR_u = H_u/(H_u+M_u),
            FAR_u = FA_u/(FA_u+CR_u)
  ) %>%
  ungroup() %>%
  mutate(HR_u = ifelse(HR_u == 0, .001, if_else(HR_u == 1, .999, HR_u ) ),
         FAR_u = ifelse(FAR_u == 0, .001, if_else(FAR_u == 1, .999, FAR_u ) )
  ) %>%
  left_join(prb_dat, by = c("subj_id", "probe_num") )

# calculate SDT for short epoch 
epoch <- c(-10, 0)

prb_dat <- dat_trl %>%
  filter(time_to_probe > epoch[1] & time_to_probe < epoch[2]) %>%
  group_by(subj_id, probe_num) %>%
  summarise(H_10 = sum(resp_type == "Hit", .5, na.rm = T), 
            M_10 = sum(resp_type == "Miss", .5, na.rm = T),
            FA_10 = sum(resp_type == "FalseAlarm", .5, na.rm = T),
            CR_10 = sum(resp_type == "CorrectReject", .5, na.rm = T),
            HR_10 = H_10/(H_10+M_10),
            FAR_10 = FA_10/(FA_10+CR_10)
  ) %>%
  ungroup() %>%
  mutate(dp_10 = qnorm(HR_10) - qnorm(FAR_10), cn_10 = rowMeans(cbind(qnorm(HR_10), qnorm(FAR_10)))
  ) %>%
  left_join(prb_dat, by = c("subj_id", "probe_num") )



## pupil data (calculated per trial ~30 prior probe)
dat_pup <- read.csv( "WIM_HB_SW_Pup.txt" )

# replace NaN with NA
dat_pup[dat_pup=="NaN"] <- NA

# exclude lost signal
dat_pup$Pup[dat_pup$Pup<5] <- NA

# factorise subject ID
dat_pup$subj_id <- as.factor(dat_pup$SubID)

# summary stats for scaling on subject level
sumStats <- filter(dat_pup, DistProbe > -30) %>% 
  group_by(subj_id) %>%
  summarise(muPup = mean(Pup, na.rm=T), sdPup = sd(Pup, na.rm=T)) %>%
  ungroup()

dat_pup <- left_join(dat_pup, sumStats, "subj_id") %>%
  mutate(Pup.z = (Pup-muPup)/sdPup
  )

# exclude extreme data and flag epochs missing data from > numNA trials
dat_pup$Pup.z[abs(dat_pup$Pup.z)>4]<-NA

prb_pup <- filter(dat_pup, DistProbe > -30) %>%
  group_by(subj_id, ProbeN) %>%
  summarise(probe_num = mean(ProbeN), nEst = sum(!is.na(Pup.z)), muPup = mean(Pup, na.rm=T), muPup.z = mean(Pup.z, na.rm=T) ) %>%
  ungroup() %>%
  mutate(v30 = nEst>=15) %>%
  left_join(prb_dat, by = c("subj_id", "probe_num") )

prb_pup <- filter(dat_pup, DistProbe > -10) %>%
  group_by(subj_id, ProbeN) %>%
  summarise(nEst10 = sum(!is.na(Pup.z)), muPup10.z = mean(Pup.z, na.rm=T) ) %>%
  ungroup() %>%
  mutate(v10 = nEst10>=5) %>%
  left_join(prb_pup, by = c("subj_id", "ProbeN") )


# descriptives
cntMS <- prb_dat %>% 
  group_by(subj_id, m_state, .drop = FALSE) %>%
  summarise( n = n() ) %>%
  ungroup()


# behav onto m_state
ms.dp <- gam(dp_10 ~ stim_type + m_state + s(subj_id, probe_num, bs="fs", m=1),
             family = gaussian(), method = "REML", data = filter(prb_dat, subj_id!="314") )

ms.cn <- gam(cn_10 ~ stim_type + m_state + s(subj_id, probe_num, bs="fs", m=1),
             family = gaussian(), method = "REML", data = filter(prb_dat, subj_id!="314") )

ms.rs <- gam(RSmu_10 ~ stim_type + m_state + s(subj_id, probe_num, bs="fs", m=1),
             family = gaussian(), method = "REML", data = filter(prb_dat, subj_id!="314" & vRS_10==T) )

ms.rv <- gam(RScv_10 ~ stim_type + m_state + s(subj_id, probe_num, bs="fs", m=1),
             family = Gamma(), method = "REML", data = filter(prb_dat, subj_id!="314" & vRS_10==T) )


# m_state onto time / vigil
tot.ms <- gam(list(mstat 
                   ~ stim_type + scale(probe_num) + s(subj_id, bs="re"),
                   ~ stim_type + scale(probe_num) + s(subj_id, bs="re")
                   ),
              family = multinom(K=2), method = "REML", data = filter(prb_dat, subj_id!="314") )


vig.ms <- gam(list(mstat 
                    ~ stim_type + vigil.s + s(subj_id, probe_num, bs="fs", m=1),
                    ~ stim_type + vigil.s + s(subj_id, probe_num, bs="fs", m=1)
                    ),
               family = multinom(K=2), method = "REML", data = filter(prb_dat, subj_id!="314") )


pup.ms <- gam(list(mstat 
                   ~ stim_type + muPup10.z + s(subj_id, probe_num, bs="fs", m=1),
                   ~ stim_type + muPup10.z + s(subj_id, probe_num, bs="fs", m=1)
                   ),
              family = multinom(K=2), method = "REML", data = filter(prb_pup, subj_id!="314" & v10==T) )


ibi.ms <- gam(list(mstat 
                   ~ stim_type + ln_IBcv_10.s + IBzu_10 + s(subj_id, probe_num, bs="fs", m=1),
                   ~ stim_type + ln_IBcv_10.s + IBzu_10 + s(subj_id, probe_num, bs="fs", m=1)
                   ),
              family = multinom(K=2), method = "REML", data = filter(prb_dat, subj_id!="314") )


# vigil onto IBI
ibi.vo <- gam(vigil ~ stim_type + ln_IBcv_30.s + IBzu_30 + s(subj_id, probe_num, bs="fs", m=1),
              family = ocat(R=4), method = "REML", data = prb_dat )


# pupil onto IBI
ibi.pz <- gam(muPup.z ~ stim_type + ln_IBcv_30.s + IBzu_30 + s(subj_id, probe_num, bs="fs", m=1),
              family = gaussian(), method = "REML", data = filter(prb_pup, v30==T) )


# behav onto IBI
ibi.dp <- gam(dp ~ stim_type + ln_IBcv_30.s + IBzu_30 + s(subj_id, probe_num, bs="fs", m=1),
              family = gaussian(), method = "REML", data = prb_dat )

ibi.cn <- gam(cn ~ stim_type + ln_IBcv_30.s + IBzu_30 + s(subj_id, probe_num, bs="fs", m=1),
              family = gaussian(), method = "REML", data = prb_dat )

ibi.hr <- gam(HR_u ~ stim_type + ln_IBcv_30.s + IBzu_30 + s(subj_id, probe_num, bs="fs", m=1),
              family = betar(), method = "REML", data = prb_dat )

ibi.fa <- gam(FAR_u ~ stim_type + ln_IBcv_30.s + IBzu_30 + s(subj_id, probe_num, bs="fs", m=1),
              family = betar(), method = "REML", data = prb_dat )

ibi.rs <- gam(RSmu_30 ~ stim_type + ln_IBcv_30.s + IBzu_30 + s(subj_id, probe_num, bs="fs", m=1),
              family = gaussian(), method = "REML", data = filter(prb_dat, vRS_30==T) )

ibi.rv <- gam(RScv_30 ~ stim_type + ln_IBcv_30.s + IBzu_30 + s(subj_id, probe_num, bs="fs", m=1),
              family = Gamma(), method = "REML", data = filter(prb_dat, vRS_30==T) )



### plots
p_size <- 2
s_size <- 1

# probes
col_mstate <- c(rgb(25.5, 153, 51, maxColorValue = 255), 
                rgb(255, 127.5, 25.5, maxColorValue = 255), 
                rgb(102, 165.75, 204, maxColorValue = 255) )

col_vigil <- brewer.pal(n = 9, name = "OrRd")

sm <- prb_dat %>%
  group_by(subj_id, m_state) %>%
  summarise( zu.m = mean(IBzu_30), zu.s = sd(IBzu_30)/sqrt(length(IBzu_30)) ) %>%
  ungroup()

gm <- sm %>%
  group_by(m_state) %>%
  summarise( zu.gm = mean(zu.m), zu.gs = sd(zu.m)/sqrt(length(zu.m)) ) %>%
  ungroup() 


s1 <- ggplot(sm, aes(y=m_state, x=zu.m, group=subj_id) ) + 
  geom_jitter(shape=21, size=p_size, stroke=s_size, height = .2, aes(fill=m_state, alpha = .2)) + 
  geom_errorbarh(data=gm, aes( y=m_state, x=NULL, xmin = zu.gm-zu.gs*1.96, xmax = zu.gm + zu.gs*1.96, group=NULL, height=0), linewidth=1.2) +
  geom_point(data=gm, shape=23, size=p_size+1, stroke=s_size+.5, aes(y = m_state, x = zu.gm, group = NULL, fill=m_state)) + 
  scale_fill_manual( values = col_mstate ) +
  scale_x_continuous(limits = c(-.96, .96), name="Mean IBI (z.u.)", expand = c(0.01,0.01) ) +
  scale_y_discrete(name="Attentional state", limits=rev) +
  annotate("segment", x = -.5, y = 1, xend = -.4, yend = 1) +
  annotate("segment", x = -.5, y = 3, xend = -.4, yend = 3) +
  annotate("segment", x = -.5, y = 3, xend = -.5, yend = 1) +
  annotate("text", x = -.8, y = 2, label = "***", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  annotate("segment", x = .7, y = 3, xend = .8, yend = 3) +
  annotate("segment", x = .7, y = 2, xend = .8, yend = 2) +
  annotate("segment", x = .8, y = 2, xend = .8, yend = 3) +
  annotate("text", x = .85, y = 2.5, label = ".", fontface = "bold", size = 6, vjust = 0.1, hjust = 0) +
  font_theme


sm <- prb_dat %>%
  mutate(vigil.r = as.ordered(abs(vigil-5))) %>%
  group_by(subj_id, vigil.r) %>%
  summarise( zu.m = mean(IBzu_30), zu.s = sd(IBzu_30)/sqrt(length(IBzu_30)) ) %>%
  ungroup() %>%
  filter(!is.na(vigil.r)) 

gm <- sm %>%
  group_by(vigil.r) %>%
  summarise( zu.gm = mean(zu.m), zu.gs = sd(zu.m)/sqrt(length(zu.m)) ) %>%
  ungroup() 

s2 <- ggplot(sm, aes(y=vigil.r, x=zu.m, group=subj_id) ) + 
  geom_jitter(shape=21, size=p_size, stroke=s_size, height = .2, aes(fill=vigil.r, alpha = .2)) + 
  geom_errorbarh(data=gm, aes( y=vigil.r, x=NULL, xmin = zu.gm-zu.gs*1.96, xmax = zu.gm + zu.gs*1.96, group=NULL, height=0), linewidth=1.2) +
  geom_point(data=gm, shape=23, size=p_size+1, stroke=s_size+.5, aes(y = vigil.r, x = zu.gm, group = NULL, fill=vigil.r)) + 
  scale_fill_manual( values = col_vigil[4:7] ) +
  scale_x_continuous(limits = c(-1.45, 1.45), name="Mean IBI (z.u.)", expand = c(0.01,0.01) ) +
  scale_y_discrete(name="Vigilance level" ) +
  annotate("text", x = -1.25, y = 4.1, label = "*", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  font_theme + theme(axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm") ) )


# pupil
sm <- filter(prb_pup, v30==T) %>% 
  group_by(subj_id) %>%
  summarise( pz.m = mean(muPup.z), pz.s = sd(muPup.z)/sqrt(length(muPup.z)),
             cv.m = mean(ln_IBcv_30.s), cv.s = sd(ln_IBcv_30.s)/sqrt(length(ln_IBcv_30.s)),
             zu.m = mean(IBzu_30), zu.s = sd(IBzu_30)/sqrt(length(IBzu_30)) ) %>%
  ungroup()

fs <- get_modelterm(ibi.pz, select = 1, n.grid = length(unique(prb_pup$probe_num)), as.data.frame = T)
muSmoo <- colMeans(matrix(fs$fit, nrow = length(unique(prb_pup$subj_id)), ncol = length(unique(prb_pup$probe_num))), na.rm = T )

em <- as.data.frame(ggemmeans(ibi.pz, "IBzu_30", condition = list(probe_num=which.min( abs(muSmoo) )) ) )

p1 <- ggplot(sm, aes(y=pz.m, x=zu.m, group=subj_id) ) + 
  geom_point(shape=21, size=p_size, stroke=s_size, fill="#8470FF", alpha=.7) + 
  geom_ribbon(data = em, aes(y=predicted, x=x, ymin = conf.low, ymax = conf.high, 
                             group=group), fill="#8470FF", alpha = .4) +
  geom_line(data = em, aes(y=predicted, x=x, group=group), linewidth=1.2) +
  scale_x_continuous(limits = c(-.39, .39), name="Mean IBI (z.u.)", expand = c(0,0) ) +
  scale_y_continuous(limits = c(-.17, .17), name="Pupil size (z.u.)" ) +
  annotate("text", x = -.335, y = .15, label = "*", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  font_theme


em <- as.data.frame(ggemmeans(ibi.pz, "ln_IBcv_30.s", condition = list(probe_num=which.min( abs(muSmoo) )) ) )

p2 <- ggplot(sm, aes(y=pz.m, x=cv.m, group=subj_id) ) +  
  geom_point(shape=21, size=p_size, stroke=s_size, fill="#8470FF", alpha=.7) + 
  geom_ribbon(data = em, aes(y=predicted, x=x, ymin = conf.low, ymax = conf.high, 
                             group=group), fill="#8470FF", alpha = .4) + 
  geom_line(data = em, aes(y=predicted, x=x, group=group), linewidth=1.2) +
  scale_x_continuous(limits = c(-1.75, 1.75), name="IBI coefficient of variation (z.u.)", expand = c(0,0) ) +
  scale_y_continuous(limits = c(-.17, .17), name="Pupil size (z.u.)" ) +
  annotate("text", x = -1.5, y = .15, label = "*", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  font_theme 


# SDT models
sm <- prb_dat %>%
  group_by(subj_id) %>%
  summarise( dp.m = mean(dp), dp.s = sd(dp)/sqrt(length(dp)),
             cn.m = mean(cn), cn.s = sd(cn)/sqrt(length(cn)),
             hr.m = mean(HR_u), hr.s = sd(HR_u)/sqrt(length(HR_u)),
             fa.m = mean(FAR_u), fa.s = sd(FAR_u)/sqrt(length(FAR_u)),
             cv.m = mean(ln_IBcv_30.s), cv.s = sd(ln_IBcv_30.s)/sqrt(length(ln_IBcv_30.s)),
             zu.m = mean(IBzu_30), zu.s = sd(IBzu_30)/sqrt(length(IBzu_30)) ) %>%
  ungroup()


em <- as.data.frame(ggemmeans(ibi.dp, "IBzu_30"))

d1 <- ggplot(sm, aes(y=dp.m, x=zu.m, group=subj_id) ) + 
  geom_point(shape=21, size=p_size, stroke=s_size, fill="#A5D6A7") + 
  geom_ribbon(data = em, aes(y=predicted, x=x, ymin = conf.low, ymax=conf.high, 
                             group=group), fill="#A5D6A7", alpha = .4) + 
  geom_line(data = em, aes(y=predicted, x=x, group=group), linewidth=1.2) +
  scale_x_continuous(limits = c(-.39, .39), name="Mean IBI (z.u.)", expand = c(0,0) ) +
  scale_y_continuous(limits = c(1.4, 2.6), name="Sensitivity (z.u.)" ) +
  annotate("text", x = -.335, y = 2.5, label = "*", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  font_theme


em <- as.data.frame(ggemmeans(ibi.dp, "ln_IBcv_30.s"))

d2 <- ggplot(sm, aes(y=dp.m, x=cv.m, group=subj_id) ) + 
  geom_point(shape=21, size=p_size, stroke=s_size, fill="#A5D6A7") + 
  geom_ribbon(data = em, aes(y=predicted, x=x, ymin = conf.low, ymax=conf.high, 
                             group=group), fill="#A5D6A7", alpha = .4) + 
  geom_line(data = em, aes(y=predicted, x=x, group=group), linewidth=1.2) +
  scale_x_continuous(limits = c(-1.75, 1.75), name="IBI coefficient of variation (z.u.)", expand = c(0,0) ) +
  scale_y_continuous(limits = c(1.4, 2.6), name="Sensitivity (z.u.)" ) +
  annotate("text", x = -1.5, y = 2.5, label = "*", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  font_theme


em <- as.data.frame(ggemmeans(ibi.hr, "IBzu_30"))

h1 <- ggplot(sm, aes(y=hr.m, x=zu.m, group=subj_id) ) + 
  geom_point(shape=21, size=p_size, stroke=s_size, fill="#A5D6A7") + 
  geom_ribbon(data = em, aes(y=predicted, x=x, ymin = conf.low, ymax=conf.high, 
                             group=group), fill="#A5D6A7", alpha = .4) + 
  geom_line(data = em, aes(y=predicted, x=x, group=group), linewidth=1.2) +
  scale_x_continuous(limits = c(-.39, .39), name="Mean IBI (z.u.)", expand = c(0,0) ) +
  scale_y_continuous(limits = c(.89, 1), name="Hit rate", n.breaks = 2, expand = c(0,0) ) +
  annotate("text", x = -.335, y = .987, label = "*", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  font_theme


em <- as.data.frame(ggemmeans(ibi.fa, "ln_IBcv_30.s"))

f1 <- ggplot(sm, aes(y=fa.m, x=cv.m, group=subj_id) ) + 
  geom_point(shape=21, size=p_size, stroke=s_size, fill="#A5D6A7") + 
  geom_ribbon(data = em, aes(y=predicted, x=x, ymin = conf.low, ymax = conf.high, 
                             group=group), fill="#A5D6A7", alpha = .4) + 
  geom_line(data = em, aes(y=predicted, x=x, group=group), linewidth=1.2) +
  scale_x_continuous(limits = c(-1.75, 1.75), name="IBI coefficient of variation (z.u.)", expand = c(0,0) ) +
  scale_y_continuous(limits = c(.1, .7), name="False alarm rate" ) +
  annotate("text", x = -1.5, y = .65, label = "***", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  font_theme


em <- as.data.frame(ggemmeans(ibi.cn, "ln_IBcv_30.s"))

c1 <- ggplot(sm, aes(y=cn.m, x=cv.m, group=subj_id) ) + 
  geom_point(shape=21, size=p_size, stroke=s_size, fill="#A5D6A7") + 
  geom_ribbon(data = em, aes(y=predicted, x=x, ymin = conf.low, ymax = conf.high, 
                             group=group), fill="#A5D6A7", alpha = .4) + 
  geom_line(data = em, aes(y=predicted, x=x, group=group), linewidth=1.2) +
  scale_x_continuous(limits = c(-1.75, 1.75), name="IBI coefficient of variation (z.u.)", expand = c(0,0) ) +
  scale_y_continuous(limits = c(.4, 1), name="Bias (z.u.)" ) +
  annotate("text", x = -1.5, y = .97, label = "**", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  font_theme


# reaction speed
sm <- filter(prb_dat, vRS_30==T) %>%
  group_by(subj_id) %>%
  summarise( rs.m = mean(RSmu_30), rs.s = sd(RSmu_30)/sqrt(length(RSmu_30)),
             cv.m = mean(ln_IBcv_30.s), cv.s = sd(ln_IBcv_30.s)/sqrt(length(ln_IBcv_30.s)),
             zu.m = mean(IBzu_30), zu.s = sd(IBzu_30)/sqrt(length(IBzu_30)) ) %>%
  ungroup()


fs <- get_modelterm(ibi.rs, select = 1, n.grid = length(unique(prb_dat$probe_num)), as.data.frame = T)
muSmoo <- colMeans(matrix(fs$fit, nrow = length(unique(prb_dat$subj_id)), ncol = length(unique(prb_dat$probe_num))), na.rm = T )

em <- as.data.frame(ggemmeans(ibi.rs, "ln_IBcv_30.s", condition = list(probe_num=which.min(abs(muSmoo)) )) )

r1 <- ggplot(sm, aes(y=rs.m, x=cv.m, group=subj_id) ) + 
  geom_point(shape=21, size=p_size, stroke=s_size, fill="#A5D6A7") + 
  geom_ribbon(data = em, aes(y=predicted, x=x, ymin = conf.low, ymax = conf.high, 
                             group=group), fill="#A5D6A7", alpha = .4) + 
  geom_line(data = em, aes(y=predicted, x=x, group=group), linewidth=1.2) +
  scale_x_continuous(limits = c(-1.75, 1.75), name="IBI coefficient of variation (z.u.)", expand = c(0,0) ) +
  scale_y_continuous(limits = c(1.55, 2.45), name=expression(bold( paste("Speed ", (s^{-1}) ) )), breaks = c(1.6, 2, 2.4) ) +
  annotate("text", x = -1.5, y = 2.4, label = "***", fontface = "bold", size = 6, vjust = .75, hjust = 0) +
  font_theme



# arrange in grid for figure 2
plot_grid(NULL,NULL, NULL,
          s1, NULL, s2,
          NULL, NULL, NULL,
          p1, NULL, p2, 
          NULL, NULL, NULL,
          d1, NULL, d2, 
          h1, NULL, f1,
          c1, NULL, r1,
          labels = c("", "", "", "A", "", "", "B", "", "", "", "", "", "C"), vjust = .7,
          rel_heights = c(.03, 1, .01, 1, .1, 1, 1, 1), align = "v", axis = "lr",
          rel_widths = c(1, .1, 1),
          nrow = 8)

ggsave('figure2.png', height = 25, width = 14, units = "cm", dpi=600)


