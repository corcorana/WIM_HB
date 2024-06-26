library(tidyverse)
library(mgcv)
library(mgcViz)
library(itsadug)
library(emmeans)
library(ggeffects)

# global settings
options(contrasts = c("contr.treatment", "contr.poly"))
emm_options(rg.limit = 1e6)

theme_set(theme_classic())
fsize = 10
font_theme <- theme(
  axis.title.x = element_text(size = fsize+2, face="bold"),
  axis.title.y = element_text(size = fsize+2, face="bold"),
  axis.text.x = element_text(size = fsize),
  axis.text.y = element_text(size = fsize),
  strip.text.x = element_text(size = fsize+2, face="bold"),
  strip.text.y = element_text(size = fsize+2, face="bold"),
  legend.position = "none",  
  panel.border = element_rect(fill=NA),
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




### plot figure

s <- summary(ibi.ms)
df <- data.frame( dv = c( "MW", "MW", "MB", "MB" ), 
                  iv = c( "m", "v", "m", "v" ), 
                  grp = c ( "s", "s", "s", "s" ),
                  est = c( s$p.coeff["IBzu_10"], s$p.coeff["ln_IBcv_10.s"], s$p.coeff["IBzu_10.1"], s$p.coeff["ln_IBcv_10.s.1"]), 
                  se = c( s$se["IBzu_10"], s$se["ln_IBcv_10.s"], s$se["IBzu_10.1"], s$se["ln_IBcv_10.s.1"] ) )
s <- summary(ibi.vo)
df <- df %>%
  add_row( dv = "vigil", iv = c("m", "v"), grp = c("s", "s"),
           est = c( s$p.coeff["IBzu_30"], s$p.coeff["ln_IBcv_30.s"]), 
           se = c( s$se["ln_IBcv_30.s"], s$se["ln_IBcv_30.s"] ) )
s <- summary(ibi.pz)
df <- df %>%
  add_row( dv = "pupil", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["IBzu_30"], s$p.coeff["ln_IBcv_30.s"]), 
           se = c( s$se["ln_IBcv_30.s"], s$se["ln_IBcv_30.s"] ) )
s <- summary(ibi.dp)
df <- df %>%
  add_row( dv = "dp", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["IBzu_30"], s$p.coeff["ln_IBcv_30.s"]), 
           se = c( s$se["ln_IBcv_30.s"], s$se["ln_IBcv_30.s"] ) )
s <- summary(ibi.cn)
df <- df %>%
  add_row( dv = "cn", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["IBzu_30"], s$p.coeff["ln_IBcv_30.s"]), 
           se = c( s$se["ln_IBcv_30.s"], s$se["ln_IBcv_30.s"] ) )
s <- summary(ibi.hr)
df <- df %>%
  add_row( dv = "hr", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["IBzu_30"], s$p.coeff["ln_IBcv_30.s"]), 
           se = c( s$se["ln_IBcv_30.s"], s$se["ln_IBcv_30.s"] ) )
s <- summary(ibi.fa)
df <- df %>%
  add_row( dv = "far", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["IBzu_30"], s$p.coeff["ln_IBcv_30.s"]), 
           se = c( s$se["ln_IBcv_30.s"], s$se["ln_IBcv_30.s"] ) )
s <- summary(ibi.rs)
df <- df %>%
  add_row( dv = "mRS", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["IBzu_30"], s$p.coeff["ln_IBcv_30.s"]), 
           se = c( s$se["ln_IBcv_30.s"], s$se["ln_IBcv_30.s"] ) )
s <- summary(ibi.rv)
df <- df %>%
  add_row( dv = "RSv", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["IBzu_30"], s$p.coeff["ln_IBcv_30.s"]), 
           se = c( s$se["ln_IBcv_30.s"], s$se["ln_IBcv_30.s"] ) )


df$dv <- factor(df$dv, 
                levels = c("RSv", "mRS", "far", "hr", "cn", "dp", "pupil", "vigil", "MB", "MW"), 
                labels = c("RS variability", "Response speed", "False alarm rate", "Hit rate", "Bias (c)", "Sensitivity (d')", "Pupil size", "Vigilance", "Mind-blanking", "Mind-wandering") )
df$iv <- factor(df$iv,
                levels = c("m", "v"),
                labels = c("IBI mean", "IBI variability") )
df$grp <- factor(df$grp,
                levels = c("s", "o") )


df_scales <- data.frame(
  Panel = c("m", "v"),
  xmin = c(-.9, -.7),
  xmax = c(.9, .7)
)
df_scales <- split(df_scales, df_scales$Panel)

scalesx <- lapply(df_scales, function(x) {
  scale_x_continuous(limits = c(x$xmin, x$xmax), name = "" )
})


levs <- levels(df$dv[df$grp=="s"])
df_labs <- data.frame(
  Panel = c(rep("m", length(levs)), rep("v", length(levs))),
  nb = c( levs, rep("NULL", length(levs)) )
)
df_labs <- split(df_labs, df_labs$Panel)

scalesy <- lapply(df_labs, function(x) {
  scale_y_discrete(breaks = x$nb, name = "", labels = levs )
})

subj <- ggplot(filter(df, grp=="s"), aes(y = dv, x = est, xmin = est-(1.96*se), xmax = est+(1.96*se), col = dv ) ) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) + 
  geom_point(size = 2) + geom_errorbarh(height=0) + 
  scale_color_manual(values = c("Mind-wandering" = "#ff801a", "Mind-blanking" = "#7aa6cc", "Vigilance" = "#d11141")) +
  facet_wrap(.~iv, scales = "free") +
  ggh4x::facetted_pos_scales(x = scalesx, y = scalesy) +
  font_theme



df_scales <- data.frame(
  Panel = c("m", "v"),
  xmin = c(-.27, -.27),
  xmax = c(.27, .27)
)
df_scales <- split(df_scales, df_scales$Panel)

scalesx <- lapply(df_scales, function(x) {
  scale_x_continuous(limits = c(x$xmin, x$xmax), name = expression(paste(beta, " coefficient")) )
})

levs <- levels(df$dv[df$grp=="o"])
df_labs <- data.frame(
  Panel = c(rep("m", length(levs)), rep("v", length(levs))),
  nb = c( levs, rep("NULL", length(levs)) )
)
df_labs <- split(df_labs, df_labs$Panel)

scalesy <- lapply(df_labs, function(x) {
  scale_y_discrete(breaks = x$nb, name = "", labels = levs )
})

obj <- ggplot(filter(df, grp=="o"), aes(y = dv, x = est, xmin = est-(1.96*se), xmax = est+(1.96*se), col = dv ) ) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) + 
  geom_point(size = 2) + geom_errorbarh(height=0) + 
  scale_color_manual(values = c("Pupil size" = "#8470FF", "Sensitivity (d')" = "#A5D6A7", "Bias (c)" = "#A5D6A7", 
                                "Hit rate" = "#A5D6A7", "False alarm rate" = "#A5D6A7", "Response speed" = "#ffc425", "RS variability" = "#ffc425")) +
  facet_wrap(.~iv, scales = "free") +
  ggh4x::facetted_pos_scales(x = scalesx, y = scalesy) +
  font_theme + theme(strip.text.x = element_blank())
  


cowplot::plot_grid(NULL,
                   subj, 
                   obj,
                   labels = c("", "A", "B"), vjust = .7,
                   rel_heights = c(.01, .5, .7), align = "v", axis = "lr",
                   nrow = 3)

ggsave('figure2_betas.png', height = 12, width = 14, units = "cm", dpi=600, path = 'plots')






