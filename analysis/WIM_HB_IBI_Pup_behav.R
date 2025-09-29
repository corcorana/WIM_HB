if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, mgcv, mgcViz, itsadug, emmeans, MASS, ggeffects, ggh4x, cowplot)

# global settings
options(contrasts = c("contr.helmert", "contr.sdif"))
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
  panel.border = element_rect(fill=NA),
  panel.background = element_rect(fill = "transparent", colour = NA),
  plot.background = element_rect(fill = "transparent", colour = NA),
  legend.background = element_rect(fill = "transparent")
)



# import & recode variables
dat <- read.csv("WIM_HB_IBI_pup_behav.csv") %>% 
  mutate(HR = (H+.5)/(H+M+1), # go trials
         FAR = (FA+.5)/(FA+CR+1), # nogo trials
         Dp = qnorm(HR) - qnorm(FAR), 
         Cn = rowMeans(cbind(qnorm(HR), qnorm(FAR))),
         site = as.factor(ifelse(subj_id<100, "PBI", "MBI")),
         subj_id = as.factor(subj_id),
         stim_type = factor(stim_type, levels = c("D", "F"), labels = c("Digit", "Face")),
         prPup = ifelse(is.nan(prPup), 0, prPup),
         state.n = ifelse(state>3, 2, state-1),
         state = factor(ifelse(state>3, 3, state), levels = c("3", "2", "1"), labels = c("MB", "MW", "ON")),
         probe_num.s = scale(probe_num),
         vigil = abs(vigil-5), # reverse score so V4 = Extremely Alert
         vigil.s = scale(vigil),
         ln.cvIBI.s = ( log(cvIBI)-mean(log(cvIBI)) ) / sd(log(cvIBI))
         ) %>%
  group_by(subj_id) %>%
  mutate(totMS = sum(!is.na(state))
         ) %>%
  ungroup

# replace NaN with NA
dat[dat=="NaN"] <- NA


# drop subjects who lack state variance
cntMS <- dat %>%
  group_by(subj_id, state, totMS, .drop = FALSE) %>%
  summarise( nMS = n() ) %>%
  ungroup() %>%
  filter(!is.na(state))

excl <- cntMS$subj_id[which(cntMS$nMS==cntMS$totMS)] 
msdat <- dat[!dat$subj_id %in% excl, ]

cntMS$totMS[is.na(cntMS$totMS)] <- 60


# thresholds for epoch-level data
rThresh <- 5 # min number of responses
pThresh <- .5 # min proportion of pupil data


# state onto time / vigil / pupil
tot.ms.m <- gam(list(state.n 
                     ~ stim_type + probe_num.s + s(subj_id, bs="re"),
                     ~ stim_type + probe_num.s + s(subj_id, bs="re")
                     ),
                family = multinom(K=2), method = "REML", data = filter(msdat, site=="MBI") )

tot.ms.p <- gam(list(state.n 
                     ~ stim_type + probe_num.s + s(subj_id, bs="re"),
                     ~ stim_type + probe_num.s + s(subj_id, bs="re")
                     ),
                family = multinom(K=2), method = "REML", data = filter(msdat, site=="PBI") )


vig.ms.m <- gam(list(state.n 
                     ~ stim_type + vigil.s + s(subj_id, probe_num, bs="fs", m=1),
                     ~ stim_type + vigil.s + s(subj_id, probe_num, bs="fs", m=1)
                     ),
                family = multinom(K=2), method = "REML", data = filter(msdat, site=="MBI") )

vig.ms.p <- gam(list(state.n 
                     ~ stim_type + vigil.s + s(subj_id, probe_num, bs="fs", m=1),
                     ~ stim_type + vigil.s + s(subj_id, probe_num, bs="fs", m=1)
                     ),
                family = multinom(K=2), method = "REML", data = filter(msdat, site=="PBI") )


pup.ms.m <- gam(list(state.n 
                     ~ stim_type + zuPup + s(subj_id, probe_num, bs="fs", m=1),
                     ~ stim_type + zuPup + s(subj_id, probe_num, bs="fs", m=1)
                     ),
                family = multinom(K=2), method = "REML", data = filter(msdat, prPup >= pThresh & site=="MBI") )

pup.ms.p <- gam(list(state.n 
                     ~ stim_type + zuPup + s(subj_id, probe_num, bs="fs", m=1),
                     ~ stim_type + zuPup + s(subj_id, probe_num, bs="fs", m=1)
                     ),
                family = multinom(K=2), method = "REML", data = filter(msdat, prPup >= pThresh & site=="PBI") )


# behav onto state
ms.dp.m <- gam(Dp ~ stim_type + state + s(subj_id, probe_num, bs="fs", m=1),
               family = gaussian(), method = "REML", data = filter(msdat, site=="MBI") )

ms.dp.p <- gam(Dp ~ stim_type + state + s(subj_id, probe_num, bs="fs", m=1),
               family = gaussian(), method = "REML", data = filter(msdat, site=="PBI") )


ms.cn.m <- gam(Cn ~ stim_type + state + s(subj_id, probe_num, bs="fs", m=1),
               family = gaussian(), method = "REML", data = filter(msdat, site=="MBI") )

ms.cn.p <- gam(Cn ~ stim_type + state + s(subj_id, probe_num, bs="fs", m=1),
               family = gaussian(), method = "REML", data = filter(msdat, site=="PBI") )


ms.rs.m <- gam(muRS ~ stim_type + state + s(subj_id, probe_num, bs="fs", m=1),
               family = gaussian(), method = "REML", data = filter(msdat, site=="MBI" & nRS >= rThresh) )

ms.rs.p <- gam(muRS ~ stim_type + state + s(subj_id, probe_num, bs="fs", m=1),
               family = gaussian(), method = "REML", data = filter(msdat, site=="PBI" & nRS >= rThresh) )


ms.rv.m <- gam(cvRS ~ stim_type + state + s(subj_id, probe_num, bs="fs", m=1),
               family = tw(), method = "REML", data = filter(msdat, site=="MBI" & nRS >= rThresh) )

ms.rv.p <- gam(cvRS ~ stim_type + state + s(subj_id, probe_num, bs="fs", m=1),
             family = tw(), method = "REML", data = filter(msdat, site=="PBI" & nRS >= rThresh) )


# summary stats
msdat %>% 
  group_by( site, state ) %>%
  summarise( dp = mean(Dp), dps = sd(Dp),
             cn = mean(Cn), cns = sd(Cn),
             )

filter(msdat, nRS >= rThresh) %>% 
  group_by( site, state ) %>%
  summarise( rsm = mean(muRS), rss = sd(muRS),
             cvm = mean(cvRS), cvs = sd(cvRS),
  )


# ms onto IBI
ibi.ms <- gam(list(state.n 
                   ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1),
                   ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1)
                   ),
              family = multinom(K=2), method = "REML", data = msdat )


# vigil onto IBI
ibi.vo <- gam(vigil ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1),
              family = ocat(R=4), method = "REML", data = dat )


# pupil onto IBI
ibi.pz <- gam(zuPup ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1),
              family = gaussian(), method = "REML", data = filter(dat, prPup >= pThresh) )


# behav onto IBI
ibi.dp <- gam(Dp ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1),
              family = gaussian(), method = "REML", data = dat )

ibi.cn <- gam(Cn ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1),
              family = gaussian(), method = "REML", data = dat )

ibi.hr <- gam(HR ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1),
              family = betar(), method = "REML", data = dat )

ibi.fa <- gam(FAR ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1),
              family = betar(), method = "REML", data = dat )

ibi.rs <- gam(muRS ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1),
              family = gaussian(), method = "REML", data = filter(dat, nRS >= rThresh) )

ibi.rv <- gam(cvRS ~ site + stim_type + ln.cvIBI.s + zuIBI + s(subj_id, probe_num, bs="fs", m=1),
              family = tw(), method = "REML", data = filter(dat, nRS >= rThresh) )



# plot model coefficients

s <- summary(ibi.ms)
df <- data.frame( dv = c( "MW", "MW", "MB", "MB" ), 
                  iv = c( "m", "v", "m", "v" ), 
                  grp = c ( "s", "s", "s", "s" ),
                  est = c( s$p.coeff["zuIBI"], s$p.coeff["ln.cvIBI.s"], s$p.coeff["zuIBI.1"], s$p.coeff["ln.cvIBI.s.1"]), 
                  se = c( s$se["zuIBI"], s$se["ln.cvIBI.s"], s$se["zuIBI.1"], s$se["ln.cvIBI.s.1"] ) )
s <- summary(ibi.vo)
df <- df %>%
  add_row( dv = "vigil", iv = c("m", "v"), grp = c("s", "s"),
           est = c( s$p.coeff["zuIBI"], s$p.coeff["ln.cvIBI.s"]), 
           se = c( s$se["zuIBI"], s$se["ln.cvIBI.s"] ) )
s <- summary(ibi.pz)
df <- df %>%
  add_row( dv = "pupil", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["zuIBI"], s$p.coeff["ln.cvIBI.s"]), 
           se = c( s$se["zuIBI"], s$se["ln.cvIBI.s"] ) )
s <- summary(ibi.dp)
df <- df %>%
  add_row( dv = "dp", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["zuIBI"], s$p.coeff["ln.cvIBI.s"]), 
           se = c( s$se["zuIBI"], s$se["ln.cvIBI.s"] ) )

s <- summary(ibi.cn)
df <- df %>%
  add_row( dv = "cn", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["zuIBI"], s$p.coeff["ln.cvIBI.s"]), 
           se = c( s$se["zuIBI"], s$se["ln.cvIBI.s"] ) )

s <- summary(ibi.hr)
df <- df %>%
  add_row( dv = "hr", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["zuIBI"], s$p.coeff["ln.cvIBI.s"]), 
           se = c( s$se["zuIBI"], s$se["ln.cvIBI.s"] ) )

s <- summary(ibi.fa)
df <- df %>%
  add_row( dv = "far", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["zuIBI"], s$p.coeff["ln.cvIBI.s"]), 
           se = c( s$se["zuIBI"], s$se["ln.cvIBI.s"] ) )

s <- summary(ibi.rs)
df <- df %>%
  add_row( dv = "mRS", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["zuIBI"], s$p.coeff["ln.cvIBI.s"]), 
           se = c( s$se["zuIBI"], s$se["ln.cvIBI.s"] ) )

s <- summary(ibi.rv)
df <- df %>%
  add_row( dv = "RSv", iv = c("m", "v"), grp = c("o", "o"),
           est = c( s$p.coeff["zuIBI"], s$p.coeff["ln.cvIBI.s"]), 
           se = c( s$se["zuIBI"], s$se["ln.cvIBI.s"] ) )


df$dv <- factor(df$dv, 
                levels = c("RSv", "mRS", "far", "hr", "cn", "dp", "pupil", "vigil", "MB", "MW"), 
                labels = c("Response speed\n variability", "Response speed", "False alarm rate", "Hit rate", "Bias (c)", "Sensitivity (d')", "Pupil size", "Vigilance", "Mind-blanking", "Mind-wandering") )
df$iv <- factor(df$iv,
                levels = c("m", "v"),
                labels = c("IBI duration", "IBI variability") )
df$grp <- factor(df$grp,
                levels = c("s", "o") )


df_scales <- data.frame(
  Panel = c("m", "v"),
  xmin = c(-.87, -.87),
  xmax = c(.87, .87)
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
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#53565b") + 
  geom_point(size = 2) + geom_errorbarh(height=0) + 
  scale_color_manual(values = c("Mind-wandering" = "#ff801a", "Mind-blanking" = "#7aa6cc", "Vigilance" = "#d11141")) +
  facet_wrap(.~iv, scales = "free") +
  facetted_pos_scales(x = scalesx, y = scalesy) +
  font_theme + 
  theme(legend.position = "none", 
        plot.margin = unit(c(.5,3,0,1), 'cm'))

df_scales <- data.frame(
  Panel = c("m", "v"),
  xmin = c(-.24, -.24),
  xmax = c(.24, .24)
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
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#53565b") + 
  geom_point(size = 2) + geom_errorbarh(height=0) + 
  scale_color_manual(values = c("Pupil size" = "#8470FF", "Sensitivity (d')" = "#A5D6A7", "Bias (c)" = "#A5D6A7", 
                                "Hit rate" = "#A5D6A7", "False alarm rate" = "#A5D6A7", "Response speed" = "#ffc425", "Response speed\n variability" = "#ffc425")) +
  facet_wrap(.~iv, scales = "free") +
  facetted_pos_scales(x = scalesx, y = scalesy) +
  font_theme + 
  theme(strip.text.x = element_blank(), legend.position = "none", 
        plot.margin = unit(c(.5,3,.5,1), 'cm'))


plot_grid(subj, obj,
          labels = c("A", "B"), 
          rel_heights = c(.45, .7), align = "v", axis = "lr",
          nrow = 2)


ggsave('figure2_betas.tiff', height = 14, width = 18, units = "cm", dpi=600, path = 'plots')


### NoGo IBI analysis

# import & recode variables
dat <- read.csv("WIM_HB_IBI_NoGo.csv") %>% 
  mutate(site = as.factor(ifelse(subj_id<300, "PBI", "MBI")),
         subj_id = as.factor(subj_id),
         stim_type = factor(stim_type, levels = c("D", "F"), labels = c("Digit", "Face")),
         resp_type = factor(resp_type, levels = c("C", "F"), labels = c("CR", "FA")),
         state = factor(ifelse(state>3, 3, state), levels = c("3", "2", "1"), labels = c("MB", "MW", "ON"))
  ) %>% 
  pivot_longer(starts_with("IBI"), names_to = "index", values_to = "zuIBI") %>%
  mutate(index = factor(index, levels = c("IBI_3", "IBI_2", "IBI_1", "IBI_0", "IBI.1", "IBI.2", "IBI.3"),
                        labels = c("IBI-3", "IBI-2", "IBI-1", "IBI 0", "IBI+1", "IBI+2", "IBI+3") ),
         index.o = as.ordered(index)
  )

# replace NaN with NA
dat[dat=="NaN"] <- NA

# drop subjects who lack state variance
cntMS <- dat %>% 
  group_by(subj_id, state, .drop = FALSE) %>%
  summarise( n = n() ) %>%
  ungroup()

excl <- cntMS$subj_id[which(cntMS$n==60)] 
msdat <- dat[!dat$subj_id %in% excl, ]


# IBI duration onto IBI sequence, response, state
nogo <- gam(zuIBI ~ site + stim_type + index.o * resp_type * state + s(subj_id, probe_num, bs="fs", m=1),
            family = gaussian(), method = "REML", data = msdat )


# stats
anova.gam(nogo)


# plot IBIs
df <- as.data.frame(ggemmeans(nogo, c("index.o", "resp_type", "state")))
df$IBI <- recode(df$x, 'IBI-3'='-3', 'IBI-2'='-2', 'IBI-1'='-1', 'IBI 0'='0', 'IBI+1'='+1', 'IBI+2'='+2', 'IBI+3'='+3')
df$group <- recode(df$group, CR = "Correct rejection", FA = "False alarm")
df$facet <- recode(df$facet, ON = "On-task", MW = "Mind-wandering", MB = "Mind-blanking" )
strip <- strip_themed(background_x = elem_list_rect(fill = c("#199933", "#ff801a", "#7aa6cc")))

ggplot( df, aes(x=IBI, y=predicted, ymax=predicted+std.error, ymin=predicted-std.error, col=group, shape=group) ) +
  geom_point(size=2, position=position_dodge(width=0.4)) + 
  geom_errorbar(width=0, position=position_dodge(width=0.4)) +
  labs(y = "IBI duration (z.u.)", 
       x = "IBI order (relative to NoGo trial)") +
  scale_colour_brewer(palette = "Set1", direction = -1) +
  ylim(-.2, .6) +
  facet_wrap2(.~factor(facet, levels = c('On-task', 'Mind-wandering', 'Mind-blanking')), strip=strip) +
  font_theme +
  theme(legend.position = c(0.15, 0.85), legend.title=element_blank() )

ggsave('figure3_emms.png', height = 7, width = 14, units = "cm", dpi = 600, path = 'plots')


# boxplots for Comms Bio
subj_dat <- msdat %>% drop_na() %>%
  group_by(subj_id, state, index.o, resp_type) %>% 
  summarise(muIBI = mean(zuIBI, na.rm=T) ) %>%
  mutate(IBI = recode(index.o, 'IBI-3'='-3', 'IBI-2'='-2', 'IBI-1'='-1', 'IBI 0'='0', 'IBI+1'='+1', 'IBI+2'='+2', 'IBI+3'='+3'),
         group = recode(resp_type, CR = "Correct rejection", FA = "False alarm"),
         facet = recode(state, ON = "On-task", MW = "Mind-wandering", MB = "Mind-blanking" )
         )

ggplot( subj_dat, aes(x=IBI, y=muIBI, colour=group, shape=group) ) +
  geom_boxplot(position = position_dodge(width=.5), alpha = .5) + 
  geom_point(data = subset(subj_dat, muIBI>5)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#53565b") + 
  labs(y = "IBI duration (z.u.)", 
       x = "IBI order (relative to NoGo trial)" ) +
  scale_colour_brewer(palette = "Set1", direction = -1) +
  ylim(-1.7, 3.4) +
  facet_wrap2(.~factor(facet, levels = c('On-task', 'Mind-wandering', 'Mind-blanking')), strip=strip)  +
  font_theme +
  theme(legend.position = c(0.105, 0.84), legend.title=element_blank() )

ggsave('figure3_bps.tiff', height = 7, width = 18, units = "cm", dpi = 600, path = 'plots')
