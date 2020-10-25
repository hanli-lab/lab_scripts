library(tidyverse)
library(ggsci)
library(cowplot)

# read data and convert to tidy format
# add name labels
# correct for dilution and absorbance blank
df <- read_csv('reads_wide.csv') %>%
    pivot_longer(!read_num, names_to='time', values_to='abs') %>%
    mutate(enzyme = case_when(
        read_num %in% seq(1, 3) ~ 'XenA',
        read_num %in% seq(4, 6) ~ 'Blank',
        read_num %in% seq(7, 9) ~ 'Lp NOX WT',
        read_num %in% seq(10, 12) ~ 'Lp NOX I158S-G178Q',
        read_num > 12 ~ 'Library A')) %>%
    mutate(abs_corr = (abs - 0.039) * 10 * 2.82) %>%
    mutate(time = as.numeric(time))


# agg data over replicates to get mean and std at each timepoint
df_grp <- df %>%
    group_by(enzyme, time) %>%
    summarize(abs_corr_mean = mean(abs_corr), abs_corr_std = sd(abs_corr))

# save cleaned data
write.csv(df, 'clean_timeseries.csv', row.names=FALSE)
write.csv(df_grp, 'clean_timeseries_agg.csv', row.names=FALSE)

##################################################################################


# plot aggregated traces on single plot for direct comparison
p <- df_grp %>%
  ggplot(aes(x=time, y=abs_corr_mean, color=enzyme)) +
  geom_point() +
  geom_line(size=0.8) +
  geom_errorbar(aes(ymin = abs_corr_mean - abs_corr_std, ymax = abs_corr_mean + abs_corr_std), width=1.2) +
  theme_half_open(12) +
  background_grid() + 
  #theme(legend.position = 'bottom') +
  labs(title='Lp NOX Growth Selection', subtitle='Library A', x='Time (hr)', y=expression("OD"[600]),
       col='Enzyme') +
  scale_color_npg()
 
ggsave('df_grp.pdf', width=6, height=4)

# plot aggregrated traces overlaid on each raw time series
# facet by enzyme
g <- df %>%
  ggplot(aes(x=time, y=abs_corr)) +
  geom_line(aes(group=read_num), color='black', alpha=0.5) +
  geom_point(data = df_grp, aes(x=time, y=abs_corr_mean, color=enzyme)) +
  geom_line(data = df_grp, aes(x=time, y=abs_corr_mean, color=enzyme), size=1) + 
  geom_errorbar(data=df_grp, aes(x=time, y=abs_corr_mean, ymin=abs_corr_mean - abs_corr_std, ymax = abs_corr_mean + abs_corr_std,
                                 color=enzyme), width=3) + 
  facet_wrap(~enzyme) +
  theme_half_open(15) +
  background_grid() +
  theme(legend.position='none') +
  theme(strip.background = element_rect(fill='gray92')) + 
  theme(plot.margin = margin(10, 20, 10, 10)) + 
  labs(title='Lp NOX Growth Selection', subtitle='Raw Traces', x='Time (hr)', y=expression("OD"[600])) +
  scale_color_npg() 

ggsave('df_raw.pdf', width=10, height=6.5)

    