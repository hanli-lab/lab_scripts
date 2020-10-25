library(tidyverse)
library(XML)
library(cowplot)
library(gtools)
library(broom)

parse_xml <- function(xml_path) {

    xml <- xmlParse(xml_path)
    
    abs <- xpathApply(xml, '//Well/RawData', xmlValue)
    time <- xpathApply(xml, '//Well/TimeData', xmlValue)[[1]][1]
    temp <- xpathApply(xml, '//TemperatureData', xmlValue)
    names <- xpathApply(xml, '//Well')
    names <- unlist(map(names, xmlGetAttr, 'Name'))
    
    abs <- str_split(abs, ' ')
    time <- str_split(time, ' ')
    temp <- str_split(temp, ' ')
    
    df <- tibble(temp = unlist(temp),
                 time = unlist(time))
    
    abs_df <- map(abs, tibble) %>%
        bind_cols() %>%
        set_names(names) %>%
        bind_cols(df) %>%
        pivot_longer(!c('temp', 'time'), values_to='abs', names_to='well') %>%
        mutate(time = as.numeric(time), abs = as.numeric(abs)) 
        
    
    return (abs_df)

}

df <-parse_xml('T2.xml')
df$well <- factor(df$well, levels=unique(mixedsort(df$well)))

# plot traces
g <- df %>%
    ggplot(aes(x = time, y=abs, color=well)) +
    #geom_point() +
    geom_line(size=0.9) +
    labs(x='Time (s)', y=expression('OD'[340]), color='Well', title='Cofactor Reduction Assay', subtitle='Trial 2') +
    theme_half_open(12) +
    background_grid() 
    
ggsave('trace.pdf', g, width=6, height=4)

# get slopes from first 60s
slopes <- df%>%
    filter(time < 60) %>%
    nest(data = -well) %>%
    mutate(
        fit = map(data, ~lm(abs ~ time, data=.x)),
        tidied = map(fit, tidy)
    ) %>%
    unnest(tidied) %>%
    filter(term == 'time') %>%
    select(well, estimate) %>%
    rename(slope_a_s = estimate) %>%
    mutate(slope_mA_min = slope_a_s * 1000 * 60)

write.csv(slopes, 'slopes.csv', row.names=FALSE)

###########################################################################

# michaelis-menten fit
rxn_vol <- 100  # ul
enz_vol <- 10
enz_mw <- 35533.5  # g/mol
ext <- 6.22  # mm^-1 cm^-1
dil <- 42.7
conc <- 12.7  # mg/ml
cofa_concs = rep(c(20, 10, 5, 2.5, 1.25, 0.625), 2)

slopes <- bind_cols(slopes, s=cofa_concs)

# abs to velocity
enz_rxn_conc_g_l <- conc * 1/dil * enz_vol/rxn_vol
enz_rxn_uM <- enz_rxn_conc_g_l * 1/enz_mw * 10^6
enz_rxn_mM <- enz_rxn_uM / 1000
pathlen <- rxn_vol / 400

slopes$v0_mM_s = slopes$slope_a_s / (ext*pathlen)

# fit
mm <- formula(v0_mM_s ~ vm * s / (km + s))
fit <- nls(mm, slopes, start=list(km = max(slopes$v0_mM_s)/2, vm = max(slopes$v0_mM_s)))
km <- as.numeric(coef(fit)['km'])
kcat <- as.numeric(coef(fit)['vm']) / enz_rxn_mM

pred_cofa <- seq(min(slopes$s), max(slopes$s), by=0.1)
mm_trace <- predict(fit, newdata=list(s = pred_cofa))
plot_mm <- tibble(trace = mm_trace / enz_rxn_mM,
                  s = pred_cofa)

g2 <- plot_mm %>%
    ggplot(aes(x=s, y=trace)) + 
    geom_line(color='dodgerblue3', size=0.9) +
    labs(x='Substrate (mM)', y=expression(paste('Rate (s'^'-1', ')')), title='Michaelis-Menten Kinetics') +
    theme_half_open(12) +
    background_grid() +
    geom_vline(xintercept=km, color='orange', alpha=0.6, size=0.7) +
    geom_point(data=slopes, aes(x=s, y =v0_mM_s / enz_rxn_mM), color='black', pch=21, fill='green', size=2.2) +
    annotate('text', x=15, y=0.15, label=paste('italic(k)[cat] ==', round(kcat, 3)), size=6, parse=TRUE) +
    annotate('text', x=15, y=0.11, label=paste('italic(K)[m] ==', round(km, 3)), size=6, parse=TRUE)
    
ggsave('kinetics.pdf', g2, width=5.5, height=4)
