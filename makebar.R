qpcr=read.csv(file='/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_6Oct22.csv',
              header=T)
dim(qpcr)
qpcr_sub=subset(qpcr, Buckley > 1.74)
dim(qpcr_sub)

library(tidyverse)
make_bar_relabun <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))
  ggplot(df_out, aes(x = Buckley, y = qPCR_RMAX)) +
    scale_fill_manual(values=c('gray50',"#a44f9a")) +
    scale_colour_manual(values=c("gray50", "#a44f9a")) +
    facet_grid(~TP, scale='free_x', space="free", shrink=TRUE) +
    geom_point(aes(colour=factor(MC), 
                   fill = factor(MC)), size=2) + 
    geom_line(aes(colour=factor(MC))) + 
    labs(x="Buoyant Density", y = "Ratio of Maximum Quantity") +
    theme_bw() +  scale_x_continuous(limits = c(1.76, 1.83)) +
    theme(axis.text = element_text(color ='black',size=12, hjust=1,vjust=1),
          axis.text.x=element_text(size=8, angle=90),
          axis.title = element_text(color='black',face='bold',size=14),
          legend.position = "none",
          panel.grid=element_blank(),
          legend.text = element_text(size=6),
          legend.key.size = unit(0.25,'cm'),
          #strip.text.x = element_text(
          # size = 10, color = "black", face = "bold"),
          #strip.background = element_rect(
          # color="white", fill="white", size=1, linetype="solid"),
          #panel.spacing = unit(0.05, "lines")
    )
}

make_bar_relabun(qpcr_sub, "MC2")
make_bar_relabun(qpcr_sub, "MC3")




make_bar_relabun <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))
  ggplot(df_out, aes(x = Buckley, y = qPCR_RMAX)) +
    scale_fill_manual(values=colls) +
    scale_colour_manual(values=colls) +
    facet_grid(~TP, scale='free_x', space="free", shrink=TRUE) +
    geom_point(aes(colour=MC, 
                   fill = MC), size=1.3) + 
    geom_line(aes(colour=MC)) + 
    labs(x="Buoyant Density (g/mL)", y = "Ratio of Maximum Quantity") +
    theme_bw() +  scale_x_continuous(limits = c(1.76, 1.83)) +
    theme(axis.text = element_text(color ='black',size=6, hjust=1,vjust=1),
          axis.text.x=element_text(size=6, angle=90),
          axis.title = element_text(color='black',face='bold',size=8),
          legend.position = "right",
          panel.grid=element_blank(),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.25,'cm'),
          strip.text.x = element_text(
           size = 10, color = "black", face = "bold"),
          strip.background = element_rect(
           color="gray80", fill="gray97", size=1, linetype="solid"),
          panel.spacing = unit(0.05, "lines")
    )
}


colls=c('gray50',"#a44f9a")
m2=make_bar_relabun(qpcr_sub, "MC2")
mm=make_bar_relabun(qpcr_sub, "MC2")

colls=c('gray50',"#6870c8")
m3=make_bar_relabun(qpcr_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=make_bar_relabun(qpcr_sub, "MC4")
colls=c('gray50',"#af953c")
m5=make_bar_relabun(qpcr_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=make_bar_relabun(qpcr_sub, "MC6")

library(cowplot)
pdf(file='bd_data.pdf', width=10,height=12)
plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
dev.off()




dim(qpcr)
qpcr_sub=subset(qpcr_sub, TP!="3.5")
dim(qpcr_sub)


make_bar_relabun <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))
  ggplot(df_out, aes(x = Buckley, y = qPCR_conc)) +
    scale_fill_manual(values=colls) +
    scale_colour_manual(values=colls) +
    facet_grid(~TP, scales="free_y", space="free", shrink=TRUE) +
    geom_point(aes(colour=MC, 
                   fill = MC), size=1.3) + 
    geom_line(aes(colour=MC)) + 
    labs(x="Buoyant Density (g/mL)", y = "Ratio of Maximum Quantity") +
    theme_bw() +  scale_x_continuous(limits = c(1.76, 1.83)) +
    theme(axis.text = element_text(color ='black',size=6, hjust=1,vjust=1),
          axis.text.x=element_text(size=6, angle=90),
          axis.title = element_text(color='black',face='bold',size=8),
          legend.position = "right",
          panel.grid=element_blank(),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.25,'cm'),
          strip.text.x = element_text(
            size = 10, color = "black", face = "bold"),
          strip.background = element_rect(
            color="gray80", fill="gray97", size=1, linetype="solid"),
          panel.spacing = unit(0.05, "lines")
    )
}


colls=c('gray50',"#a44f9a")
m2=make_bar_relabun(qpcr_sub, "MC2")
colls=c('gray50',"#6870c8")
m3=make_bar_relabun(qpcr_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=make_bar_relabun(qpcr_sub, "MC4")
colls=c('gray50',"#af953c")
m5=make_bar_relabun(qpcr_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=make_bar_relabun(qpcr_sub, "MC6")


pdf(file='bd_data2.pdf', width=10,height=12)
plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
dev.off()

plot_grid(m2, mm, align='hv',ncol=1)



