qpcr=read.csv(file='/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_27Sept22.csv',
              header=T)
dim(qpcr)
qpcr_sub=subset(qpcr, Buckley > 1.74)
dim(qpcr_sub)

s=subset(qpcr_sub, MC=="MC1" | MC== "MC2")

blah=data.frame(tax_table(phyo_abund))
colors=c('gray70',"#a44f9a")
nan=cbind(as.factor(s$Sample), (colors))



dim(s)
library(ggplot2)
p_rgn = ggplot(s, aes(x = Buckley, y = qPCR_RMAX)) +
  scale_fill_manual(values=c('gray50',"#a44f9a")) +
  scale_colour_manual(values=c("gray50", "#a44f9a")) +
  facet_grid(~TP, scale='free_x', space="free", shrink=TRUE) +
  geom_point(aes(colour=factor(MC), 
                 fill = factor(MC)), size=2) + 
  geom_line(aes(colour=factor(MC))) + 
labs(x="Buoyant Density", y = "Ratio of Maximum Quantity") +
  theme_bw() +
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
p_rgn


## goood
p_rgn = ggplot(qpcr_sub, aes(x = Buckley, y = qPCR_RMAX)) +
  scale_fill_manual(values=c('gray50',"#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f")) +
  scale_colour_manual(values=c('gray50',"#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f")) +
  facet_grid(TP ~ MC, scale='fixed', space="free", shrink=TRUE) +
  geom_point(aes(colour=factor(MC), 
                 fill = factor(MC)), size=1.3) + 
  geom_line(aes(colour=factor(MC))) + 
  labs(x="Buoyant Density", y = "Ratio of Maximum Quantity") +
  theme_bw() + scale_x_continuous(limits = c(1.76, 1.83)) +
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
p_rgn

## goood

p_rgn = ggplot(qpcr_sub, aes(x = Buckley, y = qPCR_conc)) +
  scale_fill_manual(values=c('gray50',"#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f")) +
  scale_colour_manual(values=c('gray50',"#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f")) +
  facet_grid(TP ~ MC, scale='fixed', space="free", shrink=TRUE) +
  geom_point(aes(colour=factor(MC), 
                 fill = factor(MC)), size=1.3) + 
  geom_line(aes(colour=factor(MC))) + 
  labs(x="Buoyant Density", y = "Ratio of Maximum Quantity") +
  theme_bw() + scale_x_continuous(limits = c(1.76, 1.83)) +
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
p_rgn




mm <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))
  ggplot(df_out, aes(x = Buckley, y = Ribo_RMAX)) +
    scale_fill_manual(values=colls) +
    scale_colour_manual(values=colls) +
    facet_grid(~TP, scale='free_x', space="free", shrink=TRUE) +
    geom_point(aes(colour=factor(MC), 
                   fill = factor(MC)), size=2) + 
    geom_line(aes(colour=factor(MC))) + 
    labs(x="Buoyant Density (g/mL)", y = "Ratio of Maximum Quantity") +
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
colls=c('gray50',"#af953c")
mm(qpcr_sub, "MC5")

colls=c('gray50',"#ba4a4f")
make_bar_relabun(qpcr_sub,"MC6")
mm(qpcr_sub, "MC6")




