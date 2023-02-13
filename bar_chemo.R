qpcr=read.csv(file='/Users/oliviaahern/Documents/R/Exp4/Dec22/qpcr_chemostat.csv',
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
                   fill = MC), size=1.1) + 
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
            color="gray70", fill="white", size=1, linetype="solid"),
          panel.spacing = unit(0.05, "lines")
    )
}


colls=c('gray50',"#a44f9a")
m2=make_bar_relabun(qpcr_sub, "MC2")
mm=make_bar_relabun(qpcr_sub, "MC2")
mm
colls=c('gray50',"#6870c8")
m3=make_bar_relabun(qpcr_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=make_bar_relabun(qpcr_sub, "MC4")
colls=c('gray50',"#af953c")
m5=make_bar_relabun(qpcr_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=make_bar_relabun(qpcr_sub, "MC6")

lab1=c( "MC2-Met", "MC3-Eth", "MC4-Ace",
       "MC5-Glu", "MC6-Xyl")

col1=c("#a44f9a","#6870c8","#56ae6c","#af953c","#ba4a4f")

library(cowplot)
pdf(file='bd_data_chemo.pdf', width=12,height=12)
plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
dev.off()



library(cowplot)
pdf(file='bd_data_chemo3.pdf', width=12,height=12)
plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1, labels=lab1, 
          label_fontface='plain')
dev.off()




make_bar_relabun1 <- function(df, selection){
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
    theme(axis.text.y = element_text(color ='black',size=8, hjust=1,vjust=1),
         # axis.text.x=element_text(size=6, angle=90),
          axis.title = element_text(color='black',face='bold',size=8),
          legend.position = "right",
          panel.grid=element_blank(),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.25,'cm'),
          strip.text.x = element_text(
            size = 10, color = "black", face = "bold"),
          strip.background = element_rect(
            color="gray80", fill="gray97", size=1, linetype="solid"),
          panel.spacing = unit(0.05, "lines"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
    )
}

make_bar_relabun2 <- function(df, selection){
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
    theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
          axis.text.x=element_text(size=8, angle=90),
          axis.title = element_text(color='black',face='bold',size=8),
          legend.position = "right",
          panel.grid=element_blank(),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.25,'cm'),
          strip.text.x = element_text(
            size = 10, color = "black", face = "bold"),
          strip.background = element_rect(
            color="gray99", fill="gray99", size=1, linetype="solid"),
          panel.spacing = unit(0.05, "lines")
    )
}






colls=c('gray50',"#a44f9a")
m2=make_bar_relabun1(qpcr_sub, "MC2")
m2

colls=c('gray50',"#6870c8")
m3=make_bar_relabun2(qpcr_sub, "MC3")
m3
colls=c('gray50',"#56ae6c")
m4=make_bar_relabun2(qpcr_sub, "MC4")
colls=c('gray50',"#af953c")
m5=make_bar_relabun2(qpcr_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=make_bar_relabun2(qpcr_sub, "MC6")

library(cowplot)
#pdf(file='bd_data_chemo2.pdf', width=8,height=12)
plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
dev.off()


lala <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c(selection))}
c12=subset(qpcr_sub, MC=="MC1")
met=subset(qpcr_sub, MC=="MC2")

lala=subset(qpcr_sub, MC!="MC1")
colls=c('gray50',"#a44f9a","#6870c8","#56ae6c","#af953c","#ba4a4f")

gg1=ggplot() + 
  geom_line(data = c12, aes(x = Buckley, y = qPCR_RMAX, colour=MC)) +
  geom_point(data = c12, aes(x = Buckley, y = qPCR_RMAX, colour=MC), size=1.4) +
  geom_line(data = lala, aes(x = Buckley, y = qPCR_RMAX, colour=MC)) +
  geom_point(data = lala, aes(x = Buckley, y = qPCR_RMAX, colour=MC), size=1.4) +
  scale_fill_manual(values=c('gray50',"#a44f9a","#6870c8","#56ae6c","#af953c","#ba4a4f")) +
  scale_colour_manual(values=c('gray50',"#a44f9a","#6870c8","#56ae6c","#af953c","#ba4a4f")) +
  facet_grid(~TP, scale='free_x', space="free", shrink=TRUE) +
  #geom_point(aes(colour=colls, 
  #               fill = colls), size=1.3) + 
  #geom_line(aes(colour=colls)) + 
  labs(x="Buoyant Density (g/mL)", y = "Ratio of Maximum Quantity") +
  theme_bw() +  scale_x_continuous(limits = c(1.76, 1.83)) +
  theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=8),
        legend.position = "right",
        panel.grid=element_blank(),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.25,'cm'),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
        strip.background = element_rect(
          color="gray99", fill="gray99", size=1, linetype="solid"),
        panel.spacing = unit(0.05, "lines"))
gg1







make_bar_relabun2 <- function(df, selection){
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
    theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
          axis.text.x=element_text(size=8, angle=90),
          axis.title = element_text(color='black',face='bold',size=8),
          legend.position = "right",
          panel.grid=element_blank(),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.25,'cm'),
          strip.text.x = element_text(
            size = 10, color = "black", face = "bold"),
          strip.background = element_rect(
            color="gray99", fill="gray99", size=1, linetype="solid"),
          panel.spacing = unit(0.05, "lines")
    )
}

lala=subset(qpcr_sub, MC!="MC1")

make_bar_relabun2 <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))
  ggplot(df_out, aes(x = Buckley, y = qPCR_RMAX)) +
    scale_fill_manual(values=colls) +
    scale_colour_manual(values=colls) +
    facet_grid(~TP+MC, scale='free_x', space="free", shrink=TRUE) +
    geom_point(aes(colour=MC, 
                   fill = MC), size=1.3) + 
    geom_line(aes(colour=MC)) + 
    labs(x="Buoyant Density (g/mL)", y = "Ratio of Maximum Quantity") +
    theme_bw() +  scale_x_continuous(limits = c(1.76, 1.83)) +
    theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
          axis.text.x=element_text(size=8, angle=90),
          axis.title = element_text(color='black',face='bold',size=8),
          legend.position = "right",
          panel.grid=element_blank(),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.25,'cm'),
          strip.text.x = element_text(
            size = 10, color = "black", face = "bold"),
          strip.background = element_rect(
            color="gray99", fill="gray99", size=1, linetype="solid"),
          panel.spacing = unit(0.05, "lines")
    )
}
colls=c('gray50',"#a44f9a","#6870c8","#56ae6c","#af953c","#ba4a4f")

make_bar_relabun2(qpcr, colls)

# try regular plotting
sub <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))}
f=sub(qpcr_sub, "MC2")    
eth=sub(qpcr_sub, "MC3")    
ace=sub(qpcr_sub, "MC4") 

par(mfrow=c(5,6), mar=c(2,2,1,1))
{plot(subset(f,MC=="MC1" & TP=="36.5")$Buckley, subset(f, MC=="MC1" & TP=="36.5")$qPCR_RMAX,
     type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
lines(subset(f,MC=="MC2" & TP=="36.5")$Buckley, subset(f, MC=="MC2" & TP=="36.5")$qPCR_RMAX,
      col="#a44f9a",lwd=1.4)
axis(2, at=c(0,0.5,1))


plot(subset(f,MC=="MC1" & TP=="37.5")$Buckley, subset(f, MC=="MC1" & TP=="37.5")$qPCR_RMAX,
     type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
lines(subset(f,MC=="MC2" & TP=="37.5")$Buckley, subset(f, MC=="MC2" & TP=="37.5")$qPCR_RMAX,
      col="#a44f9a",lwd=1.4)
axis(2, at=c(0,0.5,1))

plot(subset(f,MC=="MC1" & TP=="38.5")$Buckley, subset(f, MC=="MC1" & TP=="38.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
lines(subset(f,MC=="MC2" & TP=="38.5")$Buckley, subset(f, MC=="MC2" & TP=="38.5")$qPCR_RMAX,
      col="#a44f9a",lwd=1.4)
axis(2, at=c(0,0.5,1))

plot(subset(f,MC=="MC1" & TP=="40.5")$Buckley, subset(f, MC=="MC1" & TP=="40.5")$qPCR_RMAX,
     type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
lines(subset(f,MC=="MC2" & TP=="40.5")$Buckley, subset(f, MC=="MC2" & TP=="40.5")$qPCR_RMAX,
      col="#a44f9a",lwd=1.4)
axis(2, at=c(0,0.5,1))

plot(subset(f,MC=="MC1" & TP=="43.5")$Buckley, subset(f, MC=="MC1" & TP=="43.5")$qPCR_RMAX,
     type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
lines(subset(f,MC=="MC2" & TP=="43.5")$Buckley, subset(f, MC=="MC2" & TP=="43.5")$qPCR_RMAX,
      col="#a44f9a",lwd=1.4)
axis(2, at=c(0,0.5,1))

plot(subset(f,MC=="MC1" & TP=="47.5")$Buckley, subset(f, MC=="MC1" & TP=="47.5")$qPCR_RMAX,
     type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
lines(subset(f,MC=="MC2" & TP=="47.5")$Buckley, subset(f, MC=="MC2" & TP=="47.5")$qPCR_RMAX,
      col="#a44f9a",lwd=1.4)
axis(2, at=c(0,0.5,1))

}


## MC3

{plot(subset(eth,MC=="MC1" & TP=="36.5")$Buckley, subset(eth, MC=="MC1" & TP=="36.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="36.5")$Buckley, subset(eth, MC=="MC3" & TP=="36.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  
  plot(subset(eth,MC=="MC1" & TP=="37.5")$Buckley, subset(eth, MC=="MC1" & TP=="37.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="37.5")$Buckley, subset(eth, MC=="MC3" & TP=="37.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  
  plot(subset(eth,MC=="MC1" & TP=="38.5")$Buckley, subset(eth, MC=="MC1" & TP=="38.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="38.5")$Buckley, subset(eth, MC=="MC3" & TP=="38.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  
  plot(subset(eth,MC=="MC1" & TP=="40.5")$Buckley, subset(eth, MC=="MC1" & TP=="40.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="40.5")$Buckley, subset(eth, MC=="MC3" & TP=="40.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  
  plot(subset(eth,MC=="MC1" & TP=="43.5")$Buckley, subset(eth, MC=="MC1" & TP=="43.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="43.5")$Buckley, subset(eth, MC=="MC3" & TP=="43.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  plot(subset(eth,MC=="MC1" & TP=="47.5")$Buckley, subset(eth, MC=="MC1" & TP=="47.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="47.5")$Buckley, subset(eth, MC=="MC3" & TP=="47.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  
  }

## ace

{plot(subset(ace,MC=="MC1" & TP=="36.5")$Buckley, subset(ace, MC=="MC1" & TP=="36.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="36.5")$Buckley, subset(ace, MC=="MC4" & TP=="36.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  
  plot(subset(ace,MC=="MC1" & TP=="37.5")$Buckley, subset(ace, MC=="MC1" & TP=="37.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="37.5")$Buckley, subset(ace, MC=="MC4" & TP=="37.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  
  plot(subset(ace,MC=="MC1" & TP=="38.5")$Buckley, subset(ace, MC=="MC1" & TP=="38.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="38.5")$Buckley, subset(ace, MC=="MC4" & TP=="38.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  
  plot(subset(ace,MC=="MC1" & TP=="40.5")$Buckley, subset(ace, MC=="MC1" & TP=="40.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="40.5")$Buckley, subset(ace, MC=="MC4" & TP=="40.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  
  plot(subset(ace,MC=="MC1" & TP=="43.5")$Buckley, subset(ace, MC=="MC1" & TP=="43.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="43.5")$Buckley, subset(ace, MC=="MC4" & TP=="43.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  plot(subset(ace,MC=="MC1" & TP=="47.5")$Buckley, subset(ace, MC=="MC1" & TP=="47.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="47.5")$Buckley, subset(ace, MC=="MC4" & TP=="47.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1))
  }

