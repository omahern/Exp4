---
title: "Experiment 4"
output: 
  html_document:
    keep_md: true
    toc: true
    toc_float: true
    toc_depth: 6
    code_folding: hide
    number_sections: false
    theme: cosmo

knit: (function(input_file, encoding) {
  out_dir <- 'docs';
  rmarkdown::render(input_file,
 encoding=encoding,
 output_file=file.path(dirname(input_file), out_dir, 'index.html'))})
---


<video width="600" height="300" controls autoplay muted>
  <source src="/Users/oliviaahern/Documents/movie.mov" type="video/mov">
</video>

# Experimental Overview

Full experimental information, including previous experiments and timeline, can be found at [the MCP website](http://ecosystems.mbl.edu/MEP-FoodWeb/Experiments/Exp4/index.html)

RNA-SIP protocol can be found [here on github](https://omahern.github.io/RNA_SIP_Protocol/). Raw files used in this page can be found [here on github](https://github.com/omahern/Exp4/tree/main/data). 



```r
library(phyloseq)
library(reshape2)
library(tidyverse)
library(vegan)
library(HTSSIP)
library(ape)
library(CoDaSeq)
library(philr)
library(ggtree)
library(cowplot)
library(ggplot2)
library(viridis)
library(phytools)
```

# State Data

## Gas Data

Full gas data for the entire experiment can be found at [the Experiment 4 highcharts website](http://ecosystems.mbl.edu/MEP-FoodWeb/Experiments/Exp4/Exp4_HighCharts.html)


```r
x<-read.csv(file='/Users/oliviaahern/Documents/R/Exp4/22Apr22/Exp4_Chemostat.csv',header=TRUE)
timec=x$Time..d.
mc1cO2c=x$CO2.MC1....
mc2cO2c=x$CO2.MC2....
mc3cO2c=x$CO2.MC3....
mc4cO2c=x$CO2.MC4....
mc5cO2c=x$CO2.MC5....
mc6cO2c=x$CO2.MC6....
mc1O2c=x$O2.MC1....
mc2O2c=x$O2.MC2....
mc3O2c=x$O2.MC3....
mc4O2c=x$O2.MC4....
mc5O2c=x$O2.MC5....
mc6O2c=x$O2.MC6....


x<-read.csv(file='/Users/oliviaahern/Documents/R/Exp4/22Apr22/Exp4_Batch.csv',header=TRUE)
time=x$Time..d.

mc1cO2=x$CO2.MC1....
mc2cO2=x$CO2.MC2....
mc3cO2=x$CO2.MC3....
mc4cO2=x$CO2.MC4....
mc5cO2=x$CO2.MC5....
mc6cO2=x$CO2.MC6....
mc1O2=x$O2.MC1....
mc2O2=x$O2.MC2....
mc3O2=x$O2.MC3....
mc4O2=x$O2.MC4....
mc5O2=x$O2.MC5....
mc6O2=x$O2.MC6....
{
par(mfrow=c(2,2),mar=c(5,5,1,1))
layout(matrix(c(1,2,2,3,4,4), nrow = 2, ncol = 3, byrow = TRUE))

plot(time,mc1O2, type='l',pch=21, col='gray70',ylim=c(20.8,21.2),
     xlab="Time (days)", lwd=1, cex.axis=1.5,cex.lab = 1.5,
     ylab="Oxygen (%)", yaxt='n')
axis(2, at=c(20.8, 21, 21.2), cex.axis=1.5)
lines(time, mc2O2, type='l', pch=21, col='#a44f9a')
lines(time, mc3O2, type='l', pch=21, col='#6870c8')
lines(time, mc4O2, type='l', pch=21, col='#56ae6c')
lines(time, mc5O2, type='l', pch=21, col='#af953c')
lines(time, mc6O2, type='l', pch=21, col='#ba4a4f')

plot(timec,mc1O2c, type='l',pch=21, col='gray70',ylim=c(20.5,21),
     xlab="Time (days)", lwd=1,cex.axis=1.5,cex.lab = 1.5,
     ylab="Oxygen (%)", yaxt="n")
axis(2, at=c(20.6, 20.8,21), cex.axis=1.5)
lines(timec, mc2O2c, type='l', pch=21, col='#a44f9a')
lines(timec, mc3O2c, type='l', pch=21, col='#6870c8')
lines(timec, mc4O2c, type='l', pch=21, col='#56ae6c')
lines(timec, mc5O2c, type='l', pch=21, col='#af953c')
lines(timec, mc6O2c, type='l', pch=21, col='#ba4a4f')

plot(time,mc1cO2, type='l',pch=21, col='gray70',ylim=c(0,0.2),
     xlab="Time (days)", lwd=1,cex.axis=1.5,cex.lab = 1.5,
     ylab="Carbon Dioxide (%)", yaxt='n')
axis(2, at=c(0,0.1,0.2), cex.axis=1.5)
lines(time, mc2cO2, type='l', pch=21, col='#a44f9a')
lines(time, mc3cO2, type='l', pch=21, col='#6870c8')
lines(time, mc4cO2, type='l', pch=21, col='#56ae6c')
lines(time, mc5cO2, type='l', pch=21, col='#af953c')
lines(time, mc6cO2, type='l', pch=21, col='#ba4a4f')

plot(timec,mc1cO2c, type='l',pch=21, col='gray70',ylim=c(0.18,0.3),
     xlab="Time (days)", lwd=1,cex.axis=1.5,cex.lab = 1.5,
     ylab="Carbon Dioxide (%)", yaxt="n")
axis(2, at=c(0.2, 0.25,0.3), cex.axis=1.5)
lines(timec, mc2cO2c, type='l', pch=21, col='#a44f9a')
lines(timec, mc3cO2c, type='l', pch=21, col='#6870c8')
lines(timec, mc4cO2c, type='l', pch=21, col='#56ae6c')
lines(timec, mc5cO2c, type='l', pch=21, col='#af953c')
lines(timec, mc6cO2c, type='l', pch=21, col='#ba4a4f')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/gas-1.png)<!-- -->

## Atomic % 13C

Atomic % 13C was measured from dried GFF filters. Calculations and measurements were carried out by the [MBL Stable Isotope](https://www.mbl.edu/research/research-centers/ecosystems-center/stable-isotope-laboratory)


```r
read=read.csv(file="/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/SIP_Spins/atomic_C_30Aug22.csv",header=T)

read1=subset(read, B_C=="B")

MC1=subset(read1, MC=="MC1")
MC2=subset(read1, MC=="MC2")
MC3=subset(read1, MC=="MC3")
MC4=subset(read1, MC=="MC4")
MC5=subset(read1,MC=="MC5")
MC6=subset(read1,MC=="MC6")

read2=subset(read, B_C=="C")

MC1c=subset(read2, MC=="MC1")
MC2c=subset(read2, MC=="MC2")
MC3c=subset(read2, MC=="MC3")
MC4c=subset(read2, MC=="MC4")
MC5c=subset(read2,MC=="MC5")
MC6c=subset(read2,MC=="MC6")

{par(mar=c(5,5,1,1), mfrow=c(1,2))
plot(MC1$TP,MC1$Balanced_Atomic_.13C, type='o',  ylim=c(0,30),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1.5, cex.lab=1.4)
lines(MC1$TP, MC2$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)
lines(MC1$TP, MC3$Balanced_Atomic_.13C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.5)
lines(MC1$TP, MC4$Balanced_Atomic_.13C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.5)
lines(MC1$TP, MC5$Balanced_Atomic_.13C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.5)
lines(MC1$TP, MC6$Balanced_Atomic_.13C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.5)
legend('topleft',legend=c("MC1-C12", "MC2-Met", "MC3-Eth", "MC4-Ace",
                          "MC5-Glu", "MC6-Xyl"), bty='n', pch =21,
       pt.bg=c("gray70", "#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f"))


plot(MC1c$TP,MC1c$Balanced_Atomic_.13C, type='o',  ylim=c(0,30),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1.5, cex.lab=1.4)
lines(MC1c$TP, MC2c$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)
lines(MC1c$TP, MC3c$Balanced_Atomic_.13C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.5)
lines(MC1c$TP, MC4c$Balanced_Atomic_.13C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.5)
lines(MC1c$TP, MC5c$Balanced_Atomic_.13C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.5)
lines(MC1c$TP, MC6c$Balanced_Atomic_.13C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.5)
legend('topleft',legend=c("MC1-C12", "MC2-Met", "MC3-Eth", "MC4-Ace",
                          "MC5-Glu", "MC6-Xyl"), bty='n', pch =21,
       pt.bg=c("gray70", "#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f"),ncol=2)

}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/atomc2-1.png)<!-- -->


## Atomic % 13C + POC

Atomic % 13C was measured from dried GFF filters. Calculations and measurements were carried out by the [MBL Stable Isotope](https://www.mbl.edu/research/research-centers/ecosystems-center/stable-isotope-laboratory). POC was normalized to the amount of water filtered through each sample. 


```r
{#par(mar=c(5,5,1,1), mfrow=c(2,2))
  layout(matrix(c(1,2,2,3,4,4), nrow = 2, ncol = 3, byrow = TRUE))
  
  plot(MC1$TP,MC1$Balanced_Atomic_.13C, type='o',  ylim=c(0,30),pch=21,bg='gray70',
       xlab="Time (days)", ylab="Atomic %13C", cex=1.2, cex.lab=1.2,
       yaxt='n')
  axis(2,at=c(0,10,20,30))
  lines(MC1$TP, MC2$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.2)
  lines(MC1$TP, MC3$Balanced_Atomic_.13C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.2)
  lines(MC1$TP, MC4$Balanced_Atomic_.13C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.2)
  lines(MC1$TP, MC5$Balanced_Atomic_.13C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.2)
  lines(MC1$TP, MC6$Balanced_Atomic_.13C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.2)
 # legend('topleft',legend=c("MC1-C12", "MC2-Met", "MC3-Eth", "MC4-Ace",
  #                          "MC5-Glu", "MC6-Xyl"), bty='n', pch =21,
  #       pt.bg=c("gray70", "#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f"))
  
  
  plot(MC1c$TP,MC1c$Balanced_Atomic_.13C, type='o',  ylim=c(0,30),pch=21,bg='gray70',
       xlab="Time (days)", ylab="Atomic %13C", cex=1.2, cex.lab=1.2,
       yaxt='n')
  axis(2,at=c(0,10,20,30))
  lines(MC1c$TP, MC2c$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.2)
  lines(MC1c$TP, MC3c$Balanced_Atomic_.13C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.2)
  lines(MC1c$TP, MC4c$Balanced_Atomic_.13C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.2)
  lines(MC1c$TP, MC5c$Balanced_Atomic_.13C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.2)
  lines(MC1c$TP, MC6c$Balanced_Atomic_.13C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.2)
  legend('topleft',legend=c("MC1-C12", "MC2-Met", "MC3-Eth", "MC4-Ace",
                            "MC5-Glu", "MC6-Xyl"), bty='n', pch =21,
         pt.bg=c("gray70", "#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f"),ncol=2)

  
  plot(MC1$TP,MC1$POC.._M., type='o',  ylim=c(0,500),pch=21,bg='gray70',
       xlab="Time (days)", ylab="POC (uM)", cex=1.2, cex.lab=1.2,
       yaxt='n')
  axis(2,at=c(0,250,500))
  lines(MC1$TP, MC2$POC.._M., type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.2)
  lines(MC1$TP, MC3$POC.._M., type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.2)
  lines(MC1$TP, MC4$POC.._M., type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.2)
  lines(MC1$TP, MC5$POC.._M., type="o", col="#af953c", bg="#af953c", pch=21,cex=1.2)
  lines(MC1$TP, MC6$POC.._M., type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.2)

  
  plot(MC1c$TP,MC1c$POC.._M., type='o',  ylim=c(0,1100),pch=21,bg='gray70',
       xlab="Time (days)", ylab="POC (uM)", cex=1.2, cex.lab=1.2,
       yaxt='n')
  axis(2,at=c(0,500,1000))
  lines(MC1c$TP, MC2c$POC.._M., type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.2)
  lines(MC1c$TP, MC3c$POC.._M., type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.2)
  lines(MC1c$TP, MC4c$POC.._M., type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.2)
  lines(MC1c$TP, MC5c$POC.._M., type="o", col="#af953c", bg="#af953c", pch=21,cex=1.2)
  lines(MC1c$TP, MC6c$POC.._M., type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.2)
  
  
  }
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/atomc3-1.png)<!-- -->


# SIP Spins

SIP gradients were setup using a CsTFA gradient centered at ~1.80 g/mL, formamide, and between 750-500 ng of RNA. Samples were spun using an ultracentrifuge spinning at 65,000 rpm at 20C for 65 hours. 

Correction factors include Lueder's Buoyant density (Lueders et al, 2010) and the [Buckley Lab RNA SIP protocol](https://github.com/buckleylab/Buckley_Lab_SIP_project_protocols/blob/master/RNA_SIP/RNA_SIP.md).

Standard deviation between SIP runs (as of Sept 9th, 2022) were:

* Lueders BD and correction factor: 0.004110266 +/- 0.001662171 g/mL (~40% standard deviation)
* Buckley BD and correction factor: 0.00309178 +/- 0.001010027 g/mL (~33% standard deviation)

C13 peak was about 0.016 g/mL higher than the 13C peak. 

## Batch 

```r
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_6Oct22.csv',header=T)
```

### Ratio of Maximum Quantity vs. Buckley Buoyant Density 

Average peak 12C density for MC1-C12

* Lueders BD and correction factor: 1.7811 +/- 0.0090 g/mL 
* Buckley BD and correction factor: 1.7860 +/- 0.0023 g/mL


```r
qpcr=read.csv(file='/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_6Oct22.csv',
              header=T)
# dim(qpcr)
qpcr_sub=subset(qpcr, Buckley > 1.74)

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
plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_a-1.png)<!-- -->

### Ratio of Maximum Quantity vs. Leuders Buoyant Density 



```r
qpcr_leud_rmax <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))
  ggplot(df_out, aes(x = Corrected, y = qPCR_RMAX)) +
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
m2=qpcr_leud_rmax(qpcr_sub, "MC2")
colls=c('gray50',"#6870c8")
m3=qpcr_leud_rmax(qpcr_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=qpcr_leud_rmax(qpcr_sub, "MC4")
colls=c('gray50',"#af953c")
m5=qpcr_leud_rmax(qpcr_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=qpcr_leud_rmax(qpcr_sub, "MC6")

plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_b-1.png)<!-- -->

### RNA Concentration vs. Buckley Buoyant Density 



```r
qpcr_buck_conc <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))
  ggplot(df_out, aes(x = Buckley, y = qPCR_conc)) +
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
m2=qpcr_buck_conc(qpcr_sub, "MC2")
colls=c('gray50',"#6870c8")
m3=qpcr_buck_conc(qpcr_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=qpcr_buck_conc(qpcr_sub, "MC4")
colls=c('gray50',"#af953c")
m5=qpcr_buck_conc(qpcr_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=qpcr_buck_conc(qpcr_sub, "MC6")

plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_a_c-1.png)<!-- -->


### RNA Concentration vs. Leuders Buoyant Density 



```r
qpcr_leud_conc <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))
  ggplot(df_out, aes(x = Corrected, y = qPCR_conc)) +
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
m2=qpcr_leud_conc(qpcr_sub, "MC2")
colls=c('gray50',"#6870c8")
m3=qpcr_leud_conc(qpcr_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=qpcr_leud_conc(qpcr_sub, "MC4")
colls=c('gray50',"#af953c")
m5=qpcr_leud_conc(qpcr_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=qpcr_leud_conc(qpcr_sub, "MC6")

plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_a_d-1.png)<!-- -->


## Chemostat RT-qPCR

#### Ratio of Maximum Quantity vs. Buckley Buoyant Density 


```r
qpcr_chemo=read.csv(file='/Users/oliviaahern/Documents/R/Exp4/Dec22/qpcr_chemostat.csv',
              header=T)
qpcr_chemo_sub=subset(qpcr_chemo, Buckley > 1.74)

colls=c('gray50',"#a44f9a")
m2=make_bar_relabun(qpcr_chemo_sub, "MC2")
colls=c('gray50',"#6870c8")
m3=make_bar_relabun(qpcr_chemo_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=make_bar_relabun(qpcr_chemo_sub, "MC4")
colls=c('gray50',"#af953c")
m5=make_bar_relabun(qpcr_chemo_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=make_bar_relabun(qpcr_chemo_sub, "MC6")
library(cowplot)
plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_chemo_BRMAX-1.png)<!-- -->


#### Ratio of Maximum Quantity vs. Leuders Buoyant Density 


```r
colls=c('gray50',"#a44f9a")
m2=qpcr_leud_rmax(qpcr_chemo_sub, "MC2")
colls=c('gray50',"#6870c8")
m3=qpcr_leud_rmax(qpcr_chemo_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=qpcr_leud_rmax(qpcr_chemo_sub, "MC4")
colls=c('gray50',"#af953c")
m5=qpcr_leud_rmax(qpcr_chemo_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=qpcr_leud_rmax(qpcr_chemo_sub, "MC6")

plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_chemo_LRMAX-1.png)<!-- -->


#### RNA Concentration vs. Buckley Buoyant Density 



```r
colls=c('gray50',"#a44f9a")
m2=qpcr_buck_conc(qpcr_chemo_sub, "MC2")
colls=c('gray50',"#6870c8")
m3=qpcr_buck_conc(qpcr_chemo_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=qpcr_buck_conc(qpcr_chemo_sub, "MC4")
colls=c('gray50',"#af953c")
m5=qpcr_buck_conc(qpcr_chemo_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=qpcr_buck_conc(qpcr_chemo_sub, "MC6")

plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_chemo_BRNA-1.png)<!-- -->


#### RNA Concentration vs. Leuders Buoyant Density 



```r
colls=c('gray50',"#a44f9a")
m2=qpcr_leud_conc(qpcr_chemo_sub, "MC2")
colls=c('gray50',"#6870c8")
m3=qpcr_leud_conc(qpcr_chemo_sub, "MC3")
colls=c('gray50',"#56ae6c")
m4=qpcr_leud_conc(qpcr_chemo_sub, "MC4")
colls=c('gray50',"#af953c")
m5=qpcr_leud_conc(qpcr_chemo_sub, "MC5")
colls=c('gray50',"#ba4a4f")
m6=qpcr_leud_conc(qpcr_chemo_sub, "MC6")

plot_grid(m2,m3,m4,m5,m6,aligh='hv',ncol=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_chemo_LRNA-1.png)<!-- -->

# 16S Community Data


```r
library(phyloseq)
x<-read.csv(file='/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/asv-table.csv',header=TRUE,row.names=1)
OTU = otu_table(x, taxa_are_rows=T)
taxa<-read.csv(file="/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/taxonomy.csv",header=TRUE,row.names=1)
t<-as.matrix(taxa)
tax2<-tax_table(t)
map<-import_qiime_sample_data("/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/map_exp4.txt")
tree=read.tree('/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/tree.nwk')
phyo = phyloseq(OTU, tax2,map,tree)
phyo
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2102 taxa and 96 samples ]
## sample_data() Sample Data:       [ 96 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 2102 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 2102 tips and 2100 internal nodes ]
```

```r
phyo1 = subset_taxa(phyo, !Order=="Chloroplast")
phyo2 = subset_taxa(phyo1, !Family=="Mitochondria")
phyo=phyo2
phyo
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1921 taxa and 96 samples ]
## sample_data() Sample Data:       [ 96 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 1921 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1921 tips and 1919 internal nodes ]
```


## Alpha Diversity 

### ASV Barplot


```r
phyo_abund=transform_sample_counts(phyo, function(x) (x/sum(x)))
pd <- psmelt(phyo_abund)

colors=readRDS('/Users/oliviaahern/Documents/GitHub/Exp4/colors_asvs.rds')
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Strain)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) +
  labs(x=" ", y = "ASV Relative Abundance") +
  theme_bw() +
  theme(axis.text = element_text(color ='black',size=12, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=14),
        legend.position = "none",
        panel.grid=element_blank(),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
        strip.background = element_rect(
          color="white", fill="white", size=1, linetype="solid"),
        panel.spacing = unit(0.05, "lines"))

d_rgn
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/bar1-1.png)<!-- -->



### interactive ASV Barplot > 0.001%
Under construction


```r
phyo_nop=subset_samples(phyo, B_C!='pond')

library(plyr)
# get abundance in %
phy <- transform_sample_counts(phyo_nop, function(x) x/sum(x))
# create dataframe from phyloseq object
dat <- psmelt(phy)
# convert Phylum to a character vector from a factor because R
dat$OTU <- as.character(dat$OTU)
# group dataframe by Phylum, calculate median rel. abundance
medians <- ddply(dat, ~OTU, function(x) c(median=median(x$Abundance)))
# find Phyla whose rel. abund. is less than 0.1%
remainder <- medians[medians$median <= 0.00001,]$OTU
# change their name to "Other"
dat[dat$OTU %in% remainder,]$OTU <- 'Other <0.1%'
# boxplot

blah=data.frame(tax_table(phyo_abund))
nan=cbind(row.names(blah), (colors))
nano = data.frame(nan)
colnames(nano) = c("OTU", "Colors")
subs=levels(as.factor(dat$OTU))
sub=subset(nano, OTU %in% subs)
color2=c(sub$Colors, 'gray80')


p_rgn = ggplot(dat, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = OTU)) +
  scale_fill_manual(values=as.character(t((color2)))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) +
  labs(x=" ", y = "ASV Relative Abundance") +
  theme_bw() +
  theme(axis.text = element_text(color ='black',size=12, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=14),
        legend.position = "bottom",
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
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/bar2-1.png)<!-- -->

```r
#library(plotly)
#fig <- ggplotly(p_rgn)
#fig
```


### ASV Greater than 5%

```r
x<-read.csv(file='/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/otu_tab_abund.csv',header=TRUE,row.names=1)
OTU = otu_table(x, taxa_are_rows=T)
k=c(0,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0)
c2=c('gray70','#4cc65a','#e78eb9','#4c902c','#87bb37','#914574','#886b2c','#bb6da1','#e273d2','#934993','#d4a771','#627037','#dd9f54','#e18c2b','#b344b8','#bf88da','#9662dc','#b17235','#9f4c28','#e68061','#82bf6f','#817bd6','#c5392e','#4b71e5','#4d65a5','#a3ae68','#df7f81','#48b1da','#4fc2b3','#52a27c','#a14d51','#d7408e','#e03e5d','#2b6f4d','#b6395d','#49c385','#b1b23a','#cca333','#3f8747')
legs=c('Other','Gammaproteobacteria RS62 marine group ASV 05c57','Bacteroidia Sphingobacteriales ASV 093af','Bacteroidia uncultured ASV 095f9','Gammaproteobacteria Marinomonas ASV 0be80','Alphaproteobacteria Rhizobium ASV 0d102','Gammaproteobacteria  ASV 0e8b3','Alphaproteobacteria Marivita ASV 0ed10','Verrucomicrobiae Diplosphaera ASV 1a969','Alphaproteobacteria Rhizobium ASV 2db5b','Gammaproteobacteria Methyloversatilis ASV 3cd78','Gammaproteobacteria Aeromonas ASV 3fd6c','Gammaproteobacteria Vibrio ASV 46b34','Alphaproteobacteria Hyphomonas ASV 67ca5','Bacteroidia Cryomorphaceae ASV 696c1','Gammaproteobacteria Hydrogenophaga ASV 6f94b','Gammaproteobacteria Methylotenera ASV 72ceb','Alphaproteobacteria Rhodobacteraceae ASV 7744d','Gammaproteobacteria Cellvibrio ASV 78ba2','Campylobacteria Arcobacter ASV 7b5c2','Alphaproteobacteria Rhodobacteraceae ASV 7e7d9','Alphaproteobacteria Azospirillum ASV 88fc6','Alphaproteobacteria Rhizobium ASV 923fc','Gammaproteobacteria Xanthomonadaceae ASV 96c75','Gammaproteobacteria Rhodocyclaceae ASV 9906c','Gammaproteobacteria Methylotenera ASV c4deb','Alphaproteobacteria Hyphomonas ASV c7d5a','Alphaproteobacteria  ASV c896f','Bacteroidia Cryomorphaceae ASV cae83','Alphaproteobacteria Pseudorhodobacter ASV daa79','Gammaproteobacteria Methylotenera ASV dfa40','Gammaproteobacteria Methylophilus ASV e96be','Gammaproteobacteria Zoogloea ASV efc35','Alphaproteobacteria Rhodobacteraceae ASV f19ba','Gammaproteobacteria Oceanobacter ASV f7e8d','Bacteroidia NS11-12 marine group ASV fbc54','Bacteroidia Candidatus Aquirestis ASV fbc59','Gammaproteobacteria Methylophilus ASV fc169','Campylobacteria Arcobacter ASV fe8ec')
par(mar=c(15,5,2,1))
barplot(OTU, las=2, cex.names=0.5, col=c2,space=k,
        ylab="ASV Relative Abundance")
box(which='plot')
legend(0,-12, ncol=3,
       bty='n',xpd=T, pch=22, pt.bg=c2, legend=legs,cex=0.6)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/bar_joe-1.png)<!-- -->


### interactive Genera Barplot
Working on this plot 


```r
phyo_nop_g=tax_glom(phyo_nop, taxrank="Genus")


library(plyr)
# get abundance in %
phy <- transform_sample_counts(phyo_nop_g, function(x) x/sum(x))
# create dataframe from phyloseq object
dat <- psmelt(phy)
# convert Phylum to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
dat$Genus[dat$Genus==""] <- NA

# group dataframe by Phylum, calculate median rel. abundance
medians <- ddply(dat, ~Genus, function(x) c(median=median(x$Abundance)))
# find Phyla whose rel. abund. is less than 0.1%
remainder <- medians[medians$median <= 0.00001,]$Genus
# change their name to "Other"
dat[dat$Genus %in% remainder,]$Genus <- 'Other <0.1%'
# boxplot

blah=data.frame(tax_table(phyo_abund))

nan=cbind(blah$Genus, (colors))
nano = data.frame(nan)
colnames(nano) = c("Genus", "Colors")
subs=levels(as.factor(dat$Genus))
sub=subset(nano, Genus %in% subs)
color4=c(sub$Colors, 'gray80')


p_rgn = ggplot(dat, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Genus)) +
  scale_fill_manual(values=as.character(t((color4)))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) +
  labs(x=" ", y = "Genus Relative Abundance") +
  theme_bw() +
  theme(axis.text = element_text(color ='black',size=12, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=14),
        legend.position = "bottom",
        panel.grid=element_blank(),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.25, 'cm')
        #strip.text.x = element_text(
        # size = 10, color = "black", face = "bold"),
        #strip.background = element_rect(
        # color="white", fill="white", size=1, linetype="solid"),
        #panel.spacing = unit(0.05, "lines")
  )
p_rgn
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/bar5-1.png)<!-- -->

```r
#library(plotly)
#fig <- ggplotly(p_rgn)
#fig
```



### Class Barplot


```r
micro=microbiome::aggregate_taxa(phyo, "Class")

phyo_abund=transform_sample_counts(micro, function(x) (x/sum(x)))
pd <- psmelt(phyo_abund)

colors=readRDS("/Users/oliviaahern/Documents/GitHub/Exp4/c_class.rds")
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Class)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) +
  labs(x=" ", y = "Class Relative Abundance") +
  theme_bw() +
  theme(axis.text = element_text(color ='black',size=12, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=14),
       legend.position = "bottom",
        panel.grid=element_blank(),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
       legend.text = element_text(size=6),
        legend.key.size = unit(0.25, 'cm')
        #strip.background = element_rect(
        #  color="white", fill="white", size=1, linetype="solid"),
       # panel.spacing = unit(0.05, "lines"))
)
d_rgn
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/cbar5-1.png)<!-- -->



## Beta Diversity 

### Aitchison's Distance

Euclidean distance of centered log ratio transformed data 

#### PCoA with Pond 

Transform data - remove ASVs that are less than 0.01% of total dataset and minimum reads of 2. 
Aitchinson's distance - Euclidean distance of centered log ratio data. 



```r
library(CoDaSeq)
library(compositions)
input=t(data.frame(otu_table(phyo)))
dim(input)
```

```
## [1]   96 1921
```

```r
d.subset <- codaSeq.filter(input, 
                        samples.by.row=T,min.reads=2,min.prop =0.0001)

dim(d.subset)
```

```
## [1] 1145   96
```

```r
log_rats <- (compositions::clr(t(d.subset)))
OTU=otu_table(t(log_rats),taxa_are_rows = T)
comp_clr=phyloseq(OTU,map,tax2)


otu=t(otu_table(comp_clr))
euc=vegdist(otu,'euc')
p=prcomp(euc)
# summary(p)
{
par(mar=c(7,8,1,10))
plot(p$x[,1],p$x[,2], pch = sample_data(comp_clr)$pch,
     bg=sample_data(comp_clr)$col,cex.lab=1.3,
     xlab= 'PCoA1 62.10%', ylab='PCoA2 13.71%',cex=1.5)
legend(60,60, legend=c("MC1-C12 Batch","MC1-C12 Chemostat",
                            "MC2-Met Batch", "MC2-Met Chemostat",
                            "MC3-Eth Batch", "MC3-Eth Chemostat",
                            "MC4-Ace Batch", "MC4-Ace Chemostat",
                            "MC5-Glu Batch", "MC5-Glu Chemostat",
                            "MC6-Xyl Batch", "MC6-Xyl Chemostat",
                            "Pond"),xpd=T,
       pt.bg=c("gray70","gray70","#a44f9a","#a44f9a","#6870c8","#6870c8",
             "#56ae6c","#56ae6c","#af953c","#af953c","#ba4a4f","#ba4a4f",
             "black"),
       pch=c(21,23,21,23,21,23,21,23,21,23,21,23, 24),
       bty='n',cex=1)
ordiellipse(p, groups=sample_data(comp_clr)$Timepoint)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/pcoa_all_samples-1.png)<!-- -->

#### PCoA no Pond
Transform data - remove ASVs that are less than 0.01% of total dataset and minimum reads of 2. 
Aitchinson's distance - Euclidean distance of centered log ratio data. 

Variables influencing distance matrix:

* Treatment (AKA MC) R2 = 0.07573, F = 1.7734, p =0.007
* Timepoint R2 = 0.01842, F = 2.1569, p = 0.038
* Batch or Chemostat R2 = 0.04348, F = 3.4599, p = 0.007
* POC R2 = 0.05118, F = 5.5874, p = 0.001
* Atomic %13C R2 = 0.02144, F = 2.3410, p = 0.027


```r
data1=subset_samples(phyo, Treatment !="Pond")
input=t(data.frame(otu_table(data1)))
dim(input)
```

```
## [1]   94 1921
```

```r
d.subset <- codaSeq.filter(input, 
                           samples.by.row=T,min.reads=2,min.prop =0.00001)

dim(d.subset)
```

```
## [1] 1326   94
```

```r
log_rats <- (compositions::clr(t(d.subset)))
OTU=otu_table(t(log_rats),taxa_are_rows = T)
comp_clr=phyloseq(OTU,map,tax2)


otu=t(otu_table(comp_clr))
euc=vegdist(otu,'euc')
p=prcomp(euc)
#summary(p)
{
par(mar=c(7,8,1,10))
plot(p$x[,1],p$x[,2], pch = sample_data(data1)$pch,
     bg=sample_data(data1)$col,cex.lab=1.3,
     xlab= 'PCoA1 63.61%', ylab='PCoA2 12.26%',cex=1.5,
     yaxt='n',xaxt='n')
  axis(2, at = c(-40,-20,0,20,40,60))
  axis(1,at=c(-80,-40,0,40))
legend(60,40, legend=c("MC1-C12 Batch","MC1-C12 Chemostat",
                            "MC2-Met Batch", "MC2-Met Chemostat",
                            "MC3-Eth Batch", "MC3-Eth Chemostat",
                            "MC4-Ace Batch", "MC4-Ace Chemostat",
                            "MC5-Glu Batch", "MC5-Glu Chemostat",
                            "MC6-Xyl Batch", "MC6-Xyl Chemostat"),xpd=T,
       pt.bg=c("gray70","gray70","#a44f9a","#a44f9a","#6870c8","#6870c8",
             "#56ae6c","#56ae6c","#af953c","#af953c","#ba4a4f","#ba4a4f"),
       pch=c(21,23,21,23,21,23,21,23,21,23,21,23),
       bty='n',cex=1)
ordiellipse(p, groups=sample_data(comp_clr)$Timepoint)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/pcoa_nopond-1.png)<!-- -->

```r
sub=subset_samples(comp_clr, POC !="NA")

otu=t(otu_table(sub))
euc=vegdist(otu,'euc')
set.seed(1234)
adonis2(euc~ sample_data(sub)$Treatment + sample_data(sub)$Timepoint + sample_data(sub)$B_C + sample_data(sub)$POC + sample_data(sub)$Atomic_C, by='margin', permutations = 9999)
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 9999
## 
## adonis2(formula = euc ~ sample_data(sub)$Treatment + sample_data(sub)$Timepoint + sample_data(sub)$B_C + sample_data(sub)$POC + sample_data(sub)$Atomic_C, permutations = 9999, by = "margin")
##                            Df SumOfSqs      R2      F Pr(>F)    
## sample_data(sub)$Treatment  5   2209.8 0.07555 1.6902 0.0060 ** 
## sample_data(sub)$Timepoint  1    535.7 0.01831 2.0488 0.0351 *  
## sample_data(sub)$B_C        1    831.9 0.02844 3.1816 0.0031 ** 
## sample_data(sub)$POC        1   1078.2 0.03686 4.1232 0.0005 ***
## sample_data(sub)$Atomic_C   1    717.5 0.02453 2.7441 0.0082 ** 
## Residual                   62  16212.1 0.55424                  
## Total                      71  29251.0 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#adonis2(euc~  sample_data(sub)$Timepoint + sample_data(sub)$B_C + sample_data(sub)$POC + sample_data(sub)$Atomic_C, by='margin', permutations = 9999, strata=sample_data(sub)$Treatment)
```

#### Batch PCoA

Ellipses represent 80% confidence intervals around the centroid of the Timepoint.

Data is significantly tied to the timepoint sampled (adonis2, R2=0.23, F=9.4, p<0.001) not the Treatment/MC (adonis2, R2=0.08, F=0.7, p=0.973).


```r
data1=subset_samples(comp_clr, B_C =="Batch")
otu=t(otu_table(data1))
euc=vegdist(otu,'euc')

p=prcomp(euc)
#summary(p)
{
  par(mar=c(7,8,1,10))
  plot(p$x[,1],p$x[,2], pch = sample_data(data1)$pch,
       bg=sample_data(data1)$col,cex.lab=1.3,
       xlab= 'PCoA1 54.63%', ylab='PCoA2 14.00%',cex=1.5)
  ordiellipse(p, groups=sample_data(data1)$Timepoint,
              label=FALSE)

  legend(40,5, legend=c("MC1-C12 Batch",
                            "MC2-Met Batch", 
                            "MC3-Eth Batch",
                            "MC4-Ace Batch", 
                            "MC5-Glu Batch", 
                            "MC6-Xyl Batch"),xpd=T,
       pt.bg=c("gray70","#a44f9a","#6870c8",
             "#56ae6c","#af953c","#ba4a4f"),
       pch=21,
       bty='n',cex=1)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/pcoa_batch-1.png)<!-- -->

```r
adonis2(euc~sample_data(data1)$Timepoint + sample_data(data1)$Treatment + sample_data(data1)$POC + sample_data(data1)$Atomic_C,
        by='margin', perm=999)
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = euc ~ sample_data(data1)$Timepoint + sample_data(data1)$Treatment + sample_data(data1)$POC + sample_data(data1)$Atomic_C, permutations = 999, by = "margin")
##                              Df SumOfSqs      R2      F Pr(>F)    
## sample_data(data1)$Timepoint  1   1647.7 0.12289 6.3448  0.001 ***
## sample_data(data1)$Treatment  5   1343.5 0.10020 1.0347  0.395    
## sample_data(data1)$POC        1   1490.0 0.11113 5.7374  0.001 ***
## sample_data(data1)$Atomic_C   1    183.6 0.01370 0.7071  0.718    
## Residual                     27   7011.8 0.52297                  
## Total                        35  13407.5 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```



#### Chemostat PCoA

Ellipses represent 80% confidence intervals around the centroid of the Treatment. 

Chemostat data is significantly associated with both Timepoint (adonis2, R2=0.06, F=5.5, p<0.001) and Treatment/MC (adonis2, R2=0.39, F=6.8, p<0.001).


```r
data1=subset_samples(comp_clr, B_C =="Chemostat")
otu=t(otu_table(data1))
euc=vegdist(otu,'euc')

p=prcomp(euc)
#summary(p)
{
  par(mar=c(7,8,1,10))
  plot(p$x[,1],p$x[,2], pch = sample_data(data1)$pch,
       bg=sample_data(data1)$col,cex.lab=1.3,
       xlab= 'PCoA1 34.58%', ylab='PCoA2 15.20%',cex=1.5)
  ordiellipse(p, groups=sample_data(data1)$Treatment,label=TRUE)
legend(43,5, legend=c("MC1-C12 Chemostat",
                            "MC2-Met Chemostat", 
                            "MC3-Eth Chemostat",
                            "MC4-Ace Chemostat", 
                            "MC5-Glu Chemostat", 
                            "MC6-Xyl Chemostat"),xpd=T,
       pt.bg=c("gray70","#a44f9a","#6870c8",
             "#56ae6c","#af953c","#ba4a4f"),
       pch=23,
       bty='n',cex=1)
}


adonis2(euc~sample_data(data1)$Timepoint + sample_data(data1)$Treatment + sample_data(data1)$POC + sample_data(data1)$Atomic_C,
        by='margin', perm=999)
anosim(euc, sample_data(data1)$Timepoint)
```

### Phylogenetic Distance (philr)

#### PCoA with Pond

```r
library(philr)
# filter low abudnace
library(CoDaSeq)
library(compositions)
input=t(data.frame(otu_table(phyo)))
dim(input)
```

```
## [1]   96 1921
```

```r
d.subset <- codaSeq.filter(input, 
                        samples.by.row=T,min.reads=2,min.prop =0.0001)

dim(d.subset)
```

```
## [1] 1145   96
```

```r
log_rats <- (compositions::clr(t(d.subset)))
OTU=otu_table(t(log_rats),taxa_are_rows = T)
comp_clr=phyloseq(OTU,map,tax2)

tax=data.frame(tax_table(comp_clr))
otus=tax$Strain
my_subset <- subset(otu_table(phyo), rownames(otu_table(phyo)) %in% otus)
new_physeq <- merge_phyloseq(my_subset, tax_table(phyo), sample_data(phyo), phy_tree(phyo))


GP <- transform_sample_counts(new_physeq, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "Strain_31041e0311df2ecb7962b142422de261/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
treefr <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, treefr, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')


gp.dist <- dist(gp.philr, method="euclidean")
p <- prcomp(gp.dist)




{
par(mar=c(7,8,1,10))
plot(p$x[,1],p$x[,2], pch = sample_data(phyo)$pch,
     bg=sample_data(phyo)$col,cex.lab=1.3,
     xlab= 'PCoA1 78.92%', ylab='PCoA2 9.00%',cex=1.5)
legend(80,45, legend=c("MC1-C12 Batch","MC1-C12 Chemostat",
                            "MC2-Met Batch", "MC2-Met Chemostat",
                            "MC3-Eth Batch", "MC3-Eth Chemostat",
                            "MC4-Ace Batch", "MC4-Ace Chemostat",
                            "MC5-Glu Batch", "MC5-Glu Chemostat",
                            "MC6-Xyl Batch", "MC6-Xyl Chemostat",
                            "Pond"),xpd=T,
       pt.bg=c("gray70","gray70","#a44f9a","#a44f9a","#6870c8","#6870c8",
             "#56ae6c","#56ae6c","#af953c","#af953c","#ba4a4f","#ba4a4f",
             "black"),
       pch=c(21,23,21,23,21,23,21,23,21,23,21,23, 24),
       bty='n',cex=1)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/pcoa_phil1-1.png)<!-- -->




#### PCoA no Pond

Samples significantly tied to 

* Treatment (aka MC) R2 = 0.06269, F = 2.6671, p = 0.007
* Timepoint R2 = 0.01441, F = 3.0641, p = 0.046
* Batch or Chemostat R2 = 0.02973, F = 6.3246, p = 0.005
* POC (uM) R2 = 0.01730, F = 3.6802, p = 0.027
* Atomic %13C R2 = 0.02860, F = 6.0829, p = 0.004


```r
phyo1=subset_samples(phyo, B_C !="Pond")
input=t(data.frame(otu_table(phyo1)))
dim(input)
```

```
## [1]   90 1921
```

```r
d.subset <- codaSeq.filter(input, 
                        samples.by.row=T,min.reads=2,min.prop =0.00001)

dim(d.subset)
```

```
## [1] 1272   90
```

```r
log_rats <- (compositions::clr(t(d.subset)))
OTU=otu_table(t(log_rats),taxa_are_rows = T)
comp_clr=phyloseq(OTU,map,tax2)

tax=data.frame(tax_table(comp_clr))
otus=tax$Strain
my_subset <- subset(otu_table(phyo1), rownames(otu_table(phyo1)) %in% otus)
new_physeq <- merge_phyloseq(my_subset, tax_table(phyo1), sample_data(phyo1), phy_tree(phyo1))




GP <- transform_sample_counts(new_physeq, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "Species_/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
treefr <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, treefr, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')


gp.dist <- dist(gp.philr, method="euclidean")
p <- prcomp(gp.dist)




{
par(mar=c(7,8,1,10))
plot(p$x[,1],p$x[,2], pch = sample_data(phyo1)$pch,
     bg=sample_data(phyo1)$col,cex.lab=1.3,
     xlab= 'PCoA1 77.59%', ylab='PCoA2 10.00%',cex=1.5)
legend(120,45, legend=c("MC1-C12 Batch","MC1-C12 Chemostat",
                            "MC2-Met Batch", "MC2-Met Chemostat",
                            "MC3-Eth Batch", "MC3-Eth Chemostat",
                            "MC4-Ace Batch", "MC4-Ace Chemostat",
                            "MC5-Glu Batch", "MC5-Glu Chemostat",
                            "MC6-Xyl Batch", "MC6-Xyl Chemostat"),xpd=T,
       pt.bg=c("gray70","gray70","#a44f9a","#a44f9a","#6870c8","#6870c8",
             "#56ae6c","#56ae6c","#af953c","#af953c","#ba4a4f","#ba4a4f"),
       pch=c(21,23,21,23,21,23,21,23,21,23,21,23),
       bty='n',cex=1)
ordiellipse(p, groups=sample_data(phyo1)$Timepoint)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/pcoa_phil2-1.png)<!-- -->

```r
sub=subset_samples(phyo1, POC!="NA")
GP <- transform_sample_counts(sub, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "Species_/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
treefr <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, treefr, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')
gp.dist <- dist(gp.philr, method="euclidean")

set.seed(1234)
adonis2(gp.dist~ sample_data(sub)$Treatment + sample_data(sub)$Timepoint + sample_data(sub)$B_C + sample_data(sub)$POC + sample_data(sub)$Atomic_C, by='margin', permutations = 999)
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = gp.dist ~ sample_data(sub)$Treatment + sample_data(sub)$Timepoint + sample_data(sub)$B_C + sample_data(sub)$POC + sample_data(sub)$Atomic_C, permutations = 999, by = "margin")
##                            Df SumOfSqs      R2      F Pr(>F)   
## sample_data(sub)$Treatment  5    959.0 0.06269 2.6671  0.007 **
## sample_data(sub)$Timepoint  1    220.3 0.01441 3.0641  0.046 * 
## sample_data(sub)$B_C        1    454.8 0.02973 6.3246  0.005 **
## sample_data(sub)$POC        1    264.6 0.01730 3.6802  0.027 * 
## sample_data(sub)$Atomic_C   1    437.4 0.02860 6.0829  0.004 **
## Residual                   62   4458.5 0.29147                 
## Total                      71  15296.5 1.00000                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# adonis2(formula = gp.dist ~ sample_data(sub)$Timepoint + sample_data(sub)$B_C + sample_data(sub)$POC + sample_data(sub)$Atomic_C, permutations = 999, by = "margin", strata = sample_data(sub)$Treatment)
```




#### Batch PCoA

Samples significantly tied to 

* Treatment (aka MC) R2 = 0.04936, F = 0.7705, p = 0.710
* Timepoint R2 = 0.11537, F = 9.0043, p = 0.002
* POC (uM) R2 = 0.25811, F = 20.1457, p = 0.001
* Atomic %13C R2 = 0.01423, F = 1.1105, p = 0.333



```r
phyo1=subset_samples(new_physeq, B_C =="Batch")

GP <- transform_sample_counts(phyo1, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "Species_/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
treefr <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, treefr, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')


gp.dist <- dist(gp.philr, method="euclidean")
p <- prcomp(gp.dist)




{
  par(mar=c(7,8,1,10))
  plot(p$x[,1],p$x[,2], pch = sample_data(phyo1)$pch,
       bg=sample_data(phyo1)$col,cex.lab=1.3,
       xlab= 'PCoA1 62.21%', ylab='PCoA2 19.91%',cex=1.5)
  ordiellipse(p, groups=sample_data(phyo1)$Timepoint,
              label=TRUE, cex=1.5)

  legend(90,20, legend=c("MC1-C12 Batch",
                            "MC2-Met Batch", 
                            "MC3-Eth Batch",
                            "MC4-Ace Batch", 
                            "MC5-Glu Batch", 
                            "MC6-Xyl Batch"),xpd=T,
       pt.bg=c("gray70","#a44f9a","#6870c8",
             "#56ae6c","#af953c","#ba4a4f"),
       pch=21,
       bty='n',cex=1)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/pcoa_phil3-1.png)<!-- -->

```r
set.seed(1234)
adonis2(gp.dist~ sample_data(phyo1)$Treatment + sample_data(phyo1)$Timepoint + sample_data(phyo1)$POC + sample_data(phyo1)$Atomic_C, by='margin', permutations = 999)
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = gp.dist ~ sample_data(phyo1)$Treatment + sample_data(phyo1)$Timepoint + sample_data(phyo1)$POC + sample_data(phyo1)$Atomic_C, permutations = 999, by = "margin")
##                              Df SumOfSqs      R2       F Pr(>F)    
## sample_data(phyo1)$Treatment  5   1016.8 0.05178  0.8473  0.615    
## sample_data(phyo1)$Timepoint  1   1974.2 0.10053  8.2250  0.002 ** 
## sample_data(phyo1)$POC        1   5464.1 0.27823 22.7652  0.001 ***
## sample_data(phyo1)$Atomic_C   1    289.4 0.01473  1.2056  0.282    
## Residual                     27   6480.6 0.32999                   
## Total                        35  19638.7 1.00000                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```



#### Chemostat PCoA

Samples significantly tied to 

* Timepoint R2 = 0.07368, F = 3.2076, p = 0.021
* Treatment (AKA MC) R2 = 0.12729, F = 1.1083, p = 0.362
* POC (uM) R2 = 0.01347, F = 0.5866, p = 0.603
* Atomic %13C R2 = 0.13386, F = 5.8276, p = 0.003




```r
phyo1=subset_samples(new_physeq, B_C =="Chemostat")

GP <- transform_sample_counts(phyo1, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')

otu.table <- t(otu_table(GP))
treefr <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, treefr, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')


gp.dist <- dist(gp.philr, method="euclidean")
p <- prcomp(gp.dist)

{
  par(mar=c(7,8,1,10))
  plot(p$x[,1],p$x[,2], pch = sample_data(data1)$pch,
       bg=sample_data(data1)$col,cex.lab=1.3,
       xlab= 'PCoA1 62.23%', ylab='PCoA2 17.52%',cex=1.5)
  ordiellipse(p, groups=sample_data(data1)$Treatment,label=TRUE)
legend(180,50, legend=c("MC1-C12 Chemostat",
                            "MC2-Met Chemostat", 
                            "MC3-Eth Chemostat",
                            "MC4-Ace Chemostat", 
                            "MC5-Glu Chemostat", 
                            "MC6-Xyl Chemostat"),xpd=T,
       pt.bg=c("gray70","#a44f9a","#6870c8",
             "#56ae6c","#af953c","#ba4a4f"),
       pch=23,
       bty='n',cex=1)
}


sub=subset_samples(phyo1, POC!="NA")
GP <- transform_sample_counts(sub, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')

otu.table <- t(otu_table(GP))
treefr <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, treefr, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')
gp.dist <- dist(gp.philr, method="euclidean")
set.seed(1234)
adonis2(gp.dist~ sample_data(sub)$Treatment + sample_data(sub)$Timepoint + sample_data(sub)$POC + sample_data(sub)$Atomic_C, by='margin', permutations = 999)
```
