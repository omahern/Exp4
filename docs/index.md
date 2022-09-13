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


Full experimental information, including previous experiments and timeline, can be found at [the MCP website](http://ecosystems.mbl.edu/MEP-FoodWeb/Experiments/Exp4/index.html)



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
plot(MC1$TP,MC1$mm_atomic_C, type='o',  ylim=c(0,10542.03827),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1.5, cex.lab=1.4)
lines(MC1$TP, MC2$mm_atomic_C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)
lines(MC1$TP, MC3$mm_atomic_C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.5)
lines(MC1$TP, MC4$mm_atomic_C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.5)
lines(MC1$TP, MC5$mm_atomic_C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.5)
lines(MC1$TP, MC6$mm_atomic_C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.5)
legend('topleft',legend=c("MC1-C12", "MC2-Met", "MC3-Eth", "MC4-Ace",
                          "MC5-Glu", "MC6-Xyl"), bty='n', pch =21,
       pt.bg=c("gray70", "#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f"))


plot(MC1c$TP,MC1c$mm_atomic_C, type='o',  ylim=c(0,13424.45736),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1.5, cex.lab=1.4)
lines(MC1c$TP, MC2c$mm_atomic_C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)
lines(MC1c$TP, MC3c$mm_atomic_C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.5)
lines(MC1c$TP, MC4c$mm_atomic_C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.5)
lines(MC1c$TP, MC5c$mm_atomic_C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.5)
lines(MC1c$TP, MC6c$mm_atomic_C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.5)
legend('topleft',legend=c("MC1-C12", "MC2-Met", "MC3-Eth", "MC4-Ace",
                          "MC5-Glu", "MC6-Xyl"), bty='n', pch =21,
       pt.bg=c("gray70", "#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f"))

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
  #legend('topleft',legend=c("MC1-C12", "MC2-Met", "MC3-Eth", "MC4-Ace",
  #                          "MC5-Glu", "MC6-Xyl"), bty='n', pch =21,
   #      pt.bg=c("gray70", "#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f"))

  
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
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_30Aug22.csv',header=T)
dim(read)
```

```
## [1] 343  14
```




### Ribogreen Raw

#### MC1-C12

```r
MC1_T2.5=subset(read, read$MC=="MC1" & TP =="2.5")
MC1_T3.5=subset(read, read$MC=="MC1" & TP =="3.5")
MC1_T4.5=subset(read, read$MC=="MC1" & TP =="4.5")



## Colors
# MC1 - gray 70
# MC2 #a44f9a
# MC3 #6870c8
# MC4 #56ae6c
# MC5 #af953c
# MC6 #ba4a4f

# C12-T2.5 and 3.5
par(mfrow=c(1,2))
plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
lines(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_RMAX,type="o", pch=22, bg="gray70")
lines(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
     ylab="Normalized RNA Conc.", bg='gray70',pch=23)
legend("topright", legend=c("MC1-C12-T2.5", "MC1-C12-T3.5","MC1-C12-T4.5"), pch=c(21,22,23),
       pt.bg=c("gray70", "gray70"), bty='n', cex=1)

plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
     ylim=c(0,3.75))
lines(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_conc,type="o", pch=22, bg="gray70")
lines(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
     ylab="Normalized RNA Conc.", bg='gray70',pch=23)
legend("topright", legend=c("MC1-C12-T2.5", "MC1-C12-T3.5","MC1-C12-T4.5"), pch=c(21,22,23),
       pt.bg=c("gray70", "gray70"), bty='n', cex=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo1a-1.png)<!-- -->

#### MC1-C12 vs. MC2-Met




```r
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_30Aug22.csv',header=T)
MC2_T2.5=subset(read, read$MC=="MC2" & TP =="2.5")
MC2_T3.5=subset(read, read$MC=="MC2" & TP =="3.5")
MC2_T4.5=subset(read, read$MC=="MC2" & TP =="4.5")


{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
lines(MC2_T2.5$Lueders_BD, MC2_T2.5$Ribo_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
     ylim=c(0,3.2))
lines(MC2_T2.5$Lueders_BD, MC2_T2.5$Ribo_conc,type="o",  pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')


plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
lines(MC2_T3.5$Lueders_BD, MC2_T3.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,22),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
lines(MC2_T3.5$Lueders_BD, MC2_T3.5$Ribo_conc,type="o",  pch=22, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,22),
       pt.bg=c("gray70", "#a44f9a"), bty='n')



plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
lines(MC2_T4.5$Lueders_BD, MC2_T4.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a",
      col="#a44f9a")
legend("topleft", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,22),
       pt.bg=c("gray70", "#a44f9a"), bty='n')



plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
     ylim=c(0,4.6))
lines(MC2_T4.5$Lueders_BD, MC2_T4.5$Ribo_conc,type="o", pch=22, bg="#a44f9a",
      col="#a44f9a")
legend("topleft", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,22),
       pt.bg=c("gray70", "#a44f9a"), bty='n')



plot(MC2_T2.5$Lueders_BD, MC2_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
     ylab="Normalized RNA Conc.", bg='#a44f9a',pch=21, xlim=c(1.76,1.86),
     ylim=c(0,1), col="#a44f9a")
lines(MC2_T3.5$Lueders_BD, MC2_T3.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a",
      col="#a44f9a")
lines(MC2_T4.5$Lueders_BD, MC2_T4.5$Ribo_RMAX,type="o", pch=23, bg="#a44f9a",
      col="#a44f9a")
legend("topleft", legend=c("MC2-Met-T2.5","MC2-Met-T3.5", "MC2-Met-T4.5"), pch=c(21,22,23),
       pt.bg=c("#a44f9a", "#a44f9a"), bty='n')




plot(MC2_T2.5$Lueders_BD, MC2_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
     ylab="RNA Conc. (ng/uL)", bg='#a44f9a',pch=21, xlim=c(1.76,1.86),
     ylim=c(0,5), col="#a44f9a")
lines(MC2_T3.5$Lueders_BD, MC2_T3.5$Ribo_conc,type="o", pch=22, bg="#a44f9a",
      col="#a44f9a")
lines(MC2_T4.5$Lueders_BD, MC2_T4.5$Ribo_conc,type="o", pch=23, bg="#a44f9a",
      col="#a44f9a")
legend("topleft", legend=c("MC2-Met-T2.5","MC2-Met-T3.5", "MC2-Met-T4.5"), pch=c(21,22,23),
       pt.bg=c("#a44f9a", "#a44f9a"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo2a-1.png)<!-- -->



#### MC1-C12 vs. MC3-Eth

```r
# MC1-T2.5 and MC3-Eth-T2.5
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_30Aug22.csv',header=T)

MC3_T2.5=subset(read, read$MC=="MC3" & TP =="2.5")
MC3_T3.5=subset(read, read$MC=="MC3" & TP =="3.5")
MC3_T4.5=subset(read, read$MC=="MC3" & TP =="4.5")

{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC3_T2.5$Lueders_BD, MC3_T2.5$Ribo_RMAX,type="o", pch=21, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86), ylim=c(0,2.5))
  lines(MC3_T2.5$Lueders_BD, MC3_T2.5$Ribo_conc,type="o",  pch=21, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC3_T3.5$Lueders_BD, MC3_T3.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC3_T3.5$Lueders_BD, MC3_T3.5$Ribo_conc,type="o",  pch=22, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  
  plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC3_T4.5$Lueders_BD, MC3_T4.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8",
        col="#6870c8")
  legend("topleft", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  
  plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1.2))
  lines(MC3_T4.5$Lueders_BD, MC3_T4.5$Ribo_conc,type="o", pch=22, bg="#6870c8",
        col="#6870c8")
  legend("topleft", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  
  plot(MC3_T2.5$Lueders_BD, MC3_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='#6870c8',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1), col="#6870c8")
  lines(MC3_T3.5$Lueders_BD, MC3_T3.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8",
        col="#6870c8")
  lines(MC3_T4.5$Lueders_BD, MC3_T4.5$Ribo_RMAX,type="o", pch=23, bg="#6870c8",
        col="#6870c8")
  legend("topleft", legend=c("MC3-Eth-T2.5","MC3-Eth-T3.5", "MC3-Eth-T4.5"), pch=c(21,22,23),
         pt.bg=c("#6870c8", "#6870c8"), bty='n')
  
  
  
  
  plot(MC3_T2.5$Lueders_BD, MC3_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='#6870c8',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,5), col="#6870c8")
  lines(MC3_T3.5$Lueders_BD, MC3_T3.5$Ribo_conc,type="o", pch=22, bg="#6870c8",
        col="#6870c8")
  lines(MC3_T4.5$Lueders_BD, MC3_T4.5$Ribo_conc,type="o", pch=23, bg="#6870c8",
        col="#6870c8")
  legend("topleft", legend=c("MC3-Eth-T2.5","MC3-Eth-T3.5", "MC3-Eth-T4.5"), pch=c(21,22,23),
         pt.bg=c("#6870c8", "#6870c8"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo3a-1.png)<!-- -->

#### MC1-C12 vs. MC4-Ace

```r
######
# MC1 - T2.5 vs. MC4 - Ace 
MC4_T2.5=subset(read, read$MC=="MC4" & TP =="2.5")
MC4_T3.5=subset(read, read$MC=="MC4" & TP =="3.5")
MC4_T4.5=subset(read, read$MC=="MC4" & TP =="4.5")

{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC4_T2.5$Lueders_BD, MC4_T2.5$Ribo_RMAX,type="o", pch=21, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86), ylim=c(0,2.5))
  lines(MC4_T2.5$Lueders_BD, MC4_T2.5$Ribo_conc,type="o",  pch=21, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC4_T3.5$Lueders_BD, MC4_T3.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC4_T3.5$Lueders_BD, MC4_T3.5$Ribo_conc,type="o",  pch=22, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  
  plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC4_T4.5$Lueders_BD, MC4_T4.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c",
        col="#56ae6c")
  legend("topleft", legend=c("MC1-C12-T4.5", "MC4-Ace-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  
  plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,0.7))
  lines(MC4_T4.5$Lueders_BD, MC4_T4.5$Ribo_conc,type="o", pch=22, bg="#56ae6c",
        col="#56ae6c")
  legend("topleft", legend=c("MC1-C12-T4.5", "MC4-Ace-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  
  plot(MC4_T2.5$Lueders_BD, MC4_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='#56ae6c',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1), col="#56ae6c")
  lines(MC4_T3.5$Lueders_BD, MC4_T3.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c",
        col="#56ae6c")
  lines(MC4_T4.5$Lueders_BD, MC4_T4.5$Ribo_RMAX,type="o", pch=23, bg="#56ae6c",
        col="#56ae6c")
  legend("topleft", legend=c("MC4-Ace-T2.5","MC4-Ace-T3.5", "MC4-Ace-T4.5"), pch=c(21,22,23),
         pt.bg=c("#56ae6c", "#56ae6c"), bty='n')
  
  
  
  
  plot(MC4_T2.5$Lueders_BD, MC4_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='#56ae6c',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,5), col="#56ae6c")
  lines(MC4_T3.5$Lueders_BD, MC4_T3.5$Ribo_conc,type="o", pch=22, bg="#56ae6c",
        col="#56ae6c")
  lines(MC4_T4.5$Lueders_BD, MC4_T4.5$Ribo_conc,type="o", pch=23, bg="#56ae6c",
        col="#56ae6c")
  legend("topleft", legend=c("MC4-Ace-T2.5","MC4-Ace-T3.5", "MC4-Ace-T4.5"), pch=c(21,22,23),
         pt.bg=c("#56ae6c", "#56ae6c"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo4a-1.png)<!-- -->

#### MC1-C12 vs. MC5-Glu

```r
# MC1 vs MC5 
MC5_T2.5=subset(read, read$MC=="MC5" & TP =="2.5")
MC5_T3.5=subset(read, read$MC=="MC5" & TP =="3.5")
MC5_T4.5=subset(read, read$MC=="MC5" & TP =="4.5")
{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC5_T2.5$Lueders_BD, MC5_T2.5$Ribo_RMAX,type="o", pch=21, bg="#af953c", 
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86), ylim=c(0,2.5))
  lines(MC5_T2.5$Lueders_BD, MC5_T2.5$Ribo_conc,type="o",  pch=21, bg="#af953c", 
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  
  plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC5_T3.5$Lueders_BD, MC5_T3.5$Ribo_RMAX,type="o", pch=22, bg="#af953c", 
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC5_T3.5$Lueders_BD, MC5_T3.5$Ribo_conc,type="o",  pch=22, bg="#af953c", 
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  
  
  plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC5_T4.5$Lueders_BD, MC5_T4.5$Ribo_RMAX,type="o", pch=22, bg="#af953c",
        col="#af953c")
  legend("topleft", legend=c("MC1-C12-T4.5", "MC5-Glu-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  
  
  plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1))
  lines(MC5_T4.5$Lueders_BD, MC5_T4.5$Ribo_conc,type="o", pch=22, bg="#af953c",
        col="#af953c")
  legend("topleft", legend=c("MC1-C12-T4.5", "MC5-Glu-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  
  
  plot(MC5_T2.5$Lueders_BD, MC5_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='#af953c',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1), col="#af953c")
  lines(MC5_T3.5$Lueders_BD, MC5_T3.5$Ribo_RMAX,type="o", pch=22, bg="#af953c",
        col="#af953c")
  lines(MC5_T4.5$Lueders_BD, MC5_T4.5$Ribo_RMAX,type="o", pch=23, bg="#af953c",
        col="#af953c")
  legend("topleft", legend=c("MC5-Glu-T2.5","MC5-Glu-T3.5", "MC5-Glu-T4.5"), pch=c(21,22,23),
         pt.bg=c("#af953c", "#af953c"), bty='n')
  
  
  
  
  plot(MC5_T2.5$Lueders_BD, MC5_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='#af953c',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,2.5), col="#af953c")
  lines(MC5_T3.5$Lueders_BD, MC5_T3.5$Ribo_conc,type="o", pch=22, bg="#af953c",
        col="#af953c")
  lines(MC5_T4.5$Lueders_BD, MC5_T4.5$Ribo_conc,type="o", pch=23, bg="#af953c",
        col="#af953c")
  legend("topleft", legend=c("MC5-Glu-T2.5","MC5-Glu-T3.5", "MC5-Glu-T4.5"), pch=c(21,22,23),
         pt.bg=c("#af953c", "#af953c"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo5a-1.png)<!-- -->

#### MC1-C12 vs. MC6-Xyl

```r
MC6_T2.5=subset(read, read$MC=="MC6" & TP =="2.5")
MC6_T3.5=subset(read, read$MC=="MC6" & TP =="3.5")
MC6_T4.5=subset(read, read$MC=="MC6" & TP =="4.5")
{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC6_T2.5$Lueders_BD, MC6_T2.5$Ribo_RMAX,type="o", pch=21, bg="#ba4a4f", 
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  plot(MC1_T2.5$Lueders_BD, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86), ylim=c(0,2.5))
  lines(MC6_T2.5$Lueders_BD, MC6_T2.5$Ribo_conc,type="o",  pch=21, bg="#ba4a4f", 
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  
  plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC6_T3.5$Lueders_BD, MC6_T3.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f", 
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  plot(MC1_T3.5$Lueders_BD, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC6_T3.5$Lueders_BD, MC6_T3.5$Ribo_conc,type="o",  pch=22, bg="#ba4a4f", 
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  
  
  plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC6_T4.5$Lueders_BD, MC6_T4.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
        col="#ba4a4f")
  legend("topleft", legend=c("MC1-C12-T4.5", "MC6-Xyl-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  
  
  plot(MC1_T4.5$Lueders_BD, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1.5))
  lines(MC6_T4.5$Lueders_BD, MC6_T4.5$Ribo_conc,type="o", pch=22, bg="#ba4a4f",
        col="#ba4a4f")
  legend("topleft", legend=c("MC1-C12-T4.5", "MC6-Xyl-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  
  
  plot(MC6_T2.5$Lueders_BD, MC6_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders",
       ylab="Normalized RNA Conc.", bg='#ba4a4f',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1), col="#ba4a4f")
  lines(MC6_T3.5$Lueders_BD, MC6_T3.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
        col="#ba4a4f")
  lines(MC6_T4.5$Lueders_BD, MC6_T4.5$Ribo_RMAX,type="o", pch=23, bg="#ba4a4f",
        col="#ba4a4f")
  legend("topleft", legend=c("MC6-Xyl-T2.5","MC6-Xyl-T3.5", "MC6-Xyl-T4.5"), pch=c(21,22,23),
         pt.bg=c("#ba4a4f", "#ba4a4f"), bty='n')
  
  plot(MC6_T2.5$Lueders_BD, MC6_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Lueders",
       ylab="RNA Conc. (ng/uL)", bg='#ba4a4f',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,2.5), col="#ba4a4f")
  lines(MC6_T3.5$Lueders_BD, MC6_T3.5$Ribo_conc,type="o", pch=22, bg="#ba4a4f",
        col="#ba4a4f")
  lines(MC6_T4.5$Lueders_BD, MC6_T4.5$Ribo_conc,type="o", pch=23, bg="#ba4a4f",
        col="#ba4a4f")
  legend("topleft", legend=c("MC6-Xyl-T2.5","MC6-Xyl-T3.5", "MC6-Xyl-T4.5"), pch=c(21,22,23),
         pt.bg=c("#ba4a4f", "#ba4a4f"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo6a-1.png)<!-- -->

### Ribo - Buckley 

#### MC1- C12


```r
par(mfrow=c(1,2))
plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
lines(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", pch=22, bg="gray70")
lines(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", bg='gray70',pch=23)
legend("topright", legend=c("MC1-C12-T2.5", "MC1-C12-T3.5","MC1-C12-T4.5"), pch=c(21,22,23),
       pt.bg=c("gray70", "gray70"), bty='n', cex=1)

plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
     ylim=c(0,3.75))
lines(MC1_T3.5$Buckley, MC1_T3.5$Ribo_conc,type="o", pch=22, bg="gray70")
lines(MC1_T4.5$Buckley, MC1_T4.5$Ribo_conc,type="o",  bg='gray70',pch=23)
legend("topright", legend=c("MC1-C12-T2.5", "MC1-C12-T3.5","MC1-C12-T4.5"), pch=c(21,22,23),
       pt.bg=c("gray70", "gray70"), bty='n', cex=1)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo1b-1.png)<!-- -->
#### MC1-C12- vs. MC2-Met

```r
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_30Aug22.csv',header=T)
MC2_T2.5=subset(read, read$MC=="MC2" & TP =="2.5")
MC2_T3.5=subset(read, read$MC=="MC2" & TP =="3.5")
MC2_T4.5=subset(read, read$MC=="MC2" & TP =="4.5")


{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC2_T2.5$Buckley, MC2_T2.5$Ribo_RMAX,type="o", pch=21, bg="#a44f9a", 
        col="#a44f9a")
  legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#a44f9a"), bty='n')
  
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,3.2))
  lines(MC2_T2.5$Buckley, MC2_T2.5$Ribo_conc,type="o",  pch=21, bg="#a44f9a", 
        col="#a44f9a")
  legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#a44f9a"), bty='n')
  
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC2_T3.5$Buckley, MC2_T3.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a", 
        col="#a44f9a")
  legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#a44f9a"), bty='n')
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC2_T3.5$Buckley, MC2_T3.5$Ribo_conc,type="o",  pch=22, bg="#a44f9a", 
        col="#a44f9a")
  legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#a44f9a"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC2_T4.5$Buckley, MC2_T4.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a",
        col="#a44f9a")
  legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#a44f9a"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,4.6))
  lines(MC2_T4.5$Buckley, MC2_T4.5$Ribo_conc,type="o", pch=22, bg="#a44f9a",
        col="#a44f9a")
  legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#a44f9a"), bty='n')
  
  
  
  plot(MC2_T2.5$Buckley, MC2_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='#a44f9a',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1), col="#a44f9a")
  lines(MC2_T3.5$Buckley, MC2_T3.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a",
        col="#a44f9a")
  lines(MC2_T4.5$Buckley, MC2_T4.5$Ribo_RMAX,type="o", pch=23, bg="#a44f9a",
        col="#a44f9a")
  legend("topright", legend=c("MC2-Met-T2.5","MC2-Met-T3.5", "MC2-Met-T4.5"), pch=c(21,22,23),
         pt.bg=c("#a44f9a", "#a44f9a"), bty='n')
  
  
  
  
  plot(MC2_T2.5$Buckley, MC2_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='#a44f9a',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,5), col="#a44f9a")
  lines(MC2_T3.5$Buckley, MC2_T3.5$Ribo_conc,type="o", pch=22, bg="#a44f9a",
        col="#a44f9a")
  lines(MC2_T4.5$Buckley, MC2_T4.5$Ribo_conc,type="o", pch=23, bg="#a44f9a",
        col="#a44f9a")
  legend("topright", legend=c("MC2-Met-T2.5","MC2-Met-T3.5", "MC2-Met-T4.5"), pch=c(21,22,23),
         pt.bg=c("#a44f9a", "#a44f9a"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo2b-1.png)<!-- -->



#### MC1-C12 vs. MC3-Eth

```r
# MC1-T2.5 and MC3-Eth-T2.5
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_30Aug22.csv',header=T)

MC3_T2.5=subset(read, read$MC=="MC3" & TP =="2.5")
MC3_T3.5=subset(read, read$MC=="MC3" & TP =="3.5")
MC3_T4.5=subset(read, read$MC=="MC3" & TP =="4.5")

{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC3_T2.5$Buckley, MC3_T2.5$Ribo_RMAX,type="o", pch=21, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86), ylim=c(0,2.5))
  lines(MC3_T2.5$Buckley, MC3_T2.5$Ribo_conc,type="o",  pch=21, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC3_T3.5$Buckley, MC3_T3.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC3_T3.5$Buckley, MC3_T3.5$Ribo_conc,type="o",  pch=22, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC3_T4.5$Buckley, MC3_T4.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8",
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1.2))
  lines(MC3_T4.5$Buckley, MC3_T4.5$Ribo_conc,type="o", pch=22, bg="#6870c8",
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  
  plot(MC3_T2.5$Buckley, MC3_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='#6870c8',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1), col="#6870c8")
  lines(MC3_T3.5$Buckley, MC3_T3.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8",
        col="#6870c8")
  lines(MC3_T4.5$Buckley, MC3_T4.5$Ribo_RMAX,type="o", pch=23, bg="#6870c8",
        col="#6870c8")
  legend("topright", legend=c("MC3-Eth-T2.5","MC3-Eth-T3.5", "MC3-Eth-T4.5"), pch=c(21,22,23),
         pt.bg=c("#6870c8", "#6870c8"), bty='n')
  
  
  
  
  plot(MC3_T2.5$Buckley, MC3_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='#6870c8',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,5), col="#6870c8")
  lines(MC3_T3.5$Buckley, MC3_T3.5$Ribo_conc,type="o", pch=22, bg="#6870c8",
        col="#6870c8")
  lines(MC3_T4.5$Buckley, MC3_T4.5$Ribo_conc,type="o", pch=23, bg="#6870c8",
        col="#6870c8")
  legend("topright", legend=c("MC3-Eth-T2.5","MC3-Eth-T3.5", "MC3-Eth-T4.5"), pch=c(21,22,23),
         pt.bg=c("#6870c8", "#6870c8"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo3b-1.png)<!-- -->

#### MC1-C12 vs. MC4-Ace

```r
######
# MC1 - T2.5 vs. MC4 - Ace 
MC4_T2.5=subset(read, read$MC=="MC4" & TP =="2.5")
MC4_T3.5=subset(read, read$MC=="MC4" & TP =="3.5")
MC4_T4.5=subset(read, read$MC=="MC4" & TP =="4.5")

{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC4_T2.5$Buckley, MC4_T2.5$Ribo_RMAX,type="o", pch=21, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86), ylim=c(0,2.5))
  lines(MC4_T2.5$Buckley, MC4_T2.5$Ribo_conc,type="o",  pch=21, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC4_T3.5$Buckley, MC4_T3.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC4_T3.5$Buckley, MC4_T3.5$Ribo_conc,type="o",  pch=22, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC4_T4.5$Buckley, MC4_T4.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c",
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T4.5", "MC4-Ace-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,0.7))
  lines(MC4_T4.5$Buckley, MC4_T4.5$Ribo_conc,type="o", pch=22, bg="#56ae6c",
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T4.5", "MC4-Ace-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  
  plot(MC4_T2.5$Buckley, MC4_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='#56ae6c',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1), col="#56ae6c")
  lines(MC4_T3.5$Buckley, MC4_T3.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c",
        col="#56ae6c")
  lines(MC4_T4.5$Buckley, MC4_T4.5$Ribo_RMAX,type="o", pch=23, bg="#56ae6c",
        col="#56ae6c")
  legend("topright", legend=c("MC4-Ace-T2.5","MC4-Ace-T3.5", "MC4-Ace-T4.5"), pch=c(21,22,23),
         pt.bg=c("#56ae6c", "#56ae6c"), bty='n')
  
  
  
  
  plot(MC4_T2.5$Buckley, MC4_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='#56ae6c',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,5), col="#56ae6c")
  lines(MC4_T3.5$Buckley, MC4_T3.5$Ribo_conc,type="o", pch=22, bg="#56ae6c",
        col="#56ae6c")
  lines(MC4_T4.5$Buckley, MC4_T4.5$Ribo_conc,type="o", pch=23, bg="#56ae6c",
        col="#56ae6c")
  legend("topright", legend=c("MC4-Ace-T2.5","MC4-Ace-T3.5", "MC4-Ace-T4.5"), pch=c(21,22,23),
         pt.bg=c("#56ae6c", "#56ae6c"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo4b-1.png)<!-- -->

#### MC1-C12 vs. MC5-Glu

```r
# MC1 vs MC5 
MC5_T2.5=subset(read, read$MC=="MC5" & TP =="2.5")
MC5_T3.5=subset(read, read$MC=="MC5" & TP =="3.5")
MC5_T4.5=subset(read, read$MC=="MC5" & TP =="4.5")
{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC5_T2.5$Buckley, MC5_T2.5$Ribo_RMAX,type="o", pch=21, bg="#af953c", 
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86), ylim=c(0,2.5))
  lines(MC5_T2.5$Buckley, MC5_T2.5$Ribo_conc,type="o",  pch=21, bg="#af953c", 
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC5_T3.5$Buckley, MC5_T3.5$Ribo_RMAX,type="o", pch=22, bg="#af953c", 
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC5_T3.5$Buckley, MC5_T3.5$Ribo_conc,type="o",  pch=22, bg="#af953c", 
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC5_T4.5$Buckley, MC5_T4.5$Ribo_RMAX,type="o", pch=22, bg="#af953c",
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T4.5", "MC5-Glu-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1))
  lines(MC5_T4.5$Buckley, MC5_T4.5$Ribo_conc,type="o", pch=22, bg="#af953c",
        col="#af953c")
  legend("topright", legend=c("MC1-C12-T4.5", "MC5-Glu-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#af953c"), bty='n')
  
  
  
  plot(MC5_T2.5$Buckley, MC5_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='#af953c',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1), col="#af953c")
  lines(MC5_T3.5$Buckley, MC5_T3.5$Ribo_RMAX,type="o", pch=22, bg="#af953c",
        col="#af953c")
  lines(MC5_T4.5$Buckley, MC5_T4.5$Ribo_RMAX,type="o", pch=23, bg="#af953c",
        col="#af953c")
  legend("topright", legend=c("MC5-Glu-T2.5","MC5-Glu-T3.5", "MC5-Glu-T4.5"), pch=c(21,22,23),
         pt.bg=c("#af953c", "#af953c"), bty='n')
  
  
  
  
  plot(MC5_T2.5$Buckley, MC5_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='#af953c',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,2.5), col="#af953c")
  lines(MC5_T3.5$Buckley, MC5_T3.5$Ribo_conc,type="o", pch=22, bg="#af953c",
        col="#af953c")
  lines(MC5_T4.5$Buckley, MC5_T4.5$Ribo_conc,type="o", pch=23, bg="#af953c",
        col="#af953c")
  legend("topright", legend=c("MC5-Glu-T2.5","MC5-Glu-T3.5", "MC5-Glu-T4.5"), pch=c(21,22,23),
         pt.bg=c("#af953c", "#af953c"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo5b-1.png)<!-- -->

#### MC1-C12 vs. MC6-Xyl

```r
MC6_T2.5=subset(read, read$MC=="MC6" & TP =="2.5")
MC6_T3.5=subset(read, read$MC=="MC6" & TP =="3.5")
MC6_T4.5=subset(read, read$MC=="MC6" & TP =="4.5")
{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC6_T2.5$Buckley, MC6_T2.5$Ribo_RMAX,type="o", pch=21, bg="#ba4a4f", 
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86), ylim=c(0,2.5))
  lines(MC6_T2.5$Buckley, MC6_T2.5$Ribo_conc,type="o",  pch=21, bg="#ba4a4f", 
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC6_T3.5$Buckley, MC6_T3.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f", 
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC6_T3.5$Buckley, MC6_T3.5$Ribo_conc,type="o",  pch=22, bg="#ba4a4f", 
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.86))
  lines(MC6_T4.5$Buckley, MC6_T4.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T4.5", "MC6-Xyl-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1.5))
  lines(MC6_T4.5$Buckley, MC6_T4.5$Ribo_conc,type="o", pch=22, bg="#ba4a4f",
        col="#ba4a4f")
  legend("topright", legend=c("MC1-C12-T4.5", "MC6-Xyl-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#ba4a4f"), bty='n')
  
  
  
  plot(MC6_T2.5$Buckley, MC6_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='#ba4a4f',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,1), col="#ba4a4f")
  lines(MC6_T3.5$Buckley, MC6_T3.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
        col="#ba4a4f")
  lines(MC6_T4.5$Buckley, MC6_T4.5$Ribo_RMAX,type="o", pch=23, bg="#ba4a4f",
        col="#ba4a4f")
  legend("topright", legend=c("MC6-Xyl-T2.5","MC6-Xyl-T3.5", "MC6-Xyl-T4.5"), pch=c(21,22,23),
         pt.bg=c("#ba4a4f", "#ba4a4f"), bty='n')
  
  plot(MC6_T2.5$Buckley, MC6_T2.5$Ribo_conc,type="o", xlab="Buoyant Density - Buckley",
       ylab="RNA Conc. (ng/uL)", bg='#ba4a4f',pch=21, xlim=c(1.76,1.86),
       ylim=c(0,2.5), col="#ba4a4f")
  lines(MC6_T3.5$Buckley, MC6_T3.5$Ribo_conc,type="o", pch=22, bg="#ba4a4f",
        col="#ba4a4f")
  lines(MC6_T4.5$Buckley, MC6_T4.5$Ribo_conc,type="o", pch=23, bg="#ba4a4f",
        col="#ba4a4f")
  legend("topright", legend=c("MC6-Xyl-T2.5","MC6-Xyl-T3.5", "MC6-Xyl-T4.5"), pch=c(21,22,23),
         pt.bg=c("#ba4a4f", "#ba4a4f"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo6b-1.png)<!-- -->

### Ribo Correct Both
#### MC1-C12 vs. MC2-Met

```r
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_30Aug22.csv',header=T)
MC2_T2.5=subset(read, read$MC=="MC2" & TP =="2.5")
MC2_T3.5=subset(read, read$MC=="MC2" & TP =="3.5")
MC2_T4.5=subset(read, read$MC=="MC2" & TP =="4.5")
MC1_T4.5=subset(read, read$MC=="MC1" & TP == "4.5")
MC1_T2.5=subset(read, read$MC=="MC1" & TP == "2.5")
MC1_T3.5=subset(read, read$MC=="MC1" & TP == "3.5")

{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
plot(MC1_T2.5$Corrected, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc. (ribogreen)", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC2_T2.5$Corrected, MC2_T2.5$Ribo_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC2_T2.5$Buckley, MC2_T2.5$Ribo_RMAX,type="o",  pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')


plot(MC1_T3.5$Corrected, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC2_T3.5$Corrected, MC2_T3.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,22),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC2_T3.5$Buckley, MC2_T3.5$Ribo_RMAX,type="o",  pch=22, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,22),
       pt.bg=c("gray70", "#a44f9a"), bty='n')


plot(MC1_T4.5$Corrected, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc. ", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC2_T4.5$Corrected, MC2_T4.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,22),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC2_T4.5$Buckley, MC2_T4.5$Ribo_RMAX,type="o",  pch=22, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,22),
       pt.bg=c("gray70", "#a44f9a"), bty='n')


plot(MC2_T2.5$Corrected, MC2_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='#a44f9a',pch=21, xlim=c(1.735,1.85),
     col="#a44f9a")
lines(MC2_T3.5$Corrected, MC2_T3.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a", 
      col="#a44f9a")
lines(MC2_T4.5$Corrected, MC2_T4.5$Ribo_RMAX,type="o", pch=23, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC2-Met-T2.5", "MC2-Met-T3.5", "MC2-Met-T4.5"), pch=c(21,22,23),
       pt.bg=c("#a44f9a", "#a44f9a","#a44f9a"), bty='n')

plot(MC2_T2.5$Buckley, MC2_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='#a44f9a',pch=21, xlim=c(1.735,1.85),
     col="#a44f9a")
lines(MC2_T3.5$Buckley, MC2_T3.5$Ribo_RMAX,type="o", pch=22, bg="#a44f9a", 
      col="#a44f9a")
lines(MC2_T4.5$Buckley, MC2_T4.5$Ribo_RMAX,type="o", pch=23, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC2-Met-T2.5", "MC2-Met-T3.5", "MC2-Met-T4.5"), pch=c(21,22,23),
       pt.bg=c("#a44f9a", "#a44f9a","#a44f9a"), bty='n')

}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo2-1.png)<!-- -->



#### MC1-C12 vs. MC3-Eth

```r
# MC1-T2.5 and MC3-Eth-T2.5
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_30Aug22.csv',header=T)
MC3_T2.5=subset(read, read$MC=="MC3" & TP =="2.5")
MC3_T3.5=subset(read, read$MC=="MC3" & TP =="3.5")
MC3_T4.5=subset(read, read$MC=="MC3" & TP =="4.5")
MC1_T4.5=subset(read, read$MC=="MC1" & TP == "4.5")
MC1_T2.5=subset(read, read$MC=="MC1" & TP == "2.5")
MC1_T3.5=subset(read, read$MC=="MC1" & TP == "3.5")

{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Corrected, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
  lines(MC3_T2.5$Corrected, MC3_T2.5$Ribo_RMAX,type="o", pch=21, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
  lines(MC3_T2.5$Buckley, MC3_T2.5$Ribo_RMAX,type="o",  pch=21, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  plot(MC1_T3.5$Corrected, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
  lines(MC3_T3.5$Corrected, MC3_T3.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
  lines(MC3_T3.5$Buckley, MC3_T3.5$Ribo_RMAX,type="o",  pch=22, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  plot(MC1_T4.5$Corrected, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
  lines(MC3_T4.5$Corrected, MC3_T4.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
  lines(MC3_T4.5$Buckley, MC3_T4.5$Ribo_RMAX,type="o",  pch=22, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#6870c8"), bty='n')
  
  
  plot(MC3_T2.5$Corrected, MC3_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
       ylab="Normalized RNA Conc.", bg='#6870c8',pch=21, xlim=c(1.735,1.85),
       col="#6870c8")
  lines(MC3_T3.5$Corrected, MC3_T3.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8", 
        col="#6870c8")
  lines(MC3_T4.5$Corrected, MC3_T4.5$Ribo_RMAX,type="o", pch=23, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC3-Eth-T2.5", "MC3-Eth-T3.5", "MC3-Eth-T4.5"), pch=c(21,22,23),
         pt.bg=c("#6870c8", "#6870c8","#6870c8"), bty='n')
  
  plot(MC3_T2.5$Buckley, MC3_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc. ", bg='#6870c8',pch=21, xlim=c(1.735,1.85),
       col="#6870c8")
  lines(MC3_T3.5$Buckley, MC3_T3.5$Ribo_RMAX,type="o", pch=22, bg="#6870c8", 
        col="#6870c8")
  lines(MC3_T4.5$Buckley, MC3_T4.5$Ribo_RMAX,type="o", pch=23, bg="#6870c8", 
        col="#6870c8")
  legend("topright", legend=c("MC3-Eth-T2.5", "MC3-Eth-T3.5", "MC3-Eth-T4.5"), pch=c(21,22,23),
         pt.bg=c("#6870c8", "#6870c8","#6870c8"), bty='n')
  
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo3-1.png)<!-- -->

#### MC1-C12 vs. MC4-Ace

```r
######
# MC1 - T2.5 vs. MC4 - Ace 
MC4_T2.5=subset(read, read$MC=="MC4" & TP =="2.5")
MC4_T3.5=subset(read, read$MC=="MC4" & TP =="3.5")
MC4_T4.5=subset(read, read$MC=="MC4" & TP =="4.5")


{
  par(mar=c(5,5,1,1),mfrow=c(4,2))
  plot(MC1_T2.5$Corrected, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.74,1.84))
  lines(MC4_T2.5$Corrected, MC4_T2.5$Ribo_RMAX,type="o", pch=21, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.74,1.84))
  lines(MC4_T2.5$Buckley, MC4_T2.5$Ribo_RMAX,type="o",  pch=21, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  plot(MC1_T3.5$Corrected, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.74,1.84))
  lines(MC4_T3.5$Corrected, MC4_T3.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.74,1.84))
  lines(MC4_T3.5$Buckley, MC4_T3.5$Ribo_RMAX,type="o",  pch=22, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  plot(MC1_T4.5$Corrected, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.74,1.84))
  lines(MC4_T4.5$Corrected, MC4_T4.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T4.5", "MC4-Ace-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.74,1.84))
  lines(MC4_T4.5$Buckley, MC4_T4.5$Ribo_RMAX,type="o",  pch=22, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC1-C12-T4.5", "MC4-Ace-T4.5"), pch=c(21,22),
         pt.bg=c("gray70", "#56ae6c"), bty='n')
  
  
  plot(MC4_T2.5$Corrected, MC4_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
       ylab="Normalized RNA Conc.", bg='#56ae6c',pch=21, xlim=c(1.74,1.84),
       col="#56ae6c")
  lines(MC4_T3.5$Corrected, MC4_T3.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c", 
        col="#56ae6c")
  lines(MC4_T4.5$Corrected, MC4_T4.5$Ribo_RMAX,type="o", pch=23, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC4-Ace-T2.5", "MC4-Ace-T3.5", "MC4-Ace-T4.5"), pch=c(21,22,23),
         pt.bg=c("#56ae6c", "#56ae6c","#56ae6c"), bty='n')
  
  plot(MC4_T2.5$Buckley, MC4_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
       ylab="Normalized RNA Conc.", bg='#56ae6c',pch=21, xlim=c(1.74,1.84),
       col="#56ae6c")
  lines(MC4_T3.5$Buckley, MC4_T3.5$Ribo_RMAX,type="o", pch=22, bg="#56ae6c", 
        col="#56ae6c")
  lines(MC4_T4.5$Buckley, MC4_T4.5$Ribo_RMAX,type="o", pch=23, bg="#56ae6c", 
        col="#56ae6c")
  legend("topright", legend=c("MC4-Ace-T2.5", "MC4-Ace-T3.5", "MC4-Ace-T4.5"), pch=c(21,22,23),
         pt.bg=c("#56ae6c", "#56ae6c","#56ae6c"), bty='n')
  
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo4-1.png)<!-- -->

#### MC1-C12 vs. MC5-Glu

```r
# MC1 vs MC5 
MC5_T2.5=subset(read, read$MC=="MC5" & TP =="2.5")
MC5_T3.5=subset(read, read$MC=="MC5" & TP =="3.5")
MC5_T4.5=subset(read, read$MC=="MC5" & TP =="4.5")


{
par(mar=c(5,5,1,1),mfrow=c(3,2))

plot(MC1_T2.5$Corrected, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.7, 1.83))
lines(MC5_T2.5$Corrected, MC5_T2.5$Ribo_RMAX,type="o", pch=21, bg="#af953c",
      col="#af953c")
legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.7, 1.83))
lines(MC5_T2.5$Buckley, MC5_T2.5$Ribo_RMAX,type="o", pch=21, bg="#af953c",
      col="#af953c")
legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)


plot(MC1_T3.5$Corrected, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.7, 1.83))
lines(MC5_T3.5$Corrected, MC5_T3.5$Ribo_RMAX,type="o", pch=21, bg="#af953c",
      col="#af953c")
legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.7, 1.83))
lines(MC5_T3.5$Buckley, MC5_T3.5$Ribo_RMAX,type="o", pch=21, bg="#af953c",
      col="#af953c")
legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)


plot(MC5_T2.5$Corrected, MC5_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='#af953c',pch=21, xlim=c(1.7, 1.83),
     col="#af953c")
lines(MC5_T3.5$Corrected, MC5_T3.5$Ribo_RMAX,type="o", pch=22, bg="#af953c",
      col="#af953c")
legend("topright", legend=c("MC5-Glu-T2.5", "MC5-Glu-T3.5"), pch=c(21,22),
       pt.bg=c("#af953c", "#af953c"), bty='n', cex=0.9)

plot(MC5_T2.5$Buckley, MC5_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='#af953c',pch=21, xlim=c(1.7, 1.83),
     col="#af953c")
lines(MC5_T3.5$Buckley, MC5_T3.5$Ribo_RMAX,type="o", pch=22, bg="#af953c",
      col="#af953c")
legend("topright", legend=c("MC5-Glu-T2.5", "MC5-Glu-T3.5"), pch=c(21,22),
       pt.bg=c("#af953c", "#af953c"), bty='n', cex=0.9)

}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo5-1.png)<!-- -->

#### MC1-C12 vs. MC6-Xyl

```r
MC6_T2.5=subset(read, read$MC=="MC6" & TP =="2.5")
MC6_T3.5=subset(read, read$MC=="MC6" & TP =="3.5")
MC6_T4.5=subset(read, read$MC=="MC6" & TP =="4.5")

{
par(mar=c(5,5,1,1),mfrow=c(3,2))
plot(MC1_T2.5$Corrected, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.7, 1.83))
lines(MC6_T2.5$Corrected, MC6_T2.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.7, 1.83))
lines(MC6_T2.5$Buckley, MC6_T2.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)


plot(MC1_T3.5$Corrected, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.7, 1.84))
lines(MC6_T3.5$Corrected, MC6_T3.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,22),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.75, 1.85))
lines(MC6_T3.5$Buckley, MC6_T3.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,22),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)

plot(MC6_T2.5$Corrected, MC6_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='#ba4a4f',pch=21, xlim=c(1.7, 1.84),
     col="#ba4a4f")
lines(MC6_T3.5$Corrected, MC6_T3.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
      col="#ba4a4f")
legend("topright", legend=c("MC6-Xyl-T2.5", "MC6-Xyl-T3.5"), pch=c(21,22),
       pt.bg=c("#ba4a4f", "#ba4a4f"), bty='n', cex=0.9)

plot(MC6_T2.5$Buckley, MC6_T2.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='#ba4a4f',pch=21, xlim=c(1.75, 1.85), 
     col="#ba4a4f")
lines(MC6_T3.5$Buckley, MC6_T3.5$Ribo_RMAX,type="o", pch=22, bg="#ba4a4f",
      col="#ba4a4f")
legend("topright", legend=c("MC6-Xyl-T2.5", "MC6-Xyl-T3.5"), pch=c(21,22),
       pt.bg=c("#ba4a4f", "#ba4a4f"), bty='n', cex=0.9)


}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribo6-1.png)<!-- -->

#### Fig for talk 

```r
{
par(mar=c(5,5,1,1))
plot(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o", xlab="Buoyant Density (g/mL)",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.84), ylim=c(0,1),
     yaxt='n',xaxt='n', cex.lab=1.2, bty='n',cex=1.4)
rect(xleft=1.787, xright=1.794, ybottom=0, ytop=1.03, col= rgb(0.9960784,0.9411765,0.2745098,alpha=0.8),lty=1, lwd=0)
#rect(xleft=1.807990115-0.0035, xright=1.807990115+0.0035, ybottom=0, ytop=1.03, col= rgb(0.5529412,0.08627451,0.08627451,alpha=0.7),lty=1, lwd=0)
rect(xleft=1.807990115-0.0035, xright=1.807990115+0.0035, ybottom=0, ytop=1.03, col= rgb(0.9058824,0.1803922,0.1803922,alpha=0.3),lty=1, lwd=0)

lines(MC1_T4.5$Buckley, MC1_T4.5$Ribo_RMAX,type="o",  pch=21, bg="gray70",cex=1.4)
axis(1,at=c(1.76,1.78,1.80,1.82,1.84), lwd=1)
axis(2,at=c(0,0.25,0.5,0.75,1),lwd=1)
lines(MC2_T4.5$Buckley, MC2_T4.5$Ribo_RMAX,type="o",  pch=21, bg="#a44f9a", 
      col="black",cex=1.4)
lines(MC2_T4.5$Buckley, MC2_T4.5$Ribo_RMAX,type="o",  pch=21, bg="#a44f9a", 
      col="black",cex=1.4)
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/ribome-1.png)<!-- -->

### qRT PCR
#### C12 T2.5, 3.5, & 4.5

```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected", xlim=c(1.7, 1.83),
     ylab="Normalized RNA Conc", bg='gray70',pch=21)
lines(MC1_T3.5$Corrected, MC1_T3.5$qPCR_RMAX,type="o", pch=22, bg="gray70", 
      col="black")
lines(MC1_T4.5$Corrected, MC1_T4.5$qPCR_RMAX,type="o", pch=23, bg="gray70", 
      col="black")
legend("topright", legend=c("MC1-C12-T2.5", "MC1-C12-T3.5", "MC1-C12-T4.5"), pch=c(21,22,23),
       pt.bg="gray70", bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley", xlim=c(1.73, 1.85),
     ylab="Normalized RNA Conc", bg='gray70',pch=21)
lines(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", pch=22, bg="gray70", 
      col="black")
lines(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", pch=23, bg="gray70", 
      col="black")
legend("topright", legend=c("MC1-C12-T2.5", "MC1-C12-T3.5","MC1-C12-T3.5"), pch=c(21,22,23),
       pt.bg="gray70", bty='n')

plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected", xlim=c(1.7, 1.83),ylim=c(0,20),
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21)
lines(MC1_T3.5$Corrected, MC1_T3.5$qPCR_conc,type="o", pch=22, bg="gray70", 
      col="black")
lines(MC1_T4.5$Corrected, MC1_T4.5$qPCR_conc,type="o", pch=23, bg="gray70", 
      col="black")
legend("topright", legend=c("MC1-C12-T2.5", "MC1-C12-T3.5", "MC1-C12-T4.5"), pch=c(21,22,23),
       pt.bg="gray70", bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley", xlim=c(1.73, 1.85),
     ylim=c(0,20),
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21)
lines(MC1_T3.5$Buckley, MC1_T3.5$qPCR_conc,type="o", pch=22, bg="gray70", 
      col="black")
lines(MC1_T4.5$Buckley, MC1_T4.5$qPCR_conc,type="o", pch=23, bg="gray70", 
      col="black")
legend("topright", legend=c("MC1-C12-T2.5", "MC1-C12-T3.5","MC1-C12-T3.5"), pch=c(21,22,23),
       pt.bg="gray70", bty='n')

}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_a-1.png)<!-- -->
#### T2.5


##### MC1-C12 vs. MC2-Met

```r
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_30Aug22.csv',header=T)

MC2_T2.5=subset(read, read$MC=="MC2" & TP =="2.5")
MC2_T3.5=subset(read, read$MC=="MC2" & TP =="3.5")
MC2_T4.5=subset(read, read$MC=="MC2" & TP =="4.5")
MC1_T4.5=subset(read, read$MC=="MC1" & TP == "4.5")
MC1_T2.5=subset(read, read$MC=="MC1" & TP == "2.5")
MC1_T3.5=subset(read, read$MC=="MC1" & TP == "3.5")
MC3_T2.5=subset(read, read$MC=="MC3" & TP =="2.5")
MC3_T3.5=subset(read, read$MC=="MC3" & TP =="3.5")
MC3_T4.5=subset(read, read$MC=="MC3" & TP =="4.5")
MC4_T2.5=subset(read, read$MC=="MC4" & TP =="2.5")
MC4_T3.5=subset(read, read$MC=="MC4" & TP =="3.5")
MC4_T4.5=subset(read, read$MC=="MC4" & TP =="4.5")
MC5_T2.5=subset(read, read$MC=="MC5" & TP =="2.5")
MC5_T3.5=subset(read, read$MC=="MC5" & TP =="3.5")
MC5_T4.5=subset(read, read$MC=="MC5" & TP =="4.5")
MC6_T2.5=subset(read, read$MC=="MC6" & TP =="2.5")
MC6_T3.5=subset(read, read$MC=="MC6" & TP =="3.5")
MC6_T4.5=subset(read, read$MC=="MC6" & TP =="4.5")



{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC2_T2.5$Corrected, MC2_T2.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC2_T2.5$Buckley, MC2_T2.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0.00001,20))
lines(MC2_T2.5$Corrected, MC2_T2.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0.00001,20))
lines(MC2_T2.5$Buckley, MC2_T2.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_1-1.png)<!-- -->


##### MC1-C12 vs. MC3-Eth


```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC3_T2.5$Corrected, MC3_T2.5$qPCR_RMAX,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC3_T2.5$Buckley, MC3_T2.5$qPCR_RMAX,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)


plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0,25))
lines(MC3_T2.5$Corrected, MC3_T2.5$qPCR_conc,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0,25))
lines(MC3_T2.5$Buckley, MC3_T2.5$qPCR_conc,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_2-1.png)<!-- -->


##### MC1-C12 vs. MC4-Ace

```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC4_T2.5$Corrected, MC4_T2.5$qPCR_RMAX,type="o", pch=22, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC4_T2.5$Buckley, MC4_T2.5$qPCR_RMAX,type="o", pch=22, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n', cex=0.9)

plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC4_T2.5$Corrected, MC4_T2.5$qPCR_conc,type="o", pch=22, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC4_T2.5$Buckley, MC4_T2.5$qPCR_conc,type="o", pch=22, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n', cex=0.9)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_3-1.png)<!-- -->


##### MC1 vs. MC5 - Redo Spin

REDO SPIN


```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC5_T2.5$Corrected, MC5_T2.5$qPCR_RMAX,type="o", pch=22, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC5_T2.5$Buckley, MC5_T2.5$qPCR_RMAX,type="o", pch=22, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)


plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0,7.25))
lines(MC5_T2.5$Corrected, MC5_T2.5$qPCR_conc,type="o", pch=22, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0,7.25))
lines(MC5_T2.5$Buckley, MC5_T2.5$qPCR_conc,type="o", pch=22, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_4-1.png)<!-- -->



##### MC1 vs. MC6 - Redo Spin 
REDO SPIN

```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC6_T2.5$Corrected, MC6_T2.5$qPCR_RMAX,type="o", pch=22, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85))
lines(MC6_T2.5$Buckley, MC6_T2.5$qPCR_RMAX,type="o", pch=22, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)

plot(MC1_T2.5$Corrected, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0, 40))
lines(MC6_T2.5$Corrected, MC6_T2.5$qPCR_conc,type="o", pch=22, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0, 40))
lines(MC6_T2.5$Buckley, MC6_T2.5$qPCR_conc,type="o", pch=22, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_5-1.png)<!-- -->


#### T3.5 - redo C12 qRT PCR

##### MC1-C12 vs. MC2-Met -Redo qRT PCR


```r
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_1Sept22.csv',header=T)

MC2_T2.5=subset(read, read$MC=="MC2" & TP =="2.5")
MC2_T3.5=subset(read, read$MC=="MC2" & TP =="3.5")
MC2_T4.5=subset(read, read$MC=="MC2" & TP =="4.5")
MC1_T4.5=subset(read, read$MC=="MC1" & TP == "4.5")
MC1_T2.5=subset(read, read$MC=="MC1" & TP == "2.5")
MC1_T3.5=subset(read, read$MC=="MC1" & TP == "3.5")
MC3_T2.5=subset(read, read$MC=="MC3" & TP =="2.5")
MC3_T3.5=subset(read, read$MC=="MC3" & TP =="3.5")
MC3_T4.5=subset(read, read$MC=="MC3" & TP =="4.5")
MC4_T2.5=subset(read, read$MC=="MC4" & TP =="2.5")
MC4_T3.5=subset(read, read$MC=="MC4" & TP =="3.5")
MC4_T4.5=subset(read, read$MC=="MC4" & TP =="4.5")
MC5_T2.5=subset(read, read$MC=="MC5" & TP =="2.5")
MC5_T3.5=subset(read, read$MC=="MC5" & TP =="3.5")
MC5_T4.5=subset(read, read$MC=="MC5" & TP =="4.5")
MC6_T2.5=subset(read, read$MC=="MC6" & TP =="2.5")
MC6_T3.5=subset(read, read$MC=="MC6" & TP =="3.5")
MC6_T4.5=subset(read, read$MC=="MC6" & TP =="4.5")



{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc", bg='gray70',pch=21)
lines(MC2_T3.5$Corrected, MC2_T3.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc", bg='gray70',pch=21)
lines(MC2_T3.5$Buckley, MC2_T3.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, 
     ylim=c(0.00001,20))
lines(MC2_T3.5$Corrected, MC2_T3.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, 
     ylim=c(0.00001,20))
lines(MC2_T3.5$Buckley, MC2_T3.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_1a-1.png)<!-- -->




##### MC1 vs. MC3 - Redo Spin


```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)
lines(MC3_T3.5$Corrected, MC3_T3.5$qPCR_RMAX,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)
lines(MC3_T3.5$Buckley, MC3_T3.5$qPCR_RMAX,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)


plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21,
     ylim=c(0,20))
lines(MC3_T3.5$Corrected, MC3_T3.5$qPCR_conc,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, 
     ylim=c(0,20))
lines(MC3_T3.5$Buckley, MC3_T3.5$qPCR_conc,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_2a-1.png)<!-- -->




##### MC1-C12 vs. MC4-Ace

```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)
lines(MC4_T3.5$Corrected, MC4_T3.5$qPCR_RMAX,type="o", pch=22, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)
lines(MC4_T3.5$Buckley, MC4_T3.5$qPCR_RMAX,type="o", pch=22, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n', cex=0.9)

plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21)
lines(MC4_T3.5$Corrected, MC4_T3.5$qPCR_conc,type="o", pch=22, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21)
lines(MC4_T3.5$Buckley, MC4_T3.5$qPCR_conc,type="o", pch=22, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n', cex=0.9)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_3a-1.png)<!-- -->


##### MC1-C12 vs. MC5-Glu - Redo qRT PCR


```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)
lines(MC5_T3.5$Corrected, MC5_T3.5$qPCR_RMAX,type="o", pch=22, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)
lines(MC5_T3.5$Buckley, MC5_T3.5$qPCR_RMAX,type="o", pch=22, bg="#af953c", col="#af953c")
legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)

plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21)
lines(MC5_T3.5$Corrected, MC5_T3.5$qPCR_conc,type="o", pch=22, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21)
lines(MC5_T3.5$Buckley, MC5_T3.5$qPCR_conc,type="o", pch=22, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n', cex=0.9)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_4a-1.png)<!-- -->



##### MC1-C12 vs. MC6-Xyl - redo qRT PCR


```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)
lines(MC6_T3.5$Corrected, MC6_T3.5$qPCR_RMAX,type="o", pch=22, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)
lines(MC6_T3.5$Buckley, MC6_T3.5$qPCR_RMAX,type="o", pch=22, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)

plot(MC1_T3.5$Corrected, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)

lines(MC6_T3.5$Corrected, MC6_T3.5$qPCR_conc,type="o", pch=22, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21)
lines(MC6_T3.5$Buckley, MC6_T3.5$qPCR_conc,type="o", pch=22, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n', cex=0.9)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_5a-1.png)<!-- -->



#### T4.5 

##### MC1-C12 vs. MC2-Met - redo qRT PCR


```r
{
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(MC1_T4.5$Corrected, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="Normalized RNA Conc", bg='gray70',pch=21)
lines(MC2_T4.5$Corrected, MC2_T4.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc", bg='gray70',pch=21)
lines(MC2_T4.5$Buckley, MC2_T4.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Corrected, MC1_T4.5$qPCR_conc,type="o", xlab="Buoyant Density - Lueders Corrected",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, 
     ylim=c(0.00001,10))
lines(MC2_T4.5$Corrected, MC2_T4.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21,
     ylim=c(0.00001,10))
lines(MC2_T4.5$Buckley, MC2_T4.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/qpcr_1c-1.png)<!-- -->

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
                           samples.by.row=T,min.reads=2,min.prop =0.0001)

dim(d.subset)
```

```
## [1] 748  94
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
adonis2(euc~sample_data(data1)$Timepoint + sample_data(data1)$Treatment,
        by='margin')
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = euc ~ sample_data(data1)$Timepoint + sample_data(data1)$Treatment, by = "margin")
##                              Df SumOfSqs      R2      F Pr(>F)    
## sample_data(data1)$Timepoint  1   2581.6 0.22551 9.3967  0.001 ***
## sample_data(data1)$Treatment  5    899.0 0.07853 0.6544  0.977    
## Residual                     29   7967.4 0.69596                  
## Total                        35  11448.0 1.00000                  
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
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/pcoa_chemo-1.png)<!-- -->

```r
adonis2(euc~sample_data(data1)$Timepoint + sample_data(data1)$Treatment,
        by='margin')
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = euc ~ sample_data(data1)$Timepoint + sample_data(data1)$Treatment, by = "margin")
##                              Df SumOfSqs      R2      F Pr(>F)    
## sample_data(data1)$Timepoint  1    679.1 0.06347 5.5055  0.001 ***
## sample_data(data1)$Treatment  5   4222.4 0.39465 6.8461  0.001 ***
## Residual                     47   5797.5 0.54187                  
## Total                        53  10699.0 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Phylogenetic Distance (philr)

#### PCoA with Pond

```r
library(philr)

GP <- transform_sample_counts(phyo, function(x) x+1)
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
plot(p$x[,1],p$x[,2], pch = sample_data(phyo)$pch,
     bg=sample_data(phyo)$col,cex.lab=1.3,
     xlab= 'PCoA1 79.27%', ylab='PCoA2 9.19%',cex=1.5)
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

* Treatment (aka MC) R2 = 0.16, F = 6.26, p < 0.001
* Timepoint R2 = 0.021, F = 5.33, p = 0.007
* Batch or Chemostat R2 = 0.05, F=12.80, p < 0.001


```r
phyo1=subset_samples(phyo, B_C !="Pond")

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
     xlab= 'PCoA1 78.84%', ylab='PCoA2 9.93%',cex=1.5)
legend(80,45, legend=c("MC1-C12 Batch","MC1-C12 Chemostat",
                            "MC2-Met Batch", "MC2-Met Chemostat",
                            "MC3-Eth Batch", "MC3-Eth Chemostat",
                            "MC4-Ace Batch", "MC4-Ace Chemostat",
                            "MC5-Glu Batch", "MC5-Glu Chemostat",
                            "MC6-Xyl Batch", "MC6-Xyl Chemostat"),xpd=T,
       pt.bg=c("gray70","gray70","#a44f9a","#a44f9a","#6870c8","#6870c8",
             "#56ae6c","#56ae6c","#af953c","#af953c","#ba4a4f","#ba4a4f"),
       pch=c(21,23,21,23,21,23,21,23,21,23,21,23),
       bty='n',cex=1)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/pcoa_phil2-1.png)<!-- -->

```r
adonis2(gp.dist~ sample_data(phyo1)$Treatment + sample_data(phyo1)$Timepoint + sample_data(phyo1)$B_C, by='margin')
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = gp.dist ~ sample_data(phyo1)$Treatment + sample_data(phyo1)$Timepoint + sample_data(phyo1)$B_C, by = "margin")
##                              Df SumOfSqs      R2       F Pr(>F)    
## sample_data(phyo1)$Treatment  5   2060.8 0.12580  6.2597  0.001 ***
## sample_data(phyo1)$Timepoint  1    351.2 0.02144  5.3335  0.005 ** 
## sample_data(phyo1)$B_C        1    842.8 0.05145 12.7998  0.001 ***
## Residual                     82   5399.1 0.32959                   
## Total                        89  16381.4 1.00000                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```




#### Batch PCoA

Samples significantly tied to 

* Treatment (aka MC) R2 = 0.04, F = 0.3267, p = 0.99
* Timepoint R2 = 0.21501, F = 8.3907, p = 0.002



```r
phyo1=subset_samples(phyo, B_C =="Batch")

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
       xlab= 'PCoA1 64.96%', ylab='PCoA2 18.17%',cex=1.5)
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
adonis2(gp.dist~ sample_data(phyo1)$Treatment + sample_data(phyo1)$Timepoint, by='margin')
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = gp.dist ~ sample_data(phyo1)$Treatment + sample_data(phyo1)$Timepoint, by = "margin")
##                              Df SumOfSqs      R2      F Pr(>F)    
## sample_data(phyo1)$Treatment  5    818.1 0.04186 0.3267  0.994    
## sample_data(phyo1)$Timepoint  1   4202.7 0.21501 8.3907  0.001 ***
## Residual                     29  14525.5 0.74313                  
## Total                        35  19546.3 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```



#### Chemostat PCoA

Samples significantly tied to 

* Treatment (aka MC) R2 = 0.58701, F = 16.814, p = 0.001
* Timepoint R2 = 0.08482, F = 12.148, p = 0.001



```r
phyo1=subset_samples(phyo, B_C =="Chemostat")

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
  plot(p$x[,1],p$x[,2], pch = sample_data(data1)$pch,
       bg=sample_data(data1)$col,cex.lab=1.3,
       xlab= 'PCoA1 53.93%', ylab='PCoA2 24.54%',cex=1.5)
  ordiellipse(p, groups=sample_data(data1)$Treatment,label=TRUE)
legend(85,50, legend=c("MC1-C12 Chemostat",
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
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/pcoa_phil4-1.png)<!-- -->

```r
adonis2(gp.dist~ sample_data(phyo1)$Treatment + sample_data(phyo1)$Timepoint, by='margin')
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = gp.dist ~ sample_data(phyo1)$Treatment + sample_data(phyo1)$Timepoint, by = "margin")
##                              Df SumOfSqs      R2      F Pr(>F)    
## sample_data(phyo1)$Treatment  5  12134.8 0.58701 16.814  0.001 ***
## sample_data(phyo1)$Timepoint  1   1753.5 0.08482 12.148  0.001 ***
## Residual                     47   6784.0 0.32817                  
## Total                        53  20672.3 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
