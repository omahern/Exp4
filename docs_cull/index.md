---
title: "TC Exp3 MC2 Culled"
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
  out_dir <- 'docs_cull';
  rmarkdown::render(input_file,
 encoding=encoding,
 output_file=file.path(dirname(input_file), out_dir, 'index.html'))})
---

# Workspace Setup


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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/gas-1.png)<!-- -->

## Atomic %C


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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/atomc2-1.png)<!-- -->



# RNA Data

## Sips

### Batch 

```r
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_30Aug22.csv',header=T)
dim(read)
```

```
## [1] 343  14
```




#### Ribogreen Raw

##### MC1-C12

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo1a-1.png)<!-- -->

##### MC1-C12 vs. MC2-Met




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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo2a-1.png)<!-- -->



##### MC1-C12 vs. MC3-Eth

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo3a-1.png)<!-- -->

##### MC1-C12 vs. MC4-Ace

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo4a-1.png)<!-- -->

##### MC1-C12 vs. MC5-Glu

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo5a-1.png)<!-- -->

##### MC1-C12 vs. MC6-Xyl

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo6a-1.png)<!-- -->

#### Ribo - Buckley 

##### MC1- C12


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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo1b-1.png)<!-- -->
##### MC1-C12- vs. MC2-Met

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo2b-1.png)<!-- -->



##### MC1-C12 vs. MC3-Eth

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo3b-1.png)<!-- -->

##### MC1-C12 vs. MC4-Ace

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo4b-1.png)<!-- -->

##### MC1-C12 vs. MC5-Glu

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo5b-1.png)<!-- -->

##### MC1-C12 vs. MC6-Xyl

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo6b-1.png)<!-- -->

#### Ribo Correct Both
##### MC1-C12 vs. MC2-Met

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo2-1.png)<!-- -->



##### MC1-C12 vs. MC3-Eth

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo3-1.png)<!-- -->

##### MC1-C12 vs. MC4-Ace

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo4-1.png)<!-- -->

##### MC1-C12 vs. MC5-Glu

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo5-1.png)<!-- -->

##### MC1-C12 vs. MC6-Xyl

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/ribo6-1.png)<!-- -->


#### qRT PCR
##### C12 T2.5, 3.5, & 4.5

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_a-1.png)<!-- -->
##### T2.5


###### MC1 vs. MC2

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
     ylim=c(0.00001,175))
lines(MC2_T2.5$Corrected, MC2_T2.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0.00001,175))
lines(MC2_T2.5$Buckley, MC2_T2.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_1-1.png)<!-- -->


###### MC1 vs. MC3


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
     ylim=c(0,230))
lines(MC3_T2.5$Corrected, MC3_T2.5$qPCR_conc,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.735,1.85),
     ylim=c(0,230))
lines(MC3_T2.5$Buckley, MC3_T2.5$qPCR_conc,type="o", pch=22, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n', cex=0.9)
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_2-1.png)<!-- -->


###### MC1 vs. MC4

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_3-1.png)<!-- -->


###### MC1 vs. MC5 - Bad

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_4-1.png)<!-- -->



###### MC1 vs. MC6 - Bad 
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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_5-1.png)<!-- -->


##### T3.5 

###### MC1 vs. MC2

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_1a-1.png)<!-- -->




###### MC1 vs. MC3


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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_2a-1.png)<!-- -->




###### MC1 vs. MC4

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_3a-1.png)<!-- -->


###### MC1 vs. MC5


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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_4a-1.png)<!-- -->



###### MC1 vs. MC6


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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_5a-1.png)<!-- -->



##### T4.5 

###### MC1 vs. MC2

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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs_cull/index_files/figure-html/qpcr_1c-1.png)<!-- -->