read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_1Sept22.csv',header=T)
read=subset(read, Buckley > 1.74)
dim(read)
MC2_T2.5=subset(read, read$MC=="MC2" & TP =="2.5")
MC2_T3.5=subset(read, read$MC=="MC2" & TP =="3.5")
MC2_T4.5=subset(read, read$MC=="MC2" & TP =="4.5")
MC2_T5.5=subset(read, read$MC=="MC2" & TP =="5.5")
MC2_T1.5=subset(read, read$MC=="MC2" & TP =="1.5")

MC1_T1.5=subset(read, read$MC=="MC1" & TP =="1.5")
MC1_T2.5=subset(read, read$MC=="MC1" & TP =="2.5")
MC1_T3.5=subset(read, read$MC=="MC1" & TP =="3.5")
MC1_T4.5=subset(read, read$MC=="MC1" & TP =="4.5")
MC1_T5.5=subset(read, read$MC=="MC1" & TP =="5.5")


read=read.csv(file="/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/SIP_Spins/atomic_C_30Aug22.csv",header=T)

MC1=subset(read1, MC=="MC1")
MC2=subset(read1, MC=="MC2")
MC5=subset(read1, MC=="MC5")
MC6=subset(read1, MC=="MC6")
MC4=subset(read1, MC=="MC4")
MC3=subset(read1, MC=="MC3")

par(mfrow=c(3,2), mar=c(5,5,1,1))

plot(MC1$TP,MC1$Balanced_Atomic_.13C, type='o',  ylim=c(0,20),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1, cex.lab=1)
lines(MC1$TP, MC2$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)

plot(MC1_T1.5$Buckley, MC1_T1.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T1.5$Buckley, MC2_T1.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T1.5", "MC2-Met-T1.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T2.5$Buckley, MC2_T2.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T3.5$Buckley, MC2_T3.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T4.5$Buckley, MC2_T4.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T5.5$Buckley, MC1_T5.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T5.5$Buckley, MC2_T5.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T5.5", "MC2-Met-T5.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')


dev.off(
)
###

MC5_T2.5=subset(read, read$MC=="MC5" & TP =="2.5")
MC6_T2.5=subset(read, read$MC=="MC6" & TP =="2.5")
MC4_T2.5=subset(read, read$MC=="MC4" & TP =="2.5")
MC3_T2.5=subset(read, read$MC=="MC3" & TP=="2.5")

par(mfrow=c(3,2), mar=c(5,5,1,1))

plot(MC1$TP,MC1$Balanced_Atomic_.13C, type='o',  ylim=c(0,30),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1, cex.lab=1)
lines(MC1$TP, MC5$Balanced_Atomic_.13C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.5)
lines(MC1$TP, MC4$Balanced_Atomic_.13C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.5)
lines(MC1$TP, MC3$Balanced_Atomic_.13C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.5)
lines(MC1$TP, MC6$Balanced_Atomic_.13C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.5)
lines(MC1$TP, MC2$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC2_T2.5$Buckley, MC2_T2.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC3_T2.5$Buckley, MC3_T2.5$qPCR_RMAX,type="o", pch=21, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n')


plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC4_T2.5$Buckley, MC4_T2.5$qPCR_RMAX,type="o", pch=21, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n')


plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC5_T2.5$Buckley, MC5_T2.5$qPCR_RMAX,type="o", pch=21, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n')


plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC6_T2.5$Buckley, MC6_T2.5$qPCR_RMAX,type="o", pch=21, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n')





###
par(mfrow=c(3,1), mar=c(5,5,1,1))

MC4_T4.5=subset(read, MC=="MC4" & TP=="4.5")
MC3_T4.5=subset(read, MC=="MC3" & TP=="4.5")

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T4.5$Buckley, MC2_T4.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.84))
lines(MC3_T4.5$Buckley, MC3_T4.5$qPCR_RMAX,type="o", pch=21, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n')


plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC4_T4.5$Buckley, MC4_T4.5$qPCR_RMAX,type="o", pch=21, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n')



### Bad spin
par(mfrow=c(3,2), mar=c(5,5,1,1))
plot(MC1$TP,MC1$Balanced_Atomic_.13C, type='o',  ylim=c(0,30),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1, cex.lab=1)
lines(MC1$TP, MC5$Balanced_Atomic_.13C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.5)
lines(MC1$TP, MC4$Balanced_Atomic_.13C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.5)
lines(MC1$TP, MC3$Balanced_Atomic_.13C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.5)
lines(MC1$TP, MC6$Balanced_Atomic_.13C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.5)
lines(MC1$TP, MC2$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.84))
lines(MC2_T3.5$Buckley, MC2_T3.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')



plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.84))
lines(MC3_T3.5$Buckley, MC3_T3.5$qPCR_RMAX,type="o", pch=21, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T3.5", "MC3-Eth-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n')



plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC4_T3.5$Buckley, MC4_T3.5$qPCR_RMAX,type="o", pch=21, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T3.5", "MC4-Ace-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n')

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC5_T3.5$Buckley, MC5_T3.5$qPCR_RMAX,type="o", pch=21, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T3.5", "MC5-Glu-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n')


plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC6_T3.5$Buckley, MC6_T3.5$qPCR_RMAX,type="o", pch=21, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T3.5", "MC6-Xyl-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n')





####






read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_1Sept22.csv',header=T)
read=subset(read, Buckley > 1.74)
dim(read)
MC2_T2.5=subset(read, read$MC=="MC2" & TP =="2.5")
MC2_T3.5=subset(read, read$MC=="MC2" & TP =="3.5")
MC2_T4.5=subset(read, read$MC=="MC2" & TP =="4.5")
MC2_T5.5=subset(read, read$MC=="MC2" & TP =="5.5")
MC2_T1.5=subset(read, read$MC=="MC2" & TP =="1.5")

MC1_T1.5=subset(read, read$MC=="MC1" & TP =="1.5")
MC1_T2.5=subset(read, read$MC=="MC1" & TP =="2.5")
MC1_T3.5=subset(read, read$MC=="MC1" & TP =="3.5")
MC1_T4.5=subset(read, read$MC=="MC1" & TP =="4.5")
MC1_T5.5=subset(read, read$MC=="MC1" & TP =="5.5")


read=read.csv(file="/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/SIP_Spins/atomic_C_30Aug22.csv",header=T)

read1=subset(read, B_C=="B")

MC1=subset(read1, MC=="MC1")
MC2=subset(read1, MC=="MC2")


par(mfrow=c(3,2), mar=c(5,5,1,1))

plot(MC1$TP,MC1$Balanced_Atomic_.13C, type='o',  ylim=c(0,20),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1, cex.lab=1)
lines(MC1$TP, MC2$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)

plot(MC1_T1.5$Buckley, MC1_T1.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T1.5$Buckley, MC2_T1.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T1.5", "MC2-Met-T1.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T2.5$Buckley, MC2_T2.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T3.5$Buckley, MC2_T3.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T4.5$Buckley, MC2_T4.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')




plot(MC1_T5.5$Buckley, MC1_T5.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T5.5$Buckley, MC2_T5.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T5.5", "MC2-Met-T5.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')




plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T3.5$Buckley, MC2_T3.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')



### all

par(mfrow=c(5,2), mar=c(3,3,1,1))

plot(MC1_T1.5$Buckley, MC1_T1.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC2_T1.5$Buckley, MC2_T1.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T1.5", "MC2-Met-T1.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC2_T2.5$Buckley, MC2_T2.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC3_T2.5$Buckley, MC3_T2.5$qPCR_RMAX,type="o", pch=21, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n')


plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC4_T2.5$Buckley, MC4_T2.5$qPCR_RMAX,type="o", pch=21, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n')


plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC5_T2.5$Buckley, MC5_T2.5$qPCR_RMAX,type="o", pch=21, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n')


plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC6_T2.5$Buckley, MC6_T2.5$qPCR_RMAX,type="o", pch=21, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n')




plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85),
     ylim=c(0,1),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC2_T4.5$Buckley, MC2_T4.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.84),
     ylim=c(0,1),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC3_T4.5$Buckley, MC3_T4.5$qPCR_RMAX,type="o", pch=21, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n')


plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC4_T4.5$Buckley, MC4_T4.5$qPCR_RMAX,type="o", pch=21, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T4.5", "MC4-Ace-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n')



plot(MC1_T5.5$Buckley, MC1_T5.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85),
     yaxt='n')
axis(2, at = c(0,0.5,1))
lines(MC2_T5.5$Buckley, MC2_T5.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T5.5", "MC2-Met-T5.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')


## qpcr conc


par(mfrow=c(3,2), mar=c(5,5,1,1))

plot(MC1$TP,MC1$Balanced_Atomic_.13C, type='o',  ylim=c(0,30),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1, cex.lab=1)
lines(MC1$TP, MC5$Balanced_Atomic_.13C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.5)
lines(MC1$TP, MC4$Balanced_Atomic_.13C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.5)
lines(MC1$TP, MC3$Balanced_Atomic_.13C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.5)
lines(MC1$TP, MC6$Balanced_Atomic_.13C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.5)
lines(MC1$TP, MC2$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83),
     ylim=c(0,18))
lines(MC2_T2.5$Buckley, MC2_T2.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83),
     ylim=c(0,23))
lines(MC3_T2.5$Buckley, MC3_T2.5$qPCR_conc,type="o", pch=21, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T2.5", "MC3-Eth-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n')


plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC4_T2.5$Buckley, MC4_T2.5$qPCR_conc,type="o", pch=21, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T2.5", "MC4-Ace-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n')


plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC5_T2.5$Buckley, MC5_T2.5$qPCR_conc,type="o", pch=21, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T2.5", "MC5-Glu-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n')


plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83),
     ylim=c(0,8))
lines(MC6_T2.5$Buckley, MC6_T2.5$qPCR_conc,type="o", pch=21, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n')




###
par(mfrow=c(3,2), mar=c(5,5,1,1))

MC4_T4.5=subset(read, MC=="MC4" & TP=="4.5")
MC3_T4.5=subset(read, MC=="MC3" & TP=="4.5")

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T4.5$Buckley, MC2_T4.5$qPCR_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.85),
     ylim=c(0,9))
lines(MC2_T4.5$Buckley, MC2_T4.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.84))
lines(MC3_T4.5$Buckley, MC3_T4.5$qPCR_RMAX,type="o", pch=21, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n')


plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.84),
     ylim=c(0,0.5))
lines(MC3_T4.5$Buckley, MC3_T4.5$qPCR_conc,type="o", pch=21, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T4.5", "MC3-Eth-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="Normalized RNA Conc.", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC4_T4.5$Buckley, MC4_T4.5$qPCR_RMAX,type="o", pch=21, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T4.5", "MC4-Ace-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n')


plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC4_T4.5$Buckley, MC4_T4.5$qPCR_conc,type="o", pch=21, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T4.5", "MC4-Ace-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n')

### 

par(mfrow=c(3,2), mar=c(5,5,1,1))

plot(MC1$TP,MC1$Balanced_Atomic_.13C, type='o',  ylim=c(0,20),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1, cex.lab=1)
lines(MC1$TP, MC2$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)

plot(MC1_T1.5$Buckley, MC1_T1.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.85),
     ylim=c(0,3))
lines(MC2_T1.5$Buckley, MC2_T1.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T1.5", "MC2-Met-T1.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T2.5$Buckley, MC1_T2.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.85),
     ylim=c(0,17))
lines(MC2_T2.5$Buckley, MC2_T2.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T2.5", "MC2-Met-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T3.5$Buckley, MC1_T3.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T3.5$Buckley, MC2_T3.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T3.5", "MC2-Met-T3.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T4.5$Buckley, MC1_T4.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.85),
     ylim=c(0,10))
lines(MC2_T4.5$Buckley, MC2_T4.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T4.5", "MC2-Met-T4.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T5.5$Buckley, MC1_T5.5$qPCR_conc,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.85))
lines(MC2_T5.5$Buckley, MC2_T5.5$qPCR_conc,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T5.5", "MC2-Met-T5.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')





### day 1
read=read.csv('/Users/oliviaahern/Documents/R/Exp4/25Aug22/qpcr_1Sept22.csv',header=T)
read=subset(read, Buckley > 1.74)
dim(read)
MC1_T1.5=subset(read, MC=="MC1" & TP=="1.5")
MC2_T1.5=subset(read, MC=="MC2" & TP=="1.5")

MC3_T1.5=subset(read, MC=="MC3" & TP=="1.5")
MC4_T1.5=subset(read, MC=="MC4" & TP=="1.5")
MC5_T1.5=subset(read, MC=="MC5" & TP=="1.5")
MC6_T1.5=subset(read, MC=="MC6" & TP=="1.5")

par(mfrow=c(3,2), mar=c(5,5,1,1))

plot(MC1$TP,MC1$Balanced_Atomic_.13C, type='o',  ylim=c(0,30),pch=21,bg='gray70',
     xlab="Time (days)", ylab="Atomic %13C", cex=1, cex.lab=1)
lines(MC1$TP, MC5$Balanced_Atomic_.13C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1.5)
lines(MC1$TP, MC4$Balanced_Atomic_.13C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1.5)
lines(MC1$TP, MC3$Balanced_Atomic_.13C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1.5)
lines(MC1$TP, MC6$Balanced_Atomic_.13C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1.5)
lines(MC1$TP, MC2$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1.5)

plot(MC1_T1.5$Buckley, MC1_T1.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83),
     ylim=c(0,1))
lines(MC2_T1.5$Buckley, MC2_T1.5$Ribo_RMAX,type="o", pch=21, bg="#a44f9a", 
      col="#a44f9a")
legend("topright", legend=c("MC1-C12-T1.5", "MC2-Met-T1.5"), pch=c(21,21),
       pt.bg=c("gray70", "#a44f9a"), bty='n')

plot(MC1_T1.5$Buckley, MC1_T1.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83),
     ylim=c(0,1))
lines(MC3_T1.5$Buckley, MC3_T1.5$Ribo_RMAX,type="o", pch=21, bg="#6870c8", 
      col="#6870c8")
legend("topright", legend=c("MC1-C12-T1.5", "MC3-Eth-T1.5"), pch=c(21,21),
       pt.bg=c("gray70", "#6870c8"), bty='n')


plot(MC1_T1.5$Buckley, MC1_T1.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC4_T1.5$Buckley, MC4_T1.5$Ribo_RMAX,type="o", pch=21, bg="#56ae6c", 
      col="#56ae6c")
legend("topright", legend=c("MC1-C12-T1.5", "MC4-Ace-T1.5"), pch=c(21,21),
       pt.bg=c("gray70", "#56ae6c"), bty='n')


plot(MC1_T1.5$Buckley, MC1_T1.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83))
lines(MC5_T1.5$Buckley, MC5_T1.5$Ribo_RMAX,type="o", pch=21, bg="#af953c", 
      col="#af953c")
legend("topright", legend=c("MC1-C12-T1.5", "MC5-Glu-T1.5"), pch=c(21,21),
       pt.bg=c("gray70", "#af953c"), bty='n')


plot(MC1_T1.5$Buckley, MC1_T1.5$Ribo_RMAX,type="o", xlab="Buoyant Density - Buckley",
     ylab="RNA Conc. (ng/uL)", bg='gray70',pch=21, xlim=c(1.76,1.83),
     ylim=c(0,1))
lines(MC6_T1.5$Buckley, MC6_T1.5$Ribo_RMAX,type="o", pch=21, bg="#ba4a4f", 
      col="#ba4a4f")
legend("topright", legend=c("MC1-C12-T2.5", "MC6-Xyl-T2.5"), pch=c(21,21),
       pt.bg=c("gray70", "#ba4a4f"), bty='n')


