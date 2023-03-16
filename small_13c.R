plot(MC1c$TP,MC1c$Balanced_Atomic_.13C, type='o',  ylim=c(0,25),pch=21,bg='gray50',col='gray50',
     xlab="Time (days)", ylab="Atomic %13C", cex=1, cex.lab=1.3, xaxt='n')
axis(1, at=c(33.5,34.5,35.5,36.5,37.5,38.5,40.5,43.5,47.5), las=2) 
lines(MC1c$TP, MC2c$Balanced_Atomic_.13C, type="o", col="#a44f9a", bg="#a44f9a", pch=21,cex=1)
lines(MC1c$TP, MC3c$Balanced_Atomic_.13C, type="o", col="#6870c8", bg="#6870c8", pch=21,cex=1)
lines(MC1c$TP, MC4c$Balanced_Atomic_.13C, type="o", col="#56ae6c", bg="#56ae6c", pch=21,cex=1)
lines(MC1c$TP, MC5c$Balanced_Atomic_.13C, type="o", col="#af953c", bg="#af953c", pch=21,cex=1)
lines(MC1c$TP, MC6c$Balanced_Atomic_.13C, type="o", col="#ba4a4f", bg="#ba4a4f", pch=21,cex=1)
legend('topleft',legend=c("MC1-C12", "MC2-Met", "MC3-Eth", "MC4-Ace",
                          "MC5-Glu", "MC6-Xyl"), bty='n', pch =21,
       pt.bg=c("gray70", "#a44f9a", "#6870c8", "#56ae6c", "#af953c", "#ba4a4f"),ncol=2, cex=0.9)

