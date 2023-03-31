
qpcr=read.csv(file='/Users/oliviaahern/Documents/R/Exp4/Dec22/qpcr_chemostat.csv',
              header=T)
dim(qpcr)
qpcr_sub=subset(qpcr, Buckley > 1.74)
dim(qpcr_sub)



sub <- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))}
f=sub(qpcr_sub, "MC2")    
eth=sub(qpcr_sub, "MC3")    
ace=sub(qpcr_sub, "MC4") 
glu=sub(qpcr_sub, "MC5") 
xyl=sub(qpcr_sub, "MC6")



p1=<- function(df, selection){
  df_out <- df %>% 
    filter(MC %in% c("MC1",selection))
  par(mfrow=c(1,7))
  plot(subset(f,MC=="MC1" & TP=="36.5")$Buckley, subset(f, MC=="MC1" & TP=="36.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC==XX & TP=="36.5")$Buckley, subset(f, MC==XX & TP=="36.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=T, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  #box(which='plot', col='firebrick', lwd=3)
  
  
  plot(subset(f,MC=="MC1" & TP=="37.5")$Buckley, subset(f, MC=="MC1" & TP=="37.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="37.5")$Buckley, subset(f, MC=="MC2" & TP=="37.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='navy', lwd=3)
  
  
  plot(subset(f,MC=="MC1" & TP=="38.5")$Buckley, subset(f, MC=="MC1" & TP=="38.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="38.5")$Buckley, subset(f, MC=="MC2" & TP=="38.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='firebrick', lwd=3)
  
  
  plot(subset(f,MC=="MC1" & TP=="40.5")$Buckley, subset(f, MC=="MC1" & TP=="40.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="40.5")$Buckley, subset(f, MC=="MC2" & TP=="40.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  
  plot(subset(f,MC=="MC1" & TP=="43.5")$Buckley, subset(f, MC=="MC1" & TP=="43.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="43.5")$Buckley, subset(f, MC=="MC2" & TP=="43.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  
  plot(subset(f,MC=="MC1" & TP=="47.5")$Buckley, subset(f, MC=="MC1" & TP=="47.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="47.5")$Buckley, subset(f, MC=="MC2" & TP=="47.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='firebrick', lwd=3)
  
  
}





png(file='bd_data_chemo3.tiff', width=10,height=8, res=700, units='in')

par(mfrow=c(5,8), mar=c(1.9,1.9,0.25,0.25))
{plot(subset(f,MC=="MC1" & TP=="34.5")$Buckley, subset(f, MC=="MC1" & TP=="34.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="34.5")$Buckley, subset(f, MC=="MC2" & TP=="34.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=T, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "3.2",font=2)
  
  plot(subset(f,MC=="MC1" & TP=="35.5")$Buckley, subset(f, MC=="MC1" & TP=="35.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="35.5")$Buckley, subset(f, MC=="MC2" & TP=="35.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "3.8",font=2)
  
  plot(subset(f,MC=="MC1" & TP=="36.5")$Buckley, subset(f, MC=="MC1" & TP=="36.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="36.5")$Buckley, subset(f, MC=="MC2" & TP=="36.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
 axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
 axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='navy', lwd=3)
  text(1.82, 0.9, "5.11",font=2)
  
  
  plot(subset(f,MC=="MC1" & TP=="37.5")$Buckley, subset(f, MC=="MC1" & TP=="37.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="37.5")$Buckley, subset(f, MC=="MC2" & TP=="37.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='yellow', lwd=3)
  
  
  plot(subset(f,MC=="MC1" & TP=="38.5")$Buckley, subset(f, MC=="MC1" & TP=="38.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="38.5")$Buckley, subset(f, MC=="MC2" & TP=="38.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='navy', lwd=3)
  
  
  plot(subset(f,MC=="MC1" & TP=="40.5")$Buckley, subset(f, MC=="MC1" & TP=="40.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="40.5")$Buckley, subset(f, MC=="MC2" & TP=="40.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "8.5",font=2)
  
  plot(subset(f,MC=="MC1" & TP=="43.5")$Buckley, subset(f, MC=="MC1" & TP=="43.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="43.5")$Buckley, subset(f, MC=="MC2" & TP=="43.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  
  plot(subset(f,MC=="MC1" & TP=="47.5")$Buckley, subset(f, MC=="MC1" & TP=="47.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(f,MC=="MC2" & TP=="47.5")$Buckley, subset(f, MC=="MC2" & TP=="47.5")$qPCR_RMAX,
        col="#a44f9a",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='navy', lwd=3)
  text(1.82, 0.9, "13.2",font=2)
  
  
}


## MC3

{
  plot(subset(eth,MC=="MC1" & TP=="34.5")$Buckley, subset(eth, MC=="MC1" & TP=="34.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="34.5")$Buckley, subset(eth, MC=="MC3" & TP=="34.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=T, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "10.0",font=2)
  
  plot(subset(eth,MC=="MC1" & TP=="35.5")$Buckley, subset(eth, MC=="MC1" & TP=="35.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="35.5")$Buckley, subset(eth, MC=="MC3" & TP=="35.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "12.5",font=2)
  
  plot(subset(eth,MC=="MC1" & TP=="36.5")$Buckley, subset(eth, MC=="MC1" & TP=="36.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="36.5")$Buckley, subset(eth, MC=="MC3" & TP=="36.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "13.1",font=2)
  
  plot(subset(eth,MC=="MC1" & TP=="37.5")$Buckley, subset(eth, MC=="MC1" & TP=="37.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="37.5")$Buckley, subset(eth, MC=="MC3" & TP=="37.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='yellow', lwd=3)
  
  plot(subset(eth,MC=="MC1" & TP=="38.5")$Buckley, subset(eth, MC=="MC1" & TP=="38.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="38.5")$Buckley, subset(eth, MC=="MC3" & TP=="38.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  
  plot(subset(eth,MC=="MC1" & TP=="40.5")$Buckley, subset(eth, MC=="MC1" & TP=="40.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="40.5")$Buckley, subset(eth, MC=="MC3" & TP=="40.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "16.8",font=2)
  
  
  plot(subset(eth,MC=="MC1" & TP=="43.5")$Buckley, subset(eth, MC=="MC1" & TP=="43.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="43.5")$Buckley, subset(eth, MC=="MC3" & TP=="43.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  plot(subset(eth,MC=="MC1" & TP=="47.5")$Buckley, subset(eth, MC=="MC1" & TP=="47.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(eth,MC=="MC3" & TP=="47.5")$Buckley, subset(eth, MC=="MC3" & TP=="47.5")$qPCR_RMAX,
        col="#6870c8",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "14.3",font=2)
  
  
}

## ace

{plot(subset(ace,MC=="MC1" & TP=="34.5")$Buckley, subset(ace, MC=="MC1" & TP=="34.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="34.5")$Buckley, subset(ace, MC=="MC4" & TP=="34.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=T, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "9.3",font=2)
  
  plot(subset(ace,MC=="MC1" & TP=="35.5")$Buckley, subset(ace, MC=="MC1" & TP=="35.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="35.5")$Buckley, subset(ace, MC=="MC4" & TP=="35.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
 axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "11.7",font=2)
  
  plot(subset(ace,MC=="MC1" & TP=="36.5")$Buckley, subset(ace, MC=="MC1" & TP=="36.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="36.5")$Buckley, subset(ace, MC=="MC4" & TP=="36.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "12.8",font=2)
  
  plot(subset(ace,MC=="MC1" & TP=="37.5")$Buckley, subset(ace, MC=="MC1" & TP=="37.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="37.5")$Buckley, subset(ace, MC=="MC4" & TP=="37.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='firebrick', lwd=3)
  
  
  plot(subset(ace,MC=="MC1" & TP=="38.5")$Buckley, subset(ace, MC=="MC1" & TP=="38.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="38.5")$Buckley, subset(ace, MC=="MC4" & TP=="38.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  
  plot(subset(ace,MC=="MC1" & TP=="40.5")$Buckley, subset(ace, MC=="MC1" & TP=="40.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="40.5")$Buckley, subset(ace, MC=="MC4" & TP=="40.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "14.5",font=2)
  
  
  plot(subset(ace,MC=="MC1" & TP=="43.5")$Buckley, subset(ace, MC=="MC1" & TP=="43.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="43.5")$Buckley, subset(ace, MC=="MC4" & TP=="43.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  plot(subset(ace,MC=="MC1" & TP=="47.5")$Buckley, subset(ace, MC=="MC1" & TP=="47.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(ace,MC=="MC4" & TP=="47.5")$Buckley, subset(ace, MC=="MC4" & TP=="47.5")$qPCR_RMAX,
        col="#56ae6c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "15.6",font=2)
  
}


# glu

{plot(subset(glu,MC=="MC1" & TP=="34.5")$Buckley, subset(glu, MC=="MC1" & TP=="34.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(glu,MC=="MC5" & TP=="34.5")$Buckley, subset(glu, MC=="MC5" & TP=="34.5")$qPCR_RMAX,
        col="#af953c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=T, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='yellow', lwd=3)
  text(1.82, 0.9, "13.7",font=2)
  
  
  
  plot(subset(glu,MC=="MC1" & TP=="35.5")$Buckley, subset(glu, MC=="MC1" & TP=="35.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(glu,MC=="MC5" & TP=="35.5")$Buckley, subset(glu, MC=="MC5" & TP=="35.5")$qPCR_RMAX,
        col="#af953c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "16.1",font=2)
  
  plot(subset(glu,MC=="MC1" & TP=="36.5")$Buckley, subset(glu, MC=="MC1" & TP=="36.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(glu,MC=="MC5" & TP=="36.5")$Buckley, subset(glu, MC=="MC5" & TP=="36.5")$qPCR_RMAX,
        col="#af953c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "20.1",font=2)
  
  plot(subset(glu,MC=="MC1" & TP=="37.5")$Buckley, subset(glu, MC=="MC1" & TP=="37.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(glu,MC=="MC5" & TP=="37.5")$Buckley, subset(glu, MC=="MC5" & TP=="37.5")$qPCR_RMAX,
        col="#af953c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  
  plot(subset(glu,MC=="MC1" & TP=="38.5")$Buckley, subset(glu, MC=="MC1" & TP=="38.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(glu,MC=="MC5" & TP=="38.5")$Buckley, subset(glu, MC=="MC5" & TP=="38.5")$qPCR_RMAX,
        col="#af953c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  box(which='plot', col='yellow', lwd=3)
  
  
  plot(subset(glu,MC=="MC1" & TP=="40.5")$Buckley, subset(glu, MC=="MC1" & TP=="40.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(glu,MC=="MC5" & TP=="40.5")$Buckley, subset(glu, MC=="MC5" & TP=="40.5")$qPCR_RMAX,
        col="#af953c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "24.4",font=2)
  
  
  plot(subset(glu,MC=="MC1" & TP=="43.5")$Buckley, subset(glu, MC=="MC1" & TP=="43.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(glu,MC=="MC5" & TP=="43.5")$Buckley, subset(glu, MC=="MC5" & TP=="43.5")$qPCR_RMAX,
        col="#af953c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  
  plot(subset(glu,MC=="MC1" & TP=="47.5")$Buckley, subset(glu, MC=="MC1" & TP=="47.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(glu,MC=="MC5" & TP=="47.5")$Buckley, subset(glu, MC=="MC5" & TP=="47.5")$qPCR_RMAX,
        col="#af953c",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=F)
  text(1.82, 0.9, "21.2",font=2)
  
}


## xyl
{plot(subset(xyl,MC=="MC1" & TP=="34.5")$Buckley, subset(xyl, MC=="MC1" & TP=="34.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(xyl,MC=="MC6" & TP=="34.5")$Buckley, subset(xyl, MC=="MC6" & TP=="34.5")$qPCR_RMAX,
        col="#ba4a4f",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=T, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=T,cex=0.7)
  text(1.82, 0.9, "18.2",font=2)
  
  plot(subset(xyl,MC=="MC1" & TP=="35.5")$Buckley, subset(xyl, MC=="MC1" & TP=="35.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(xyl,MC=="MC6" & TP=="35.5")$Buckley, subset(xyl, MC=="MC6" & TP=="35.5")$qPCR_RMAX,
        col="#ba4a4f",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=T,cex=0.7)
  text(1.82, 0.9, "15.5",font=2)
  
  plot(subset(xyl,MC=="MC1" & TP=="36.5")$Buckley, subset(xyl, MC=="MC1" & TP=="36.5")$qPCR_RMAX,
      type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(xyl,MC=="MC6" & TP=="36.5")$Buckley, subset(xyl, MC=="MC6" & TP=="36.5")$qPCR_RMAX,
        col="#ba4a4f",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F, cex=0.7) 
  axis(1, at=c(1.76,1.8,1.83), labels=T,cex=0.7)
  box(which='plot', col='yellow', lwd=3)
  text(1.82, 0.9, "16.2",font=2)
  
  
  
  plot(subset(xyl,MC=="MC1" & TP=="37.5")$Buckley, subset(xyl, MC=="MC1" & TP=="37.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(xyl,MC=="MC6" & TP=="37.5")$Buckley, subset(xyl, MC=="MC6" & TP=="37.5")$qPCR_RMAX,
        col="#ba4a4f",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=T,cex=0.7)
  
  plot(subset(xyl,MC=="MC1" & TP=="38.5")$Buckley, subset(xyl, MC=="MC1" & TP=="38.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(xyl,MC=="MC6" & TP=="38.5")$Buckley, subset(xyl, MC=="MC6" & TP=="38.5")$qPCR_RMAX,
        col="#ba4a4f",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=T,cex=0.7)
  
  plot(subset(xyl,MC=="MC1" & TP=="40.5")$Buckley, subset(xyl, MC=="MC1" & TP=="40.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(xyl,MC=="MC6" & TP=="40.5")$Buckley, subset(xyl, MC=="MC6" & TP=="40.5")$qPCR_RMAX,
        col="#ba4a4f",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=T,cex=0.7)
  text(1.82, 0.9, "19.4",font=2)
  
  
  plot(subset(xyl,MC=="MC1" & TP=="43.5")$Buckley, subset(xyl, MC=="MC1" & TP=="43.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(xyl,MC=="MC6" & TP=="43.5")$Buckley, subset(xyl, MC=="MC6" & TP=="43.5")$qPCR_RMAX,
        col="#ba4a4f",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=T,cex=0.7)
  box(which='plot', col='navy', lwd=3)
  
  
  plot(subset(xyl,MC=="MC1" & TP=="47.5")$Buckley, subset(xyl, MC=="MC1" & TP=="47.5")$qPCR_RMAX,
       type='l', xlab="", ylab="", xlim =c(1.76,1.83),lwd=1.4,yaxt='n',xaxt='n')
  lines(subset(xyl,MC=="MC6" & TP=="47.5")$Buckley, subset(xyl, MC=="MC6" & TP=="47.5")$qPCR_RMAX,
        col="#ba4a4f",lwd=1.4)
  axis(2, at=c(0,0.5,1), labels=F) 
  axis(1, at=c(1.76,1.8,1.83), labels=T,cex=0.7)
  #box(which='plot', col='navy', lwd=3)
  text(1.82, 0.9, "17.6",font=2)
  
  
}

library(cowplot)


dev.off()
