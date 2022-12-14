
library(plyr)
phyo_nop=subset_samples(phyo, Treatment !="Pond")
# get abundance in %
phy <- transform_sample_counts(phyo_nop, function(x) x/sum(x)*100)
# create dataframe from phyloseq object
dat <- psmelt(phy)
# convert Phylum to a character vector from a factor because R
dat$OTU <- as.character(dat$OTU)
# group dataframe by Phylum, calculate median rel. abundance
medians <- ddply(dat, ~OTU, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 0.1%
remainder <- medians[medians$median <= 0.00001,]$OTU
# change their name to "Other"
dat[dat$OTU %in% remainder,]$OTU <- 'Other <0.1%'
# boxplot

colors=readRDS('/Users/oliviaahern/Documents/GitHub/Exp4/colors_asvs.rds')
blah=data.frame(tax_table(phy))
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


AID_R=phyo_nop

write.csv(otu_table(AID_norm),'otu_tab_abund.csv')




library(phyloseq)
x<-read.csv(file='/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/otu_tab_abund.csv',header=TRUE,row.names=1)
OTU = otu_table(x, taxa_are_rows=T)
taxa<-read.csv(file="/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/taxonomy2.csv",header=TRUE,row.names=1)
t<-as.matrix(taxa)
tax2<-tax_table(t)
map<-import_qiime_sample_data("/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/map_exp4.txt")
phyo = phyloseq(OTU, tax2,map)
phyo=phyo2
phyo

x<-read.csv(file='/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/otu_tab_abund.csv',header=TRUE,row.names=1)
OTU = otu_table(x, taxa_are_rows=T)
k=c(0,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0.75,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0)
c2=c('gray70','#4cc65a','#e78eb9','#4c902c','#87bb37','#914574','#886b2c','#bb6da1','#e273d2','#934993','#d4a771','#627037','#dd9f54','#e18c2b','#b344b8','#bf88da','#9662dc','#b17235','#9f4c28','#e68061','#82bf6f','#817bd6','#c5392e','#4b71e5','#4d65a5','#a3ae68','#df7f81','#48b1da','#4fc2b3','#52a27c','#a14d51','#d7408e','#e03e5d','#2b6f4d','#b6395d','#49c385','#b1b23a','#cca333','#3f8747')
legs=c('Other','Gammaproteobacteria RS62 marine group ASV 05c57','Bacteroidia Sphingobacteriales ASV 093af','Bacteroidia uncultured ASV 095f9','Gammaproteobacteria Marinomonas ASV 0be80','Alphaproteobacteria Rhizobium ASV 0d102','Gammaproteobacteria  ASV 0e8b3','Alphaproteobacteria Marivita ASV 0ed10','Verrucomicrobiae Diplosphaera ASV 1a969','Alphaproteobacteria Rhizobium ASV 2db5b','Gammaproteobacteria Methyloversatilis ASV 3cd78','Gammaproteobacteria Aeromonas ASV 3fd6c','Gammaproteobacteria Vibrio ASV 46b34','Alphaproteobacteria Hyphomonas ASV 67ca5','Bacteroidia Cryomorphaceae ASV 696c1','Gammaproteobacteria Hydrogenophaga ASV 6f94b','Gammaproteobacteria Methylotenera ASV 72ceb','Alphaproteobacteria Rhodobacteraceae ASV 7744d','Gammaproteobacteria Cellvibrio ASV 78ba2','Campylobacteria Arcobacter ASV 7b5c2','Alphaproteobacteria Rhodobacteraceae ASV 7e7d9','Alphaproteobacteria Azospirillum ASV 88fc6','Alphaproteobacteria Rhizobium ASV 923fc','Gammaproteobacteria Xanthomonadaceae ASV 96c75','Gammaproteobacteria Rhodocyclaceae ASV 9906c','Gammaproteobacteria Methylotenera ASV c4deb','Alphaproteobacteria Hyphomonas ASV c7d5a','Alphaproteobacteria  ASV c896f','Bacteroidia Cryomorphaceae ASV cae83','Alphaproteobacteria Pseudorhodobacter ASV daa79','Gammaproteobacteria Methylotenera ASV dfa40','Gammaproteobacteria Methylophilus ASV e96be','Gammaproteobacteria Zoogloea ASV efc35','Alphaproteobacteria Rhodobacteraceae ASV f19ba','Gammaproteobacteria Oceanobacter ASV f7e8d','Bacteroidia NS11-12 marine group ASV fbc54','Bacteroidia Candidatus Aquirestis ASV fbc59','Gammaproteobacteria Methylophilus ASV fc169','Campylobacteria Arcobacter ASV fe8ec')
par(mar=c(15,5,2,1))
barplot(OTU, las=2, cex.names=0.5, col=c2,space=k,
        ylab="ASV Relative Abundance")
box(which='plot')
legend(-10,-20, ncol=3,
       bty='n',xpd=T, pch=22, pt.bg=c2, legend=legs,cex=0.8)

