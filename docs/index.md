---
title: "Test Plate"
output: 
  html_document:
    keep_md: true
    toc: true
    toc_float: true
    toc_depth: 6
    code_folding: hide
    number_sections: false
    theme: lumen

knit: (function(input_file, encoding) {
  out_dir <- 'docs';
  rmarkdown::render(input_file,
 encoding=encoding,
 output_file=file.path(dirname(input_file), out_dir, 'index.html'))})
---


Data from the test plate from IMR using "Universal Primers" from Parada et al, 2016. I followed the pipeline from [Jesse McNichols](https://github.com/jcmcnch/eASV-pipeline-for-515Y-926R) and the [Fuhrman lab protocols page](https://www.protocols.io/private/C0F9404AB3DAEC96683F142351CEF59C?step=1).

Questions: 

+ How many reads do we get per sample?
+ Is there a relationship between cDNA concentration and No. of Reads?
+ How many reads are prokaryotic and eukaryotic? 
+ What do the eukaryotic communities look like?
+ What do my T0 communities look like compared to other Siders Pond data? 



# Raw Reads

```r
read=read.csv(file='/Users/oliviaahern/Documents/GitHub/TestPlate/raw_seqs.csv',header=T,row.names=1)
data=data.frame(read[,7:8])
#barplot(data$Raw)
```

## No. Reads Community Samples

Same input for cDNA for all community samples (25ng into 20 uL reaction).


```r
sub=subset(read, Frac_Comm=="Comm")
par(mar=c(5,5,0,1))
barplot(sub$Raw,las=2,names.arg = row.names(sub),horiz=T,cex.names=0.5,log='x', 
        xlab="Log10 No. Raw Reads")
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/raw-community-1.png)<!-- -->

## No. Reads Fractions 

Fractions had different RNA inputs, so different cDNA concentrations


```r
sub=subset(read, Frac_Comm=="Frac")
par(mar=c(5,5,0,1))
barplot(sub$Raw,las=2,names.arg = row.names(sub),horiz=T,cex.names=0.5,log='x', 
        xlab="Log No. Raw Reads")
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/raw-frac-1.png)<!-- -->

## No. Reads vs. cDNA Conc. 

No real correlation between cDNA/RNA concentration and the number of reads per sample. 


```r
plot(log10(sub$Raw), sub$Frac.RNA, xlab='Log10 No. of Raw Reads', ylab='cDNA concentration', pch =21, bg='black')
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/raw-cDNA-1.png)<!-- -->


## % Proks, Euks, + Cyano

Percent abundance of prokaryotic, eukaryotic, and cyanobacterial reads after sorting raw sequences based on hits to Silva and PR2 databases

On average, 1.8% of the total sequences were Eukaryotes

### Community % Proks, Euks, + Cyano
On average 0.46% of the raw reads were eukaryotes, so about 163/47,000 reads. 



```r
reads=read.csv(file='/Users/oliviaahern/Documents/GitHub/TestPlate/pct_euk_bact.csv',header=T,row.names=1)
dim(reads)
```

```
## [1] 95 17
```

```r
data1=data.frame(reads[,13:15])

sub=subset(reads, F_C=="Comm")
data1=data.frame(sub[,13:15])


par(mar=c(7,5,1,1),xpd=T)
barplot(as.matrix(t(data1)),las=2,names.arg = row.names(sub),horiz=T,cex.names=0.5,
        xlab="% Abundance", col=c('navy','forestgreen','firebrick'),space=0)
box(which='plot')
legend(0,-10, legend=c("Prok", "Cyano", "Euks"), pch =22, pt.bg=c("firebrick", 'forestgreen','navy'),
       bty='n', ncol=3)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/seqs-bins-1.png)<!-- -->


### Fractions % Proks, Euks + Cyanos
On average, 3.4% of the fractionated sequences were eukaryotes, so about 806/24,000 reads.

Percentage of total raw reads that match to prokaryotic, eukarytoic, or cyanobacterial databases. The samples that don't reach to 100 did not have matches to the databases, therefore we can assume they are probably chimera.s 


```r
sub=subset(reads, F_C=="Frac")
data1=data.frame(sub[,13:15])
barplot(as.matrix(t(data1)),las=2,names.arg = row.names(sub),horiz=T,cex.names=0.5,
        xlab="% Abundance", col=c('navy','forestgreen','firebrick'),space=0)
box(which='plot')
legend(0,-10, legend=c("Prok", "Cyano", "Euks"), pch =22, pt.bg=c("firebrick", 'forestgreen','navy'),
       bty='n', ncol=3)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/seqs-bins2-1.png)<!-- -->



# 18S sequences

Total of 110 unique ASVs. 


## read in data


```r
library(phyloseq)
x<-read.csv(file='/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/IMR_Jan2023/18S/feature-table.biom.csv',header=TRUE,row.names=1)
OTU = otu_table(x, taxa_are_rows=T)
taxa<-read.csv(file="/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/IMR_Jan2023/18S/taxonomy.csv",header=TRUE,row.names=1)
t<-as.matrix(taxa)
tax2<-tax_table(t)
map<-import_qiime_sample_data("/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/IMR_Jan2023/18S/map.txt")
#tree=read.tree("/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/tree.nwk")
phyo = phyloseq(OTU, tax2, map)
phyo
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 110 taxa and 84 samples ]
## sample_data() Sample Data:       [ 84 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 110 taxa by 9 taxonomic ranks ]
```


## overall ASV diversity

After QC - Average of 340 euk reads/sample + 3.4 ASVs per sample overall 

+ Community samples: avg 91 reads/sample + 2 ASVs/sample
+ Fractions: avg 600 reads/sample + 4 ASVs per sample 



```r
library(ggplot2)
#phyo_abund=subset_samples(phyo, Comm_Frac=="Comm")
phyo_abund=transform_sample_counts(phyo, function(x) x/sum(x))
pd <- psmelt(phyo_abund)


#colors=randomcoloR::randomColor(110)
#write_rds(colors, 'colors_class.rds')
colors=readRDS('/Users/oliviaahern/Documents/GitHub/TestPlate/colors_class_18s.rds')
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = ASV)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~Treatment, scale='free_x', space="free", shrink=TRUE) +
  labs(x=" ", y = "Class Relative Abundance") +
  theme_bw() +
  theme(axis.text = element_text(color ='black',size=12, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=14),
       legend.position = "none")
        
d_rgn
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/18s-asv-1.png)<!-- -->


## overall ASV diversity w/o Fungi



```r
phy=subset_taxa(phyo, Class!="Fungi")
#phyo_abund=subset_samples(phyo, Comm_Frac=="Comm")
phyo_abund=transform_sample_counts(phy, function(x) x/sum(x))
pd <- psmelt(phyo_abund)




#colors=randomcoloR::randomColor(110)
#write_rds(colors, 'colors_class.rds')
colors=readRDS('/Users/oliviaahern/Documents/GitHub/TestPlate/colors_class_18s.rds')
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Class)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~Treatment, scale='free_x', space="free", shrink=TRUE) +
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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/18s-asv-nof-1.png)<!-- -->

```r
#sub=subset_samples(phyo_abund, Comm_Frac=="Frac")
#write.csv(otu_table(sub),'18S_nofungi.csv')

#write.csv(tax_table(sub),'18S_nofungi_tax.csv')
```


## ASVs in fractions

*Note that the bar on the right is the heaviest of the fraction.*

We can begin to see differences in abundance at the heavier vs. lighter fractions with even a few 18S samples


```r
sub=subset_samples(phyo_abund, Comm_Frac=="Frac")

control <- subset_samples(sub,SampleID == "Ace.TB5.1" | SampleID == "Ace.TB5.4" | SampleID == "Ace.TB5.5" | SampleID == "Ace.TB5.15" | SampleID == "C12.TB5.7" | SampleID == "C12.TB5.11" | SampleID == "Met.TB5.3" | SampleID == "Met.TB5.9" | SampleID == "Met.TB5.10" | SampleID == "Met.TB5.12")

pd <- psmelt(control)

#colors=randomcoloR::randomColor(110)
#write_rds(colors, 'colors_class.rds')
colors=readRDS('/Users/oliviaahern/Documents/GitHub/TestPlate/colors_class_18s.rds')
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, BD), y = Abundance, fill = Class)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~Treatment, scale='free_x', space="free", shrink=TRUE) +
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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/18s-asvs-frac-1.png)<!-- -->


## My T0 Community vs. Old Siders Data, Phyla level

Relatively the same, except my data does not have any unknowns, potentially due to under sequencing. 

Sequencing Depth (seqs/sample): 

+ MC T0s 		6.9 ± 18.4
+ Axial 		2,810
+ Siders 		41,361 ± 5,659

Note: Opisthokonta is Fungi 


```r
#colors=randomcoloR::randomColor(110)
#write_rds(colors, 'colors_class.rds')
colors <- c("#bfd3e6","#fa9fb5", "#74c476", "#fc8d59", "#807dba", "#238443","#e31a1c", "#ec7014", "#969696")

read=read.csv(file='controlotu.csv',header=T,row.names=1)
colors <- c("#bfd3e6","#fa9fb5", "#74c476", "#fc8d59", "#807dba", "#238443","#e31a1c", "#ec7014", "#969696")

par(mar=c(15,5,1,10),xpd=T)
barplot(as.matrix(read), las=2, col = colors, space =0,
        ylab="% Relative Abundance")
box(which='plot')
legend(8,.7, legend=row.names(read), pch=22,
       pt.bg=colors)
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/mvsarah-1.png)<!-- -->

# 16S Data


Total of 1,632 ASVs, 1,608 ASVs if you get rid of Mitochondria and chloroplasts. 


```r
library(phyloseq)
x<-read.csv(file='/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/IMR_Jan2023/16S/feature-table.biom.csv',header=TRUE,row.names=1)
OTU = otu_table(x, taxa_are_rows=T)
taxa<-read.csv(file="/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/IMR_Jan2023/16S/taxonomy_16s.csv",header=TRUE,row.names=1)
t<-as.matrix(taxa)
tax2<-tax_table(t)
map<-import_qiime_sample_data("/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment4/Sequences/IMR_Jan2023/18S/map.txt")
#tree=read.tree("/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/tree.nwk")
phyo = phyloseq(OTU, tax2, map)
phyo1 = subset_taxa(phyo, !Order==" Chloroplast")
phyo2 = subset_taxa(phyo1, !Family==" Mitochondria")
phyo=phyo2
```


## Class Relative Abundance

### Community Class Relative Abundance


```r
agg=tax_glom(phyo2, taxrank="Class")

t=transform_sample_counts(agg, function(x) x/sum(x))
sub=subset_samples(t, Comm_Frac=="Comm")
phyo_abund=subset_samples(t, Comm_Frac=="Comm")

pd <- psmelt(phyo_abund)


#colors=randomcoloR::randomColor(52)
#write_rds(colors, 'colors_class.rds')
colors=readRDS('/Users/oliviaahern/Documents/GitHub/TestPlate/colors_class.rds')
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Class)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~Treatment, scale='free_x', space="free", shrink=TRUE) +
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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-comm-class-1.png)<!-- -->



### Fraction Class Relative Abundance

Two samples did not amplify: 

+ Ace-TB5-2
+ Met-TB5-4


```r
phyo_abund=subset_samples(t, Comm_Frac=="Frac")

pd <- psmelt(phyo_abund)


#colors=randomcoloR::randomColor(52)
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, BD), y = Abundance, fill = Class)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~Treatment, scale='free_x', space="free", shrink=TRUE) +
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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-class-frac-1.png)<!-- -->

### Fractions Top Classes No Repeats

Most Frequent Classes: 

+ Gammaproteobacteria 90%
+ Bacteroidia 84%
+ Alphaproteobacteria 81%
+ Desulfuromonadia 65%
+ Campylobacteria 58%




```r
read1=read.csv(file='/Users/oliviaahern/Documents/GitHub/TestPlate/fracs_class.csv',header=T,row.names=1)

data=data.frame(t(read1))
control=subset(data, Treatment=="Control")
meth=subset(data, Treatment=="Methanol")
acetate=subset(data,Treatment=="Acetate")


{
par(mar=c(5,5,1,1),mfrow=c(2,2))
plot(control$BD, control$Gammaproteobacteria,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),
     xlab="Buoyant Density", ylab="% Abundance Gammaproteobacteria")
lines(meth$BD, meth$Gammaproteobacteria,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$Gammaproteobacteria,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")
legend(1.76,.3, legend=c("Control", "Methanol", "Acetate"), text.col=c('gray70',"#a44f9a","#56ae6c"), bty='n')



plot(control$BD, control$Campylobacteria,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),xlab="Buoyant Density", ylab="% Abundance Campylobacteria")
lines(meth$BD, meth$Campylobacteria,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$Campylobacteria,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")


plot(control$BD, control$Alphaproteobacteria,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),xlab="Buoyant Density", ylab="% Abundance Alphaproteobacteria")
lines(meth$BD, meth$Alphaproteobacteria,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$Alphaproteobacteria,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")


plot(control$BD, control$Bacteroidia,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),xlab="Buoyant Density", ylab="% Abundance Bacteroidia")
lines(meth$BD, meth$Bacteroidia,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$Bacteroidia,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")

}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-classes-1.png)<!-- -->

## ASV Class Relative Abundance


### Community ASV Relative Abundance

```r
abund=transform_sample_counts(phyo, function(x) x/sum(x))
phyo_abund=subset_samples(abund, Comm_Frac=="Comm")

pd <- psmelt(phyo_abund)


#colors=randomcoloR::randomColor(1610)
#write_rds(colors, 'colors_class.rds')
colors=readRDS('/Users/oliviaahern/Documents/GitHub/TestPlate/colors_class.rds')
colors_asv=c('gray70', colors)
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Class)) +
  scale_fill_manual(values=as.character(t(colors_asv))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~Treatment, scale='free_x', space="free", shrink=TRUE) +
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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-asv-class-1.png)<!-- -->



### Fraction ASV Relative Abundance




```r
abund=transform_sample_counts(phyo, function(x) x/sum(x))
phyo_abund=subset_samples(abund, Comm_Frac=="Frac")

pd <- psmelt(phyo_abund)


#colors=randomcoloR::randomColor(1610)
#write_rds(colors, 'colors_class.rds')
colors=readRDS('/Users/oliviaahern/Documents/GitHub/TestPlate/colors_class.rds')
colors_asv=c('gray70', colors)
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Class)) +
  scale_fill_manual(values=as.character(t(colors_asv))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~Treatment, scale='free_x', space="free", shrink=TRUE) +
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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-asv-class-frac-1.png)<!-- -->

### Most Frequent ASVs in Fractions
Most Frequent ASVs

+ 16e2e2408cedf8f9faae33004fb5eaa2 Gammaproteobacteria Aeromonas 83%
+ X4785f91e7f34f642ba5880ff1f2e7526 Gammaproteobacteria Cellvibrio 70%
+ ddb2bf017eeb122792d652cc207cb3b0 Bacteroidia Bacteroidetes/Chlorobi 70%
+ X639b8d71f6806556291f6114837ecc41 Gammaproteobacteria Methylotenera 67%
+ X72043b4ee056432fd25f465446b2ff27 Gammaproteobacteria Methylotenera 67%
+ d1c8f7b6346015d2d94b7db3e6250ffc Alphaproteobacteria Rhodobacteraceae 65%

#### Abundance in Community 



```r
phyo_abund=transform_sample_counts(phyo, function(x) x/sum(x))
comm1=subset_samples(phyo_abund, Comm_Frac=="Comm")

s=subset_taxa(comm1, strain=="16e2e2408cedf8f9faae33004fb5eaa2" | strain=="4785f91e7f34f642ba5880ff1f2e7526" | strain=="ddb2bf017eeb122792d652cc207cb3b0" | strain=="639b8d71f6806556291f6114837ecc41" |
strain=="72043b4ee056432fd25f465446b2ff27" |  strain =="d1c8f7b6346015d2d94b7db3e6250ffc")


pd <- psmelt(s)



colors_asv=c("#801900",
"#ffc77e",
"#525b00",
"#b8c636",
"#013c9e",
"#dc74e6")
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = strain)) +
  scale_fill_manual(values=as.character(t(colors_asv))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~Treatment, scale='free_x', space="free", shrink=TRUE) +
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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-comm-asv-abund-1.png)<!-- -->

#### Fractions ASV vs. Buoyant Density

Some facts to note

+ [Cellvibrio can grow off of Xylose and glucose, but not Acetate](https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijs.0.02271-0). Also Cellivibrio used to be Pseudomonas
+ [The "Chlorobi" strain matched perfectly on blast to Saccharicrinis aurantiacus strain HQYD1 16S](https://www.ncbi.nlm.nih.gov/nucleotide/NR_156071.1?report=genbank&log$=nuclalign&blast_rank=1&RID=0PTR50JD01N) which was isolated from a sea squirt and can [use glucose and xylose as carbon sources](https://bacdive.dsmz.de/strain/133336)
+ [Methyloterna strain 1 matched to Methylotenera versatilis 301, a facultative anaerobe, mesophilic, bacterium isolated from freshwater sediment with 366/374 98% identity  ](https://bacdive.dsmz.de/strain/23100)
+ [Methyloterna strain 2 matched to Methylotenera versatilis 301, a facultative anaerobe, mesophilic, bacterium isolated from freshwater sediment with 100% identity  ](https://bacdive.dsmz.de/strain/23100)
+ The Rhodobacter strain matched to [Cypionkella sinensis strain Y1R2-4 16S which can use Xylose, Glucose, and Acetate as C sources](https://bacdive.dsmz.de/strain/141057)


```r
read1=read.csv(file='/Users/oliviaahern/Documents/GitHub/TestPlate/16S_asv_subset.csv',header=T,row.names=1)

control=subset(read1, Treatment=="Control")
meth=subset(read1, Treatment=="Methanol")
acetate=subset(read1,Treatment=="Acetate")


{
par(mar=c(5,5,1,1),mfrow=c(2,3))
plot(control$BD, control$X16e2e2408cedf8f9faae33004fb5eaa2,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),
     xlab="Buoyant Density", ylab="% Abundance Aeromonas")
lines(meth$BD, meth$X16e2e2408cedf8f9faae33004fb5eaa2,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$X16e2e2408cedf8f9faae33004fb5eaa2,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")
legend(1.76,.8, legend=c("Control", "Methanol", "Acetate"), text.col=c('gray70',"#a44f9a","#56ae6c"), bty='n')

plot(control$BD, control$X4785f91e7f34f642ba5880ff1f2e7526,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),
     xlab="Buoyant Density", ylab="% Abundance Cellvibrio")
lines(meth$BD, meth$X4785f91e7f34f642ba5880ff1f2e7526,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$X4785f91e7f34f642ba5880ff1f2e7526,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")
legend(1.76,.8, legend=c("Control", "Methanol", "Acetate"), text.col=c('gray70',"#a44f9a","#56ae6c"), bty='n')


plot(control$BD, control$ddb2bf017eeb122792d652cc207cb3b0,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),
     xlab="Buoyant Density", ylab="% Abundance Chlorobi")
lines(meth$BD, meth$ddb2bf017eeb122792d652cc207cb3b0,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$ddb2bf017eeb122792d652cc207cb3b0,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")
legend(1.76,.8, legend=c("Control", "Methanol", "Acetate"), text.col=c('gray70',"#a44f9a","#56ae6c"), bty='n')

plot(control$BD, control$X639b8d71f6806556291f6114837ecc41,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),
     xlab="Buoyant Density", ylab="% Abundance Methylotenera 1")
lines(meth$BD, meth$X639b8d71f6806556291f6114837ecc41,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$X639b8d71f6806556291f6114837ecc41,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")
legend(1.76,.8, legend=c("Control", "Methanol", "Acetate"), text.col=c('gray70',"#a44f9a","#56ae6c"), bty='n')


plot(control$BD, control$X72043b4ee056432fd25f465446b2ff27,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),
     xlab="Buoyant Density", ylab="% Abundance Methylotenera 2")
lines(meth$BD, meth$X72043b4ee056432fd25f465446b2ff27,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$X72043b4ee056432fd25f465446b2ff27,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")
legend(1.76,.8, legend=c("Control", "Methanol", "Acetate"), text.col=c('gray70',"#a44f9a","#56ae6c"), bty='n')

plot(control$BD, control$d1c8f7b6346015d2d94b7db3e6250ffc,
     type='o', pch=21,col='gray70', bg='gray90', ylim=c(0,1), xlim=c(1.76, 1.83),
     xlab="Buoyant Density", ylab="% Abundance Rhodobacter")
lines(meth$BD, meth$d1c8f7b6346015d2d94b7db3e6250ffc,
      type="o", col='#a44f9a', pch=21,
      bg="#a44f9a")
lines(acetate$BD, acetate$d1c8f7b6346015d2d94b7db3e6250ffc,
      type="o", pch=21, bg="#56ae6c", col="#56ae6c")
legend(1.76,.8, legend=c("Control", "Methanol", "Acetate"), text.col=c('gray70',"#a44f9a","#56ae6c"), bty='n')
}
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-frac-asv-abund-1.png)<!-- -->


## Fractions w/ multiple cycles

### Alpha Diversity 

#### Species Number 

```r
sub=subset_samples(phyo, Repeats=="Yes")
otu=otu_table(sub)
library(vegan)
spec=specnumber(t(otu))

par(mar=c(5,10,1,1))
barplot(spec, horiz=T, las=2, xlab='No. of Species')
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-repeats1-1.png)<!-- -->

#### Species Number vs. Cycles

There appears to be PCR bias in the nested PCR.



```r
cycles=c(25,35,40,50,55,25,35,40,50,55,25,35,40,50,55)
par(mar=c(5,5,1,1))
plot(cycles, spec, pch=21, bg='black', 
     xlab="No. of PCR Cycles", ylab="No. of Species")
```

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-repeats2-1.png)<!-- -->

### Class Level/ASV


```r
abund=transform_sample_counts(phyo, function(x) x/sum(x))
phyo_abund=subset_samples(abund, Comm_Frac=="Frac")
sub=subset_samples(phyo_abund, Repeats=="Yes")

pd <- psmelt(sub)

colors_asv=c('gray70', colors)
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Class)) +
  scale_fill_manual(values=as.character(t(colors_asv))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~Treatment, scale='free_x', space="free", shrink=TRUE) +
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

![](/Users/oliviaahern/Documents/GitHub/Exp4/docs/index_files/figure-html/16s-repeats-1.png)<!-- -->


