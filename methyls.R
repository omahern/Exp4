s=subset_taxa(phyo, Family =="Methylophilaceae" & Order=="Methylococcales")

abund=transform_sample_counts(phyo, function(x) (x/sum(x)))
list=c('cce70c631b2db81a7d46759431028440','67946ba2bcab24f67feffe65d12fb6be','ae6739901982a2ea3424bf8090ff5a23','21f656ab80a97a6a77337e7d7c178542','27d76641734cb206f22b18adec07a35a','3518223fe974cf97df77024998254d5a','49fae5227c15913fa3e00ffd09a5cbef','70652fafc1dc64f9f0d4d4e1087700e4','d7972c27c10d043696a797f99f1635e8','e96be20c155b93c31a23ea5cc2fd01cd','fc169a62076539bcbb44364586da6a81','263ab46f41dea694af7b411e17133aa6','2907f7cc27f892f56fb9b80ea2ca5250','320063ed141bd04436f82761120387a9','35334e2fbbc14d27479474096d0fb8f9','45cc071db1d6be0adec13086e6ea46f1','72ceb780c4a419f55a1839f39a95d20a','72ee1348e6876fc3387ce405317b7359','ac5abd8b75a074d05c9f1b922e272445','c4debfffb5552abc5445848e898e9b3b','c8b292289369a4666964261b1ddb7827','d3b330d52ce6cbc3a6259ccc1444a65c','d6c8ff4a795713fac330832cee2ccdc5','dfa4024993b2772bac9172cae6b7ce33','e2c65d17b5a2497b63b32e7de74d8fde','e76a3ed55411fdced6fd74afb1c83490','f1bef73774e8ee3a5ff8009521a6d990','f94fdd288cd02629ee46e176835cafef','f345eaf38b04cf0ae07d05751f6441e8','110f699baf2b0b328410b3ad634f9e1c','81725d931c79fed5a0ecd45fe69fa66e','8fed996fbe127c0dbd9401c26801bf2b','bd5b2355e067c418bb310d3be315c958','be183a83929b1197e75d5a3641d4f53f','c1453a490ba9da64b0b208d5457c7760','c97da5c097f5790a8d5bed843427371d','dfa4d6145f2b1205c94aab72c85b2a43','e45711696a4c1da50e19910a85dc4ae8','e72869dee81a7fb6d21208cf7401d822','f1049c3750b69fdbaa37ee1c136158d5','24f8405cbb0234bccb868c8fa3db1803','320184fed9bea185a1a17f70b49d2479','3cd7863b3a24b807deaf59bcae7b4889','5d17d6edaf996823b77d05b3fd663792','65e479befe53e1f77b10166317d058ff','791430fdb1a404d801f8b4777d103256','81755d85d2ad09a9cdd59a5844e940cc','9c7537874de2afb051b566148b70306a','ebc7a5b542b40a778048496c305923d8','eeb4e9080dd338e3fb5783231223f6a7','f5eafa997442a3aaa9eb0c4d2fb3925a','80739c294da8ce6a722c2405051673c8','9f360d919f74d2232a5a24f17ff95bd1','a6fe2e2d740b1011d9f21c6ab9b99890')
my_subset <- subset(otu_table(abund), rownames(otu_table(abund)) %in% list)
new_physeq <- merge_phyloseq(my_subset, tax_table(abund), sample_data(abund))

batch=subset_samples(new_physeq, B_C=="Batch")
o=otu_table(batch)
spec=specnumber(t(o))
pd <- psmelt(batch)
b1=subset_samples(new_physeq, B_C=="Batch")
pd <- psmelt(batch)



colors=randomcoloR::randomColor(54)
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Strain)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) +
  labs(x=" ", y = "Methyl* ASV Relative Abundance") +
  theme_bw() +
  theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
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

chemo=subset_samples(new_physeq, B_C=="Chemostat")
pd2 <- psmelt(chemo)
pp = ggplot(pd2, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Strain)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) +
  labs(x=" ", y = "Methyl* ASV Relative Abundance") +
  theme_bw() +
  theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=14),
        legend.position = "none",
        panel.grid=element_blank(),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
        strip.background = element_rect(
          color="white", fill="white", size=1, linetype="solid"),
        panel.spacing = unit(0.05, "lines"))
pp




?cowplot



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
  theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
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


abund=transform_sample_counts(phyo, function(x) (x/sum(x)))
list=c('cce70c631b2db81a7d46759431028440','67946ba2bcab24f67feffe65d12fb6be','ae6739901982a2ea3424bf8090ff5a23','21f656ab80a97a6a77337e7d7c178542','27d76641734cb206f22b18adec07a35a','3518223fe974cf97df77024998254d5a','49fae5227c15913fa3e00ffd09a5cbef','70652fafc1dc64f9f0d4d4e1087700e4','d7972c27c10d043696a797f99f1635e8','e96be20c155b93c31a23ea5cc2fd01cd','fc169a62076539bcbb44364586da6a81','263ab46f41dea694af7b411e17133aa6','2907f7cc27f892f56fb9b80ea2ca5250','320063ed141bd04436f82761120387a9','35334e2fbbc14d27479474096d0fb8f9','45cc071db1d6be0adec13086e6ea46f1','72ceb780c4a419f55a1839f39a95d20a','72ee1348e6876fc3387ce405317b7359','ac5abd8b75a074d05c9f1b922e272445','c4debfffb5552abc5445848e898e9b3b','c8b292289369a4666964261b1ddb7827','d3b330d52ce6cbc3a6259ccc1444a65c','d6c8ff4a795713fac330832cee2ccdc5','dfa4024993b2772bac9172cae6b7ce33','e2c65d17b5a2497b63b32e7de74d8fde','e76a3ed55411fdced6fd74afb1c83490','f1bef73774e8ee3a5ff8009521a6d990','f94fdd288cd02629ee46e176835cafef','f345eaf38b04cf0ae07d05751f6441e8','110f699baf2b0b328410b3ad634f9e1c','81725d931c79fed5a0ecd45fe69fa66e','8fed996fbe127c0dbd9401c26801bf2b','bd5b2355e067c418bb310d3be315c958','be183a83929b1197e75d5a3641d4f53f','c1453a490ba9da64b0b208d5457c7760','c97da5c097f5790a8d5bed843427371d','dfa4d6145f2b1205c94aab72c85b2a43','e45711696a4c1da50e19910a85dc4ae8','e72869dee81a7fb6d21208cf7401d822','f1049c3750b69fdbaa37ee1c136158d5','24f8405cbb0234bccb868c8fa3db1803','320184fed9bea185a1a17f70b49d2479','3cd7863b3a24b807deaf59bcae7b4889','5d17d6edaf996823b77d05b3fd663792','65e479befe53e1f77b10166317d058ff','791430fdb1a404d801f8b4777d103256','81755d85d2ad09a9cdd59a5844e940cc','9c7537874de2afb051b566148b70306a','ebc7a5b542b40a778048496c305923d8','eeb4e9080dd338e3fb5783231223f6a7','f5eafa997442a3aaa9eb0c4d2fb3925a','80739c294da8ce6a722c2405051673c8','9f360d919f74d2232a5a24f17ff95bd1','a6fe2e2d740b1011d9f21c6ab9b99890')
my_subset <- subset(otu_table(phyo), rownames(otu_table(phyo)) %in% list)
new_physeq <- merge_phyloseq(my_subset, tax_table(phyo), sample_data(phyo))

b1=subset_samples(new_physeq, B_C=="Batch")
#phylosmith::set_sample_order(b1, sort_on = "Timepoint")
bar1=plot_richness(b1, measures="Observed") + 
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) + 
  theme_bw() + geom_bar(stat = "identity", width=0.97, color='black',
                        lwd=0.1) +
   theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=14),
        legend.position = "none",
        panel.grid=element_blank(),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
        strip.background = element_rect(
          color="white", fill="white", size=1, linetype="solid"),
        panel.spacing = unit(0.05, "lines")) +
  labs(x=" ", y = "No. Methy* ASVs") 

b2=subset_samples(new_physeq, B_C=="Chemostat")
b=phylosmith::set_sample_order(b2, sort_on = c("Timepoint","MC"))

bar2=plot_richness(b2, measures="Observed") + 
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) + 
  theme_bw() +  geom_bar(stat = "identity", width=0.97, color='black',
                         lwd=0.1) +
  theme(axis.text = element_text(color ='black',size=12, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=14),
        legend.position = "none",
        panel.grid=element_blank(),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
        strip.background = element_rect(
          color="white", fill="white", size=1, linetype="solid"),
        panel.spacing = unit(0.05, "lines")) +
  labs(x=" ", y = "No. Methy* ASVs") 
bar2

plot_grid( d_rgn,pp, align="v",axis="b", ncol=1)

#plot_grid(bar1, d_rgn, align="h", axis='h', cols=1)


1=subset_samples(new_physeq, B_C=="Batch")
pd <- psmelt(batch)



colors=randomcoloR::randomColor(54)
library(grid)
#palette(as.character(t(colors)))
d_rgn = ggplot(pd, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Strain)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) +
  labs(x=" ", y = "Methyl* ASV Relative Abundance") +
  theme_bw() +
  theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=10),
        legend.position = "none",
        panel.grid=element_blank(),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
        strip.background = element_rect(
          color="white", fill="white", size=1, linetype="solid"),
        panel.spacing = unit(0.05, "lines"))
d_rgn

chemo=subset_samples(new_physeq, B_C=="Chemostat")
pd2 <- psmelt(chemo)
pp = ggplot(pd2, aes(x = reorder(Sample, Timepoint), y = Abundance, fill = Strain)) +
  scale_fill_manual(values=as.character(t(colors))) +
  geom_bar(stat = "identity", width=0.97, color='black',
           lwd=0.1) +
  facet_grid(~MC, scale='free_x', space="free", shrink=TRUE) +
  labs(x=" ", y = "Methyl* ASV Relative Abundance") +
  theme_bw() +
  theme(axis.text = element_text(color ='black',size=8, hjust=1,vjust=1),
        axis.text.x=element_text(size=8, angle=90),
        axis.title = element_text(color='black',face='bold',size=10),
        legend.position = "none",
        panel.grid=element_blank(),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
        strip.background = element_rect(
          color="white", fill="white", size=1, linetype="solid"),
        panel.spacing = unit(0.05, "lines"))
pp



plot_grid( d_rgn,pp, align="v",axis="b", ncol=1)




