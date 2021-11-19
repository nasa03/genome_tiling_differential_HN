#!/usr/bin/Rscript
##
##  compareScores.R
##
##  EDF 20211116
##


library(dplyr)
library(ggplot2)
library(tidyr)

print("reading files...")
N2a_human = read.table("scores/PB-PY-0011.hg38.scores.min.txt",
                       header=TRUE)
head(N2a_human)
# keep_human=sample(1:nrow(N2a_human), 10^6, replace=FALSE)
# N2a_human_10M = N2a_human[keep_human,]
# write.table(N2a_human_10M,
#             "scores/PB-PY-0011.hg38.scores.min.1Mrand.txt",
#             row.names=FALSE,col.names=TRUE,
#             sep='\t',quote=FALSE)
N2a_human_1M = read.table("scores/PB-PY-0011.hg38.scores.min.1Mrand.txt",
                          header=FALSE)

N2a_mouse = read.table("scores/PB-PY-0011.mm10.scores.min.txt",
                       header=TRUE)
head(N2a_mouse)
# keep_mouse = sample(1:nrow(N2a_mouse), 10^6, replace=FALSE)
# N2a_mouse_10M = N2a_mouse[keep_mouse,]
# write.table(N2a_mouse_10M,
#             "scores/PB-PY-0011.mm10.scores.min.1Mrand.txt",
#             row.names=FALSE,col.names=TRUE,
#             sep='\t',quote=FALSE)
N2a_mouse_1M = read.table("scores/PB-PY-0011.mm10.scores.min.1Mrand.txt",
                          header=FALSE)

HepG2_human = read.table("scores/PB-PY-0012.hg38.scores.min.txt",
                         header=TRUE)
head(HepG2_human)
# HepG2_human_10M = HepG2_human %>%
#   unite(var,chr,pos,sep="_") %>%
#   filter(var %in% (N2a_human_10M %>% unite(var,chr,pos,sep="_") %>% pull(var)))
HepG2_human_1M = read.table("scores/PB-PY-0012.hg38.scores.min.1Mrand.txt",
                            header=FALSE)

# write.table(HepG2_human_10M,
#             "scores/PB-PY-0012.hg38.scores.min.1Mrand.txt",
#             row.names=FALSE,col.names=TRUE,
#             sep='\t',quote=FALSE)

HepG2_mouse = read.table("scores/PB-PY-0012.mm10.scores.min.txt",
                         header=TRUE)
head(HepG2_mouse)
# HepG2_mouse_10M = HepG2_human %>%
#   unite(var,chr,pos,sep="_") %>%
#   filter(var %in% (N2a_human_10M %>% unite(var,chr,pos,sep="_") %>% pull(var)))
# write.table(HepG2_mouse_10M,
#             "scores/PB-PY-0012.mm10.scores.min.1Mrand.txt",
#             row.names=FALSE,col.names=TRUE,
#             sep='\t',quote=FALSE)
HepG2_mouse_1M = read.table("scores/PB-PY-0012.mm10.scores.min.1Mrand.txt",
                            header=FALSE)


rbind(N2a_human_1M %>% mutate(spec='human'),
      N2a_mouse_1M %>% mutate(spec='mouse')) %>%
  ggplot(aes(spec,V1)) +
  geom_violin() +
  geom_boxplot(width=.2, outlier.color=NA) +
  theme_classic() +
  ggtitle('Model 11 (Neuro2a) sequence predictions') +
  ylab('predicted activity')
ggsave("plots/N2a_11_humanvsmouse.1Mrand.pdf")

rbind(HepG2_human_1M %>% mutate(spec='human'),
      HepG2_mouse_1M %>% mutate(spec='mouse')) %>%
  ggplot(aes(spec,V1)) +
  geom_violin() +
  geom_boxplot(width=.2, outlier.color=NA) +
  theme_classic() +
  ggtitle('Model 12 (HepG2) sequence predictions') +
  ylab('predicted activity')
ggsave("plots/HepG2_12_humanvsmouse.1Mrand.pdf")

all_long = rbind(N2a_human_1M %>% mutate(spec='human', model='N2a'),
      N2a_mouse_1M %>% mutate(spec='mouse', model='N2a')) %>%
  rbind(HepG2_human_1M %>% mutate(spec='human', model='HepG2')) %>%
  rbind(HepG2_mouse_1M %>% mutate(spec='mouse', model='HepG2'))
all_long %>%
  ggplot(aes(model,V1)) +
  geom_hline(yintercept=0) +
  geom_violin(aes(color=spec),
              position=position_dodge(width=.7)) +
  # geom_boxplot(aes(color=spec),
  #              position=position_dodge(width=.7),
  #              width=.3,
  #              outlier.color=NA) +
  theme_classic() +
  ggtitle('Model 11 + 12 sequence predictions') +
  ylab('predicted activity')
ggsave("plots/both_11_12_humanvsmouse.1Mrand.pdf")


both_human = merge(N2a_human_1M %>% rename(N2a=V1),
                   HepG2_human_1M %>% rename(HepG2=V1),
                   by=c('chr','pos')) 
both_human %>%
  ggplot(aes(HepG2,N2a)) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_bin2d(aes(fill=log10(..count..))) +
  theme_bw() +
  scale_x_continuous(breaks=c(-8,-4,0,4,8)) +
  scale_y_continuous(breaks=c(-8,-4,0,4,8)) +
  ggtitle('Model 11 (Neuro2a) vs Model 12 (HepG2) predictions (human)')
ggsave("plots/both_11_12_human.1Mrand.pdf")

both_mouse = merge(N2a_mouse_1M %>% rename(N2a=V1),
                   HepG2_mouse_1M %>% rename(HepG2=V1),
                   by=c('chr','pos'))
both_mouse %>%
  ggplot(aes(HepG2,N2a)) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_bin2d(aes(fill=log10(..count..))) +
  theme_bw() +
  scale_x_continuous(breaks=c(-8,-4,0,4,8)) +
  scale_y_continuous(breaks=c(-8,-4,0,4,8)) +
  ggtitle('Model 11 (Neuro2a) vs Model 12 (HepG2) predictions (mouse)')
ggsave("plots/both_11_12_mouse.1M.pdf")


both_both = rbind(both_mouse %>% mutate(spec='mouse'),
                  both_human %>% mutate(spec='human'))
both_both %>%
  ggplot(aes(HepG2,N2a)) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_bin2d(aes(fill=log10(..count..))) +
  theme_bw() +
  scale_x_continuous(breaks=c(-8,-4,0,4,8)) +
  scale_y_continuous(breaks=c(-8,-4,0,4,8)) +
  facet_wrap(~spec) +
  ggtitle('Model 11 (Neuro2a) vs Model 12 (HepG2) predictions')
ggsave("plots/both_11_12_both.1Mrand.pdf")

diffs = both_both %>%
  mutate(diff=HepG2-N2a)

diffs %>% 
  filter(abs(diff) > 1) %>%
  ggplot(aes(HepG2,N2a)) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_bin2d(aes(fill=log10(..count..))) +
  theme_bw() +
  scale_x_continuous(breaks=c(-8,-4,0,4,8)) +
  scale_y_continuous(breaks=c(-8,-4,0,4,8)) +
  facet_wrap(~spec) +
  ggtitle('Model 11 (Neuro2a) vs Model 12 (HepG2) predictions (diff only)')
ggsave("plots/both_11_12_both_diffonly.1Mrand.pdf")

