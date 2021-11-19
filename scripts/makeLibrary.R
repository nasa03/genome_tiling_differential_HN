#!/usr/bin/Rscript
##
##  makeLibrary.R
##
##  EDF 10/13/21
##

library(dplyr)
library(tidyr)

lib_name="PB-TL-0030"
left_ad = "ACTGGCCGCTTGACG"
right_ad = "CACTGCGGCTCCTGC"

tiled_seqs = read.table("library/all_tiled_seqs.txt",
                  header=TRUE,sep='\t')

carryover_seqs = read.table("library/carryover_seqs.txt",
                            header=TRUE, sep='\t')

tiled_seqs_table = tiled_seqs %>%
  mutate(PB_ID = paste0(lib_name,"-",
                        sprintf("%06d",as.numeric(row_number())))) %>%
  mutate(oligo = paste0(left_ad,`sequence`,right_ad),
         left_adaptor = left_ad,
         right_adaptor = right_ad,
         carryover='False',
         library=lib_name) %>%
  arrange(-abs(diff)) %>%
  distinct(substr(sequence,0,140), 
           .keep_all=TRUE)

carryover_seqs_table = carryover_seqs %>%
  mutate(sequence = substr(sequence,0,170)) %>%
  filter(nchar(sequence) == 170) %>%
  mutate(oligo = paste0(left_ad,`sequence`,right_ad),
         left_adaptor = left_ad,
         right_adaptor = right_ad,
         carryover='True')

all_seqs_table = tiled_seqs_table %>%
  rename(high_diff_pred=diff_pred) %>%
  select(PB_ID,sequence,
         left_adaptor,right_adaptor,oligo,
         carryover,library,
         region,spec,build,
         HepG2_pred=HepG2, N2a_pred=N2a, diff_pred=diff,
         high_HepG2_pred, high_N2a_pred, high_diff_pred) %>%
  bind_rows(carryover_seqs_table %>%
              select(PB_ID,sequence,
                 left_adaptor,right_adaptor,oligo,
                 carryover,library,
                 HepG2_obs, N2a_obs, diff_obs, 
                 HepG2_pred, N2a_pred, 
                 KI_high, KI_low, KI_diff, Ang_diff)) %>%
  mutate(seq_140 = substr(sequence,0,140),
         seq_100 = substr(sequence,0,100),
         seq_24 = substr(sequence,0,24),
         seq_12 = substr(sequence,0,12)) %>%
  filter(!grepl("N",sequence)) 


## Check for duplicate seqs (looks good)
nrow(all_seqs_table)
length(unique(all_seqs_table$sequence))
length(unique(all_seqs_table$seq_140))
length(unique(all_seqs_table$seq_100))
length(unique(all_seqs_table$seq_24))
length(unique(all_seqs_table$seq_12))

table(all_seqs_table %>% mutate(tmp=nchar(sequence)) %>% pull(tmp))

write.table(all_seqs_table,
            paste0("library/",lib_name,"_annotations.tsv"),
            col.names=TRUE,row.names=FALSE,
            sep='\t',quote=FALSE)


seqs_table = all_seqs_table %>%
  select(PB_ID,sequence,oligo,left_adaptor,right_adaptor,carryover)
write.table(seqs_table,
            paste0("library/",lib_name,"_seqs.tsv"),
            col.names=TRUE,row.names=FALSE,
            sep='\t',quote=FALSE)


oligos_table = seqs_table %>%
  mutate(out=paste0(">",PB_ID,"\n",oligo)) %>%
  select(out)
write.table(oligos_table,
            paste0("library/",lib_name,"_oligos.fasta"),
            row.names=FALSE,col.names=FALSE,
            quote=FALSE)

reference_table = seqs_table %>%
  mutate(out=paste0(">",PB_ID,"\n",sequence)) %>%
  select(out)
write.table(reference_table,
            paste0("library/",lib_name,"_reference.fasta"),
            row.names=FALSE,col.names=FALSE,
            quote=FALSE)
