##############################################################################################
library(tidyverse)
fastq <- list.files("/Users/arbones/Dropbox/SyncBriefcase/LAB/UK/APOE_switch/Georgia_bulk/fastq",
                    recursive = T, pattern = "fq.gz") |> 
    as_tibble()

    
##########################################

## Sample sheet for nf-core:RNAseq
path = "/scratch/jar301/georgia/fastq/"

lps <- fastq |>  
    mutate(sample = str_split_i(value,pattern = "/",1) |> 
               str_remove("_")) |> 
    mutate(value =paste0(path,value),
           fastq = rep(c("fastq_1","fastq_2"),20)) |> 
    pivot_wider(names_from = fastq, values_from = value) |> 
    mutate(strandedness = "auto") |> 
    select(sample,fastq_1,fastq_2,strandedness)


write_csv(lps,file = "20250907-sample_sheet.csv")
