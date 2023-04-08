library(dplyr)
library(tidyr)
library(bio3d)

template_file <- unlist(snakemake@input)
output_file <- unlist(snakemake@output)

templates <- read.table("pdb/template.txt", col.names = c("ID", "_P_", "fname")) %>%
    mutate(ID = substring(ID, 2))
data <- lapply(templates$fname, read.pdb) %>%
    setNames(templates$ID)
sheets <- lapply(data, `[[`, "sheet") %>%
    lapply(data.frame) %>%
    setNames(templates$ID) %>%
    bind_rows(.id = "ID") %>%
    mutate(sense = ifelse("sense" %in% names(.), as.numeric(sense), NA), sense = ifelse(sense < 0, "-", "+"))
helices <- lapply(data, `[[`, "helix") %>%
    lapply(data.frame) %>%
    setNames(templates$ID) %>%
    bind_rows(.id = "ID")
list(sheet = sheets, helix = helices) %>%
    bind_rows(.id = "ss") %>%
    filter(chain == "A") %>%
    replace_na(list(sense = ".")) %>%
    mutate(source = "SS", score = ".", frame = ".", attrib = ".") %>%
    select(ID, source, ss, start, end, score, sense, frame, attrib) %>%
    write.table(output_file, quote = F, sep = "\t", col.names = F, row.names = F)
