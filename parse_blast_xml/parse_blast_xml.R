library(tidyverse)
library(stringr)
library(rvest)
library(XML)
library(Biostrings)

# blastp with uniprot/siwssprot db
xml <- xmlParse('SJG1UGX5016-Alignment.xml')

# extract xml data from hits
hit_num <- xpathApply(xml, '//Hit/Hit_num', xmlValue)
hit_id <- xpathApply(xml, '//Hit/Hit_id', xmlValue)
hit_accession <- xpathApply(xml, '//Hit/Hit_accession', xmlValue)
hit_def <- xpathApply(xml, '//Hit/Hit_def', xmlValue)

align_len <- xpathApply(xml, '//Hit/Hit_hsps/Hsp/Hsp_align-len', xmlValue)
hsp_ident <- xpathApply(xml, '//Hit/Hit_hsps/Hsp/Hsp_identity', xmlValue)
hsp_pos <- xpathApply(xml, '//Hit/Hit_hsps/Hsp/Hsp_positive', xmlValue)

# seq length of query
query_len = as.numeric(unlist(xpathApply(xml, '//Iteration/Iteration_query-len', xmlValue)))

df <- tibble(hit_num = unlist(hit_num), 
             hit_id = unlist(hit_id),
             hit_accession = unlist(hit_accession),
             hit_def = unlist(hit_def),
             align = unlist(map(align_len, as.numeric)) / query_len,
             identity = unlist(map(hsp_ident, as.numeric)) / query_len,
             pos = unlist(map(hsp_pos, as.numeric)) / query_len)

# enzyme name
get_full_name <- function(x) {
    name = str_split(x, ';')[[1]][1]
    name = str_split(name, '=')[[1]][2]
}

# organism
# remove brackets
get_organism <- function(x) {
    last = str_split(x, ';')[[1]]
    last = last[length(last)]
    org = str_extract(last, '\\[.*\\]')
    org = str_replace(org, '\\[', '')
    org = str_replace(org, '\\]', '')
}

orgs = map(df$hit_def, get_organism)
full = map(df$hit_def, get_full_name)

df$orgs <- unlist(orgs)
df$full <- unlist(full)

# save blast output
write.csv(df, 'blast.csv', row.names = FALSE)


##################################################################################

# download fasta sequences

download_fasta <- function(accession) {
    url <- paste0('https://www.uniprot.org/uniprot/', accession, '.fasta')
    fasta <- xml2::read_html(url) %>%
        html_text()
    
    # split on first newline, remove newlines in sequence
    #fasta_split = str_split(fasta, "\n", n=2)
    #seq_desc = fasta_split[[1]][1]
    #seq = str_replace_all(fasta_split[[1]][2], "\n", '')
    
    # write fasta output
    out = paste0(accession, '.fasta')
    cat(fasta, file=out)

}

map(df$hit_accession, download_fasta)
