library(seqinr)
library(tidyverse)

##########################################################################

# fasta to dataframe

fastas <- read.fasta('gapn_hits.fasta', seqtype='AA', as.string=TRUE)

read_fasta <- function(fasta) {

    seq <- fasta[1]
    name <- attr(fasta, 'name')
    annot <- attr(fasta, 'Annot')

    # need to escape \
    # use match to get groups, extract to get direct match
    accession <- str_match(annot, '\\|(\\w+)\\|')[1, 2]
    enz_name <- str_split(str_split(annot, ' ', n=2)[[1]][2], 'OS=')[[1]][1]
    species_name <- str_match(annot, 'OS=(.*)\\sOX')[1, 2]
    gene_name <- str_match(annot, 'GN=(\\w+)')[1, 2]

    df <- tibble(name = name,
                 accession = accession,
                 enz = enz_name,
                 gene = gene_name,
                 species = species_name,
                 seq = seq)

    return (df)
}

df_fasta <- map(fastas, read_fasta) %>%
    bind_rows()# %>%
    #write.csv('gapn_df_r.csv', row.names = FALSE)

##########################################################################

# align sequences with mafft

#system('mafft gapn_hits.fasta > gapn_mafft.fasta')

##########################################################################

# remove gaps in alignment according to template sequence

gapn_template = read.fasta('gapn_mafft.fasta', seqtype='AA', as.string=TRUE)[[1]]

get_gaps <- function(fasta) {

    seq <- fasta[1]
    gap_idx <- str_locate_all(seq, '[^-]')
    return (gap_idx[[1]])
}

fastas_aln <- read.fasta('gapn_mafft.fasta', seqtype='AA', as.string=TRUE)
gap_idx <- get_gaps(gapn_template)

ungap <- function(fasta, gap_idx) {

    seq <- fasta[1]
    seq_clean <- str_sub(seq, gap_idx)
    seq_clean <- paste(seq_clean, collapse='')

    fasta[1] <- seq_clean

    return (fasta)
}

# remove first >
annot <- map(fastas_aln, attr, 'Annot') %>%
    str_sub(2, str_length(.))

seq_clean <- map(fastas_aln, ungap, gap_idx)
write.fasta(seq_clean, names=annot, file.out='gapn_ungap.fasta', nbchar = 60, as.string = TRUE)
seq_clean <- read.fasta('gapn_ungap.fasta', seqtype = 'AA', as.string=TRUE)


#########################################################################

# calculate sequence identity to template

gapn_seq <- fastas[[1]][1]

get_seqid <- function(fasta, template) {

    seq_len <- str_length(template)
    mismatches <- adist(fasta[1], template)[1, 1]
    name <- attr(fasta, 'name')

    seq_id <- 1 - (mismatches / seq_len)

    df <- tibble(name = name,
                 seq_id = seq_id)

    return (df)
}

df_seqid <- map(seq_clean, get_seqid, gapn_seq) %>%
    bind_rows() %>%
    inner_join(df_fasta) %>%
    relocate(seq_id, .after='gene') %>%
    distinct(name, .keep_all=TRUE) %>%
    arrange(desc(seq_id))

write.csv(df_seqid, 'gapn_df_r.csv', row.names=FALSE)

##########################################################################

# get residue identity at each position

get_resis <- function(fasta){
    
    seq <- fasta[1]
    seq_split <- str_extract_all(seq, boundary('character'))
    
    df <- tibble(unlist(seq_split))
    name <- attr(fasta, 'name')
    accession <- str_match(name, '\\|(\\w+)\\|')[1, 2]
    names(df) <- accession
    
    return (df)
}

each_pos_aa <- map(seq_clean, get_resis) %>%
    bind_cols()

each_pos_aa <- as_tibble(cbind(sample=names(each_pos_aa), t(each_pos_aa)))

col_names <- c('sample', seq(1:(ncol(each_pos_aa) - 1)))
names(each_pos_aa) <- col_names
each_pos_aa <- column_to_rownames(each_pos_aa, 'sample')

# to select by presence of residue, all sequences with M at position 1
# filter(each_pos_aa, `1` == 'M')

write.csv(each_pos_aa, 'gapn_each_pos_aa.csv')


##########################################################################

# residue frequencies

# slow, need to rewrite without the for loop
get_freqs <- function(each_pos_aa) {
    
    colnames <- names(each_pos_aa)

    for (i in seq_along(each_pos_aa)) {
        
        #each_pos_aa[, i] <- str_replace_all(each_pos_aa[, i], 'X', '-')
        counts <- table(each_pos_aa[, i])
        names <- attr(counts, 'dimnames')
        
        df <- tibble(residues = unlist(names),
                     counts = counts)
        names(df) <- c('residues', colnames[i])
        
        if (i == 1) {
            res <- df
        } else {
            res <- full_join(res, df)
        }
            
    }
    
    out <- res %>%
        mutate(across(everything(), ~replace_na(.x, 0))) %>%
        filter(residues != '-') %>%
        filter(residues != 'X') %>%
        mutate_if(is.numeric, list(~./sum(.)*100))
        
    return (out)
}


freqs <- get_freqs(each_pos_aa)
write.csv(freqs, 'gapn_freqs.csv', row.names=FALSE)
