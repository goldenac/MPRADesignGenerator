# CUSTOM MPRA PROCESSOR #############################################################################################

#source("processor_functions.R")

#####################################################################################################################

#' Generate Design File
#'
#' To actually generate the file containing all oligos in the library, run this function.
#'
#'
#' @param variant input path, tags path, scrambled variants path
#' @return design file
#' @export
generate = function(variant_input_path, tag_path, scrambled_path){
  # 1) LOAD FILTERED TAGS #############################################################################################
  tags <- readr::read_csv(tag_path, col_names=TRUE, col_types=cols("c"))
  
  # 2) CHECK TAGS FOR DIGESTION SITES #################################################################################
  tag_check <- tags
  tag_check$ez1 <- grepl("GGTACC", tag_check$barcode)
  tag_check$ez2 <- grepl("TCTAGA", tag_check$barcode)
  tag_check$ez3 <- grepl("GGCC.....GGCC", tag_check$barcode)
  bad_tag_rows <- which(tag_check$ez1==TRUE|tag_check$ez2==TRUE|tag_check$ez3==TRUE)
  if(length(bad_tag_rows)!=0)
  {
    tags <- tags[-c(bad_tag_rows),]
  }
  # this section should be tested before being used as our tag set had no digestion sites so no tags were removed
  
  # 2) CREATE SCRAMBLED SEQUENCES #####################################################################################
  scrambled_vars <- readr::read_csv(scrambled_path, col_names=TRUE, col_types=cols("c","c"))
  
  # separate tags to be used for scrambled sequences
  tag_number <- nrow(tags)
  scrambled_tags_needed <- nrow(scrambled_vars)*25
  scrambled_tags <- tags[1:scrambled_tags_needed,]
  tags <- tags[8601:tag_number,]
  
  fwdprimer <- "ACTGGCCGCTTCACTG"
  fwdspacer <- "TG"
  enzyme1 <- "GGTACC"
  enzyme2 <- "TCTAGA"
  revspacer <- "GGC"
  revprimer <- "AGATCGGAAGAGCGTCG"
  
  scrambled_vars <- scrambled_vars %>% dplyr::slice(rep(1:n(), each=25))
  
  # format scrambled_vars to match complete_output
  scrambled_vars$CHROM <- "0"
  scrambled_vars$POS <- "0"
  scrambled_vars$ID <- scrambled_vars$rs_number
  scrambled_vars$REF <- "ref"
  scrambled_vars$ALT <- "alt"
  scrambled_vars$seq <- scrambled_vars$seq145
  scrambled_vars$class <- "scrambled"
  scrambled_vars$tag <- scrambled_tags$barcode
  scrambled_vars$construct <- paste0(fwdprimer, fwdspacer, scrambled_vars$seq145, enzyme1, enzyme2, scrambled_tags$barcode, revspacer, revprimer)
  
  scrambled_vars <- subset(scrambled_vars, select=-c(rs_number, seq145))
  
  # 2) READ IN VARIANT INPUT ##########################################################################################
  all_variants <- readr::read_csv(variant_input_path, col_names=TRUE, col_types=cols("c","d","c","c","c"))
  
  # 3) GET SEQUENCES ##################################################################################################
  genome <- BSgenome.Hsapiens.UCSC.hg38
  all_variants$REFseq <- retrieveSeq(genome, all_variants$CHROM, all_variants$POS)
  
  # 3.5) GENERATE A LIST OF SEQUENCES THAT HAVE DIGESTION SITES #######################################################
  seqs_w_digest <- all_variants
  seqs_w_digest$KpnI <- grepl("GGTACC", seqs_w_digest$REFseq)
  seqs_w_digest$XbaI <- grepl("TCTAGA", seqs_w_digest$REFseq)
  seqs_w_digest$SfiI <- grepl("GGCC.....GGCC", seqs_w_digest$REFseq)
  seqs_w_digest <- dplyr::filter(seqs_w_digest, KpnI=="TRUE"|XbaI=="TRUE"|SfiI=="TRUE")
  seqs_w_digest <- subset(seqs_w_digest, select=c(ID, REFseq, KpnI, XbaI, SfiI))
  #write.csv(seqs_w_digest, "sequences_w_digestion_sites.csv", row.names=FALSE)
  
  # 4) DETECT AND REPAIR DIGESTION SITES ##############################################################################
  # By searching for & repairing sites here, we reduce the number of operations the program must perform,
  # but we also must assume that creating the alt sequence does not introduce new digestion sites.
  # This is currently acceptable since we have already removed those variants with this problem, but in the future
  # this check should be done for both ref and alt sequences (if only 1 has site, alt causes prob), not just ref.
  
  all_variants$REFseq <- gsub("GGTACC", "GGATCC", all_variants$REFseq)
  all_variants$REFseq <- gsub("TCTAGA", "TCATGA", all_variants$REFseq)
  
  # SfiI (will only repair one site per sequence)
  sfii <- "GGCC.....GGCC"
  all_variants$sfii <- grepl(sfii, all_variants$REFseq)
  variants_wo_sfii <- dplyr::filter(all_variants, sfii=="FALSE")
  variants_w_sfii <- dplyr::filter(all_variants, sfii=="TRUE")
  variants_w_sfii_fixed <- variants_w_sfii
  variants_w_sfii$coords <- str_locate(variants_w_sfii$REFseq, "GGCC.....GGCC")
  variants_w_sfii$coord1 <- variants_w_sfii$coords[,"start"]
  variants_w_sfii$coord2 <- variants_w_sfii$coords[,"end"]
  variants_w_sfii_fixed$REFseq <- fixSfiI(variants_w_sfii$REFseq, variants_w_sfii$coord1, variants_w_sfii$coord2)
  all_variants <- rbind(variants_wo_sfii, variants_w_sfii_fixed)
  all_variants <- subset(all_variants, select=-sfii)
  
  # 5) SEPARATE SNPS, INSERTIONS, AND DELETIONS #######################################################################
  insertions <- dplyr::filter(all_variants, REF=="-")
  deletions <- dplyr::filter(all_variants, ALT=="-")
  snps <- dplyr::filter(all_variants, REF!="-" & ALT!="-")
  
  # #) COMPARE WHAT WILL BE CHANGED TO WHAT SHOULD BE CHANGED #########################################################
  # a way of checking that the coordinates are correct
  
  deletion_coord_check <- deletions
  deletion_coord_check$deleted_bases <- delBases(deletion_coord_check$REFseq, deletions$REF)
  deletion_coord_check$bases_match <- ifelse(deletion_coord_check$REF==deletion_coord_check$deleted_bases, "YES", "NO")
  wrong_coords <- filter(deletion_coord_check, bases_match=="NO")
  print(wrong_coords)
  # there are no deletions with the wrong coords!
  
  # 6) CREATE ALT SEQ FOR DELETIONS ###################################################################################
  deletions$ALTseq <- altDel(deletions$REFseq, deletions$REF)
  
  # 7) CREATE ALT SEQ FOR INSERTIONS ##################################################################################
  # PREVIOUS VERSION WAS WRONG!!!!!! (SHORTENTED REF BEFORE MAKING ALT -> INSERTION OCCURRED AT WRONG PLACE)
  insertions$ALTseq <- altIns(insertions$REFseq, insertions$ALT)
  insertions$REFseq <- refIns(insertions$REFseq, insertions$ALT)
  
  # 8) CREATE ALT SEQ FOR SNPS ########################################################################################
  snps$ALTseq <- altSnps(snps$REFseq, snps$ALT)
  
  # 9) GENERATE REV COMP SEQ FOR DELETIONS ############################################################################
  rev_deletions <- deletions
  
  rev_deletions_ref_DNAset <- Biostrings::DNAStringSet(rev_deletions$REFseq)
  rev_deletions_ref_DNAset <- Biostrings::reverseComplement(rev_deletions_ref_DNAset)
  rev_deletions_ref_DNAset <- data.frame(REFseq=as.character(rev_deletions_ref_DNAset))
  
  rev_deletions_alt_DNAset <- Biostrings::DNAStringSet(rev_deletions$ALTseq)
  rev_deletions_alt_DNAset <- Biostrings::reverseComplement(rev_deletions_alt_DNAset)
  rev_deletions_alt_DNAset <- data.frame(ALTseq=as.character(rev_deletions_alt_DNAset))
  
  rev_deletions$REFseq <- rev_deletions_ref_DNAset$REFseq
  rev_deletions$ALTseq <- rev_deletions_alt_DNAset$ALTseq
  
  # 10) GENERATE REV COMP SEQ FOR INSERTIONS ###############################################################################
  rev_insertions <- insertions
  
  rev_insertions_ref_DNAset <- Biostrings::DNAStringSet(rev_insertions$REFseq)
  rev_insertions_ref_DNAset <- Biostrings::reverseComplement(rev_insertions_ref_DNAset)
  rev_insertions_ref_DNAset <- data.frame(REFseq=as.character(rev_insertions_ref_DNAset))
  
  rev_insertions_alt_DNAset <- Biostrings::DNAStringSet(rev_insertions$ALTseq)
  rev_insertions_alt_DNAset <- Biostrings::reverseComplement(rev_insertions_alt_DNAset)
  rev_insertions_alt_DNAset <- data.frame(ALTseq=as.character(rev_insertions_alt_DNAset))
  
  rev_insertions$REFseq <- rev_insertions_ref_DNAset$REFseq
  rev_insertions$ALTseq <- rev_insertions_alt_DNAset$ALTseq
  
  # 11) GENERATE REV COMP SEQ FOR SNPS #####################################################################################
  rev_snps <- snps
  
  rev_snps_ref_DNAset <- Biostrings::DNAStringSet(rev_snps$REFseq)
  rev_snps_ref_DNAset <- Biostrings::reverseComplement(rev_snps_ref_DNAset)
  rev_snps_ref_DNAset <- data.frame(REFseq=as.character(rev_snps_ref_DNAset))
  
  rev_snps_alt_DNAset <- Biostrings::DNAStringSet(rev_snps$ALTseq)
  rev_snps_alt_DNAset <- Biostrings::reverseComplement(rev_snps_alt_DNAset)
  rev_snps_alt_DNAset <- data.frame(ALTseq=as.character(rev_snps_alt_DNAset))
  
  rev_snps$REFseq <- rev_snps_ref_DNAset$REFseq
  rev_snps$ALTseq <- rev_snps_alt_DNAset$ALTseq
  
  # 12) BIND ALL DATA TOGETHER ########################################################################################
  fwd_deletions_ref <- subset(deletions, select=-ALTseq)
  fwd_deletions_ref$class <- "fwd_ref"
  fwd_deletions_ref <- fwd_deletions_ref %>% dplyr::rename(seq=REFseq)
  
  fwd_deletions_alt <- subset(deletions, select=-REFseq)
  fwd_deletions_alt$class <- "fwd_alt"
  fwd_deletions_alt <- fwd_deletions_alt %>% dplyr::rename(seq=ALTseq)
  
  fwd_insertions_ref <- subset(insertions, select=-ALTseq)
  fwd_insertions_ref$class <- "fwd_ref"
  fwd_insertions_ref <- fwd_insertions_ref %>% dplyr::rename(seq=REFseq)
  
  fwd_insertions_alt <- subset(insertions, select=-REFseq)
  fwd_insertions_alt$class <- "fwd_alt"
  fwd_insertions_alt <- fwd_insertions_alt %>% dplyr::rename(seq=ALTseq)
  
  fwd_snps_ref <- subset(snps, select=-ALTseq)
  fwd_snps_ref$class <- "fwd_ref"
  fwd_snps_ref <- fwd_snps_ref %>% dplyr::rename(seq=REFseq)
  
  fwd_snps_alt <- subset(snps, select=-REFseq)
  fwd_snps_alt$class <- "fwd_alt"
  fwd_snps_alt <- fwd_snps_alt %>% dplyr::rename(seq=ALTseq)
  
  
  rev_deletions_ref <- subset(rev_deletions, select=-ALTseq)
  rev_deletions_ref$class <- "rev_ref"
  rev_deletions_ref <- rev_deletions_ref %>% dplyr::rename(seq=REFseq)
  
  rev_deletions_alt <- subset(rev_deletions, select=-REFseq)
  rev_deletions_alt$class <- "rev_alt"
  rev_deletions_alt <- rev_deletions_alt %>% dplyr::rename(seq=ALTseq)
  
  rev_insertions_ref <- subset(rev_insertions, select=-ALTseq)
  rev_insertions_ref$class <- "rev_ref"
  rev_insertions_ref <- rev_insertions_ref %>% dplyr::rename(seq=REFseq)
  
  rev_insertions_alt <- subset(rev_insertions, select=-REFseq)
  rev_insertions_alt$class <- "rev_alt"
  rev_insertions_alt <- rev_insertions_alt %>% dplyr::rename(seq=ALTseq)
  
  rev_snps_ref <- subset(rev_snps, select=-ALTseq)
  rev_snps_ref$class <- "rev_ref"
  rev_snps_ref <- rev_snps_ref %>% dplyr::rename(seq=REFseq)
  
  rev_snps_alt <- subset(rev_snps, select=-REFseq)
  rev_snps_alt$class <- "rev_alt"
  rev_snps_alt <- rev_snps_alt %>% dplyr::rename(seq=ALTseq)
  
  complete_variants <- rbind(fwd_deletions_ref, fwd_deletions_alt, fwd_insertions_ref, fwd_insertions_alt, fwd_snps_ref, fwd_snps_alt, rev_deletions_ref, rev_deletions_alt, rev_insertions_ref, rev_insertions_alt, rev_snps_ref, rev_snps_alt)
  
  # 13) ASSEMBLE CONSTRUCTS ############################################################################################
  complete_variants <- complete_variants %>% dplyr::slice(rep(1:n(), each=25))
  
  fwdprimer <- "ACTGGCCGCTTCACTG"
  fwdspacer <- "TG"
  enzyme1 <- "GGTACC"
  enzyme2 <- "TCTAGA"
  revspacer <- "GGC"
  revprimer <- "AGATCGGAAGAGCGTCG"
  
  tags_needed = nrow(complete_variants)
  if(tags_needed>nrow(tags))
  {
    print("WARNING: The number of tags needed exceeds the number of tags available. Some oligos will not have tags")
  }
  tags_to_use <- tags[1:tags_needed,]
  complete_variants <- cbind(complete_variants, tags_to_use)
  
  complete_variants$construct <- paste0(fwdprimer, fwdspacer, complete_variants$seq, enzyme1, enzyme2, complete_variants$barcode, revspacer, revprimer)
  
  # 14) BIND SCRAMBLED WITH OTHER VARIANTS ############################################################################
  
  complete_variants <- complete_variants %>% dplyr::rename(tag=barcode)
  final_output <- rbind(scrambled_vars, complete_variants)
  # write.csv(final_output, "FINAL_OUTPUT_13Jan2021.csv", row.names=FALSE)
  
  # 15) FINAL CHECK FOR DIGESTION SITES (INCLUDING JUNCTIONS) #########################################################
  # will detect remaining digestion sites (e.g. there was more than 1 in a sequence) as well
  # This should be a very low number of sequences and as such can be fixed manually
  
  final_check <- subset(final_output, select=c(ID, class, construct))
  final_check$short_seq <- gsub('.{39}$', '', final_check$construct) #remove last 39 bases
  final_check <- checkDigest(final_check)
  print(final_check)

  
  write.csv(final_output, "OLIGO_LIBRARY.csv", row.names=FALSE)
  
  # # check if all tags are unique
  # test <- unique(final_output$tag)
  # print(length(test)) # should be equal to 239,800 in this case (total number of sequences generated)
  # 
  # #### CREATE FILES FOR SUBMISSIONS ################################################################################
  # 
  # # Create a unique identifier for each sequence
  # final_output_SUBMISSION <- final_output
  # final_output_SUBMISSION$row_number <- seq.int(nrow(final_output_SUBMISSION))
  # final_output_SUBMISSION$identifier <- paste0(final_output_SUBMISSION$ID,"_",final_output_SUBMISSION$class,"_",final_output_SUBMISSION$row_number)
  # final_output_SUBMISSION <- subset(final_output_SUBMISSION, select=-row_number)
  # 
  # print(final_output_SUBMISSION)
  # write.csv(final_output_SUBMISSION, "OLIGO_LIBRARY.csv", row.names=FALSE)
  # 
  # # Shuffle sequences
  # row_indices_shuffled <- sample(nrow(final_output_SUBMISSION))
  # shuffled_output <- final_output_SUBMISSION[row_indices_shuffled,]
  # 
  # write.csv(shuffled_output, "OLIGO_LIBRARY_SHUFFLED.CSV", row.names=FALSE)
}

