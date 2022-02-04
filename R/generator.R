# CUSTOM MPRA DESIGN GENERATOR #############################################################################################

#####################################################################################################################

#' Generate Design File
#'
#' To actually generate the file containing all oligos in the library, run this function.
#'
#' @param fwdprimer forward primer
#' @param revprimer reverse primer
#' @param tags_per_variant number of oligos with unique tags created for each variant
#' @param enz1 first digestion site
#' @param enz2 second digestion site
#' @param enz3 third digestion site
#' @param enz1FIX how the first digestion site should be changed if found in genomic region (use . to represent any base)
#' @param enz2FIX how the second digestion site should be changed if found in genomic region (use . to represent any base)
#' @param enz3FIX how the third digestion site should be changed if found in genomic region (use . to represent any base)
#' @param variant_input_path path to file with variant input (e.g. "Documents/input_files/variant_input.csv")
#' @param tag_path path to file containing tags
#' @param scrambled_path path to file containing scrambled sequences
#' @return design file
#' @export
generate = function(fwdprimer, revprimer, tags_per_variant, enz1, enz2, enz3, enz1FIX, enz2FIX, enz3FIX, variant_input_path, tag_path, scrambled_path){
  # 1) LOAD FILTERED TAGS #############################################################################################
  tags <- readr::read_csv(tag_path, col_names=TRUE, col_types=cols("c"))

  # 2) CHECK TAGS FOR DIGESTION SITES #################################################################################
  tag_check <- tags
  tag_check$ez1 <- grepl(enz1, tag_check$barcode)
  tag_check$ez2 <- grepl(enz2, tag_check$barcode)
  tag_check$ez3 <- grepl(enz3, tag_check$barcode)
  bad_tag_rows <- which(tag_check$ez1==TRUE|tag_check$ez2==TRUE|tag_check$ez3==TRUE)
  if(length(bad_tag_rows)!=0)
  {
    tags <- tags[-c(bad_tag_rows),]
  }
  # this section should be tested before being used as our tag set had no digestion sites so no tags were removed

  # 2) CREATE SCRAMBLED SEQUENCES #####################################################################################
  if(!missing(scrambled_path))
  {
    scrambled_vars <- readr::read_csv(scrambled_path, col_names=TRUE, col_types=cols("c","c"))
    if(nrow(tags) < nrow(scrambled_vars)*tags_per_variant)
    {
      print("ERROR: The number of tags provided is not sufficient to create all oligos containing scrambled sequences. Provide more tags or change the number of tags used per variant.")
      print("Design File generation halted")
      return(FALSE)
    }
    else
    {
      # separate tags to be used for scrambled sequences
      tag_number <- nrow(tags)
      scrambled_tags_needed <- nrow(scrambled_vars)*tags_per_variant
      scrambled_tags <- tags[1:scrambled_tags_needed,]
      start_of_next_tags <- scrambled_tags_needed + 1
      tags <- tags[start_of_next_tags:tag_number,]

      #fwdprimer <- "ACTGGCCGCTTCACTG"
      fwdspacer <- "TG"
      enzyme1 <- enz1
      enzyme2 <- enz2
      revspacer <- "GGC"
      #revprimer <- "AGATCGGAAGAGCGTCG"

      scrambled_vars <- scrambled_vars %>% dplyr::slice(rep(1:n(), each=tags_per_variant))

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
    }
  }

  # 2) READ IN VARIANT INPUT ##########################################################################################
  all_variants <- readr::read_csv(variant_input_path, col_names=TRUE, col_types=cols("c","d","c","c","c"))

  # 3) GET SEQUENCES ##################################################################################################
  genome <- BSgenome.Hsapiens.UCSC.hg38
  all_variants$REFseq <- retrieveSeq(genome, all_variants$CHROM, all_variants$POS)

  # 3.5) GENERATE A LIST OF SEQUENCES THAT HAVE DIGESTION SITES #######################################################
  seqs_w_digest <- all_variants

  seqs_w_digest$Enz1 <- grepl(enz1, seqs_w_digest$REFseq)
  seqs_w_digest$Enz2 <- grepl(enz2, seqs_w_digest$REFseq)
  seqs_w_digest$Enz3 <- grepl(enz3, seqs_w_digest$REFseq)
  seqs_w_digest <- dplyr::filter(seqs_w_digest, Enz1=="TRUE"|Enz2=="TRUE"|Enz3=="TRUE")
  seqs_w_digest <- subset(seqs_w_digest, select=c(ID, REFseq, Enz1, Enz2, Enz3))
  write.csv(seqs_w_digest, "sequences_w_digestion_sites.csv", row.names=FALSE)

  # need a function to check that digestion site isn't around SNP (calculate how far away first or last base of digestion site is from middle)
  seqs_w_digest_coords_e1 <- dplyr::filter(seqs_w_digest, Enz1=="TRUE")
  seqs_w_digest_coords_e2 <- dplyr::filter(seqs_w_digest, Enz2=="TRUE")
  seqs_w_digest_coords_e3 <- dplyr::filter(seqs_w_digest, Enz3=="TRUE")

  ### find e1 variants that are unfixable
  seqs_w_digest_coords_e1$num_occurences <- stringr::str_count(seqs_w_digest_coords_e1$REFseq, enz1)
  seq_w_multiple_e1 <- filter(seqs_w_digest_coords_e1, num_occurences>1)
  seqs_w_digest_coords_e1 <- filter(seqs_w_digest_coords_e1, num_occurences==1)

  seqs_w_digest_coords_e1$coords <- str_locate(seqs_w_digest_coords_e1$REFseq, enz1)
  seqs_w_digest_coords_e1$c1 <- seqs_w_digest_coords_e1$coords[,"start"]
  seqs_w_digest_coords_e1$c2 <- seqs_w_digest_coords_e1$coords[,"end"]

  cant_fix_e1 <- filter(seqs_w_digest_coords_e1, c1<=73 & c2>=73) # To make length adjustable, calculate middle base from user input and use here instead of 73
  cant_fix_e1 <- subset(cant_fix_e1, select=-c(coords, c1, c2))

  multiple_unfixable_e1 <- seq_w_multiple_e1
  if(nrow(seq_w_multiple_e1)>0)
  {
    seq_w_multiple_e1$fixable <- sapply(seq_w_multiple_e1$REFseq, are_repeat_sites_fixable, enzx=enz1) # adjust length here as well
    multiple_unfixable_e1 <- filter(seq_w_multiple_e1, fixable=="FALSE")
    multiple_unfixable_e1 <- subset(multiple_unfixable_e1, select=-fixable)
    multiple_fixable_e1 <- filter(seq_w_multiple_e1, fixable=="TRUE")
  }

  cant_fix_e1 <- rbind(cant_fix_e1, multiple_unfixable_e1)

  ### find e2 variants that are unfixable
  seqs_w_digest_coords_e2$num_occurences <- stringr::str_count(seqs_w_digest_coords_e2$REFseq, enz2)
  seq_w_multiple_e2 <- filter(seqs_w_digest_coords_e2, num_occurences>1)
  seqs_w_digest_coords_e2 <- filter(seqs_w_digest_coords_e2, num_occurences==1)

  seqs_w_digest_coords_e2$coords <- str_locate(seqs_w_digest_coords_e2$REFseq, enz2)
  seqs_w_digest_coords_e2$c1 <- seqs_w_digest_coords_e2$coords[,"start"]
  seqs_w_digest_coords_e2$c2 <- seqs_w_digest_coords_e2$coords[,"end"]

  cant_fix_e2 <- filter(seqs_w_digest_coords_e2, c1<=73 & c2>=73) # To make length adjustable, calculate middle base from user input and use here instead of 73
  cant_fix_e2 <- subset(cant_fix_e2, select=-c(coords, c1, c2))

  multiple_unfixable_e2 <- seq_w_multiple_e2
  if(nrow(seq_w_multiple_e2)>0)
  {
    seq_w_multiple_e2$fixable <- sapply(seq_w_multiple_e2$REFseq, are_repeat_sites_fixable, enzx=enz2) # adjust length here as well
    multiple_unfixable_e2 <- filter(seq_w_multiple_e2, fixable=="FALSE")
    multiple_unfixable_e2 <- subset(multiple_unfixable_e2, select=-fixable)
    multiple_fixable_e2 <- filter(seq_w_multiple_e2, fixable=="TRUE")
  }

  cant_fix_e2 <- rbind(cant_fix_e2, multiple_unfixable_e2)

  ### find e3 variants that are unfixable
  seqs_w_digest_coords_e3$num_occurences <- stringr::str_count(seqs_w_digest_coords_e3$REFseq, enz3)
  seq_w_multiple_e3 <- filter(seqs_w_digest_coords_e3, num_occurences>1)
  seqs_w_digest_coords_e3 <- filter(seqs_w_digest_coords_e3, num_occurences==1)

  seqs_w_digest_coords_e3$coords <- str_locate(seqs_w_digest_coords_e3$REFseq, enz3)
  seqs_w_digest_coords_e3$c1 <- seqs_w_digest_coords_e3$coords[,"start"]
  seqs_w_digest_coords_e3$c2 <- seqs_w_digest_coords_e3$coords[,"end"]

  cant_fix_e3 <- filter(seqs_w_digest_coords_e3, c1<=73 & c2>=73) # To make length adjustable, calculate middle base from user input and use here instead of 73
  cant_fix_e3 <- subset(cant_fix_e3, select=-c(coords, c1, c2))

  multiple_unfixable_e3 <- seq_w_multiple_e3
  if(nrow(seq_w_multiple_e3)>0)
  {
    seq_w_multiple_e3$fixable <- sapply(seq_w_multiple_e3$REFseq, are_repeat_sites_fixable, enzx=enz3) # adjust length here as well
    multiple_unfixable_e3 <- filter(seq_w_multiple_e3, fixable=="FALSE")
    multiple_unfixable_e3 <- subset(multiple_unfixable_e3, select=-fixable)
    multiple_fixable_e3 <- filter(seq_w_multiple_e3, fixable=="TRUE")
  }

  cant_fix_e3 <- rbind(cant_fix_e3, multiple_unfixable_e3)

  ### remove unfixable variants from all variants
  cant_fix_variants <- rbind(cant_fix_e1, cant_fix_e2, cant_fix_e3)

  if(nrow(cant_fix_variants)>0)
  {
    print("WARNING: Some variants have a digestion site around the SNP. No oligos will be generated for these variants. You can see which variants have been excluded in the file titled 'cannot_be_fixed.csv.' ")
    write.csv(cant_fix_variants, "cannot_be_fixed.csv", row.names=FALSE)

    rows_with_unfixables <- which(all_variants$ID %in% cant_fix_variants$ID)
    all_variants <- all_variants[-rows_with_unfixables,]
  }


  # 4) DETECT AND REPAIR DIGESTION SITES ##############################################################################
  # By searching for & repairing sites here, we reduce the number of operations the program must perform.
  # This is acceptable because we have already confirmed that the variants in all_variants do not have digestion sites
  # at the SNP.


  # Assuming enzyme1 and enzyme2 do not include special characters
  all_variants$REFseq <- gsub(enz1, enz1FIX, all_variants$REFseq)
  all_variants$REFseq <- gsub(enz2, enz2FIX, all_variants$REFseq)

  # Assuming enzyme3 may or may not contain special characters (put in loop so all sites are repaired)
  all_variants$enz3 <- grepl(enz3, all_variants$REFseq)
  variants_wo_enz3 <- dplyr::filter(all_variants, enz3=="FALSE")
  variants_w_enz3 <- dplyr::filter(all_variants, enz3=="TRUE")
  variants_w_enz3_fixed <- variants_w_enz3

  while("TRUE" %in% variants_w_enz3_fixed$enz3 == TRUE)
  {
    variants_w_enz3 <- dplyr::filter(variants_w_enz3_fixed, enz3=="TRUE")
    variants_w_enz3_fixed <- variants_w_enz3
    variants_w_enz3$coords <- stringr::str_locate(variants_w_enz3$REFseq, enz3)
    variants_w_enz3$coord1 <- variants_w_enz3$coords[,"start"]
    variants_w_enz3$coord2 <- variants_w_enz3$coords[,"end"]
    variants_w_enz3_fixed$REFseq <- fixDigWild(enz3, enz3FIX, variants_w_enz3$REFseq, variants_w_enz3$coord1)
    variants_w_enz3_fixed$enz3 <- grepl(enz3, variants_w_enz3_fixed$REFseq)
    completely_fixed_transfer <- dplyr::filter(variants_w_enz3_fixed, enz3=="FALSE")
    variants_wo_enz3 <- rbind(variants_wo_enz3, completely_fixed_transfer)
  }

  all_variants <- variants_wo_enz3
  all_variants <- subset(all_variants, select=-enz3)

  # 5) SEPARATE SNPS, INSERTIONS, AND DELETIONS #######################################################################
  insertions <- dplyr::filter(all_variants, REF=="-")
  deletions <- dplyr::filter(all_variants, ALT=="-")
  snps <- dplyr::filter(all_variants, REF!="-" & ALT!="-")

  # #) COMPARE WHAT WILL BE CHANGED TO WHAT SHOULD BE CHANGED #########################################################
  # a way of checking that the coordinates are correct

  # check deletions
  deletion_coord_check <- deletions
  deletion_coord_check$deleted_bases <- delBases(deletion_coord_check$REFseq, deletions$REF) # update this if updating length
  deletion_coord_check$bases_match <- ifelse(deletion_coord_check$REF==deletion_coord_check$deleted_bases, "YES", "NO")
  wrong_coords <- filter(deletion_coord_check, bases_match=="NO")
  if(nrow(wrong_coords)>0)
  {
    del_err_msg <- "The bases that will be deleted by MPRADesignGenerator (REF) were compared to the bases you have indicated should be deleted. For some variants, these bases do not match.\nThis usually indicates that the coordinates you have provided are incorrect (often the POS coordinate is off by 1). The file 'deletion_mismatch.csv' contains a list of the\nvariants where the bases marked for deletion do not match."
    writeLines(del_err_msg)
    write.csv(wrong_coords, "deletion_mismatch.csv")
  }

  # check snps
  snp_coord_check <- snps
  snp_coord_check$changed_bases <- snpBases(snp_coord_check$REFseq) # update this if updating length
  snp_coord_check$bases_match <- ifelse(snp_coord_check$REF==snp_coord_check$changed_bases, "YES", "NO")
  snp_wrong_coords <- filter(snp_coord_check, bases_match=="NO")
  if(nrow(snp_wrong_coords)>0)
  {
    snp_err_msg <- "The base that will be changed by MPRADesignGenerator (REF) was compared to the base you have indicated should be changed. For some variants, these bases do not match.\nThis usually indicates that the coordinates you have provided are incorrect (often the POS coordinate is off by 1) or that REF and ALT bases are switched.\nThe file 'snp_mismatch.csv' contains a list of the variants where the bases marked for change do not match."
    writeLines(snp_err_msg)
    write.csv(snp_wrong_coords, "snp_mismatch.csv")
  }

  # no good way to check for insertions at the moment

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

  # 10) GENERATE REV COMP SEQ FOR INSERTIONS ###########################################################################
  rev_insertions <- insertions

  rev_insertions_ref_DNAset <- Biostrings::DNAStringSet(rev_insertions$REFseq)
  rev_insertions_ref_DNAset <- Biostrings::reverseComplement(rev_insertions_ref_DNAset)
  rev_insertions_ref_DNAset <- data.frame(REFseq=as.character(rev_insertions_ref_DNAset))

  rev_insertions_alt_DNAset <- Biostrings::DNAStringSet(rev_insertions$ALTseq)
  rev_insertions_alt_DNAset <- Biostrings::reverseComplement(rev_insertions_alt_DNAset)
  rev_insertions_alt_DNAset <- data.frame(ALTseq=as.character(rev_insertions_alt_DNAset))

  rev_insertions$REFseq <- rev_insertions_ref_DNAset$REFseq
  rev_insertions$ALTseq <- rev_insertions_alt_DNAset$ALTseq

  # 11) GENERATE REV COMP SEQ FOR SNPS ##################################################################################
  rev_snps <- snps

  rev_snps_ref_DNAset <- Biostrings::DNAStringSet(rev_snps$REFseq)
  rev_snps_ref_DNAset <- Biostrings::reverseComplement(rev_snps_ref_DNAset)
  rev_snps_ref_DNAset <- data.frame(REFseq=as.character(rev_snps_ref_DNAset))

  rev_snps_alt_DNAset <- Biostrings::DNAStringSet(rev_snps$ALTseq)
  rev_snps_alt_DNAset <- Biostrings::reverseComplement(rev_snps_alt_DNAset)
  rev_snps_alt_DNAset <- data.frame(ALTseq=as.character(rev_snps_alt_DNAset))

  rev_snps$REFseq <- rev_snps_ref_DNAset$REFseq
  rev_snps$ALTseq <- rev_snps_alt_DNAset$ALTseq

  # 12) BIND ALL DATA TOGETHER ###########################################################################################
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
  complete_variants <- complete_variants %>% dplyr::slice(rep(1:n(), each=tags_per_variant))

  #fwdprimer <- "ACTGGCCGCTTCACTG"
  fwdspacer <- "TG"
  enzyme1 <- enz1
  enzyme2 <- enz2
  revspacer <- "GGC"
  #revprimer <- "AGATCGGAAGAGCGTCG"

  tags_needed = nrow(complete_variants)
  if(tags_needed>nrow(tags))
  {
    print("ERROR: The number of tags needed to create all oligos exceeds the number of tags available. Provide more tags or change the number of tags used per variant.")
    print("Design File generation halted")
    return(FALSE)
  }

  tags_to_use <- tags[1:tags_needed,]
  complete_variants <- cbind(complete_variants, tags_to_use)

  complete_variants$construct <- paste0(fwdprimer, fwdspacer, complete_variants$seq, enzyme1, enzyme2, complete_variants$barcode, revspacer, revprimer)

  # 14) BIND SCRAMBLED WITH OTHER VARIANTS ############################################################################

  complete_variants <- complete_variants %>% dplyr::rename(tag=barcode)
  if(missing(scrambled_path))
  {
    final_output <- complete_variants
  }
  else
  {
    final_output <- rbind(scrambled_vars, complete_variants)
  }

  # 15) FINAL CHECK FOR DIGESTION SITES (INCLUDING JUNCTIONS) #########################################################
  # will detect remaining digestion sites (though there should be none) as well as at junctions where components of the construct meet.
  # This should be a very low number of sequences and as such can be fixed manually

  final_check <- subset(final_output, select=c(ID, class, construct))
  final_check$short_seq <- gsub('.{39}$', '', final_check$construct) #remove last 39 bases (ASSUMES REVPRIMER IS 17BP)
  final_check <- checkDigest(final_check, enz1, enz2, enz3) #update to check for variable site

  if(nrow(final_check)>0)
  {
    print("Some oligos have a digestion site where different components (primer, tag, etc.) are joined together. The file titled oligos_w_digestion_site.csv contains a complete list of such oligos.")
    write.csv(final_check, "oligos_w_digestion_site.csv", row.names=FALSE)
  }

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
