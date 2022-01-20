#' RETRIEVE SEQUENCES
#'
#' Retrieve a 145bp sequence from BSgenome.Hsapiens.UCSC.hg19
#' (72bp downstream from specified base, 72b upstream from specified base)
#'
#' @param genome specification, chromosome number, coordinate
#' @return 145bp sequence
#' @export
retrieveSeq = function(genome, chromNum, coord)
{
  start = coord - 72
  end = coord + 72
  chr = paste0("chr", chromNum)
  sequence <- as.character(Biostrings::getSeq(genome, chr, start, end))
  return(sequence)
}

#' FIX SfiI DIGESTION SITE
#'
#' Changes SfiI site (GGCCNNNNNGGCC -> GCGCNNNNNGGCC)
#' Will only first first digestion site in sequence.
#'
#'
#' @param 145bp sequence containing site, coordinate of first base of
#' restriction site, coordinate of second base of restriction site.
#' @return corrected sequence
#' @export
fixSfiI = function(sequence, coord1, coord2)
{
  temp <- substr(sequence, coord1, coord2)
  substr(temp, 2, 2) <- 'C'
  substr(temp, 3, 3) <- 'G'
  str_sub(sequence, coord1, coord2) <- temp
  return(sequence)
}

#' CREATE ALTERNATIVE SEQUENCES FOR DELETIONS
#'
#' Delete specified number of bases from ref seq. (Counts the number of bases
#' in deletion string, removes bases starting at 73)
#'
#'
#' @param 145bp reference sequence, string containing bases to delete.
#' @return alt sequence with specified bases removed (counts the number of bases
#' in string, removes bases starting at 73)
#' @export
altDel = function(reference_sequence, deletion_bases)
{
  num_bases_to_delete = nchar(deletion_bases)
  str_sub(reference_sequence, 73, 73+num_bases_to_delete-1) <- ""
  return(reference_sequence)
}

#' CREATE ALTERNATIVE SEQUENCES FOR INSERTIONS
#'
#' Insert specified bases after base 73 into provided ref seq. Shortens sequence
#' from both sides so total length is <= 145bp.
#'
#'
#' @param 145bp reference sequence, string containing bases to insert.
#' @return 145bp alt sequence with specified bases inserted
#' @export
altIns = function(reference_sequence, insertion_bases)
{
  part_one <- substr(reference_sequence, 1, 73)
  part_two <- substr(reference_sequence, 74, 145)
  alt_sequence <- paste0(part_one, insertion_bases, part_two)
  # shorten sequences to 145 bases or less
  to_remove <- ceiling((nchar(alt_sequence)-145)/2) 
  alt_sequence <- substr(alt_sequence, to_remove+1, nchar(alt_sequence)-to_remove)
  return(alt_sequence)
}

#' SHORTEN REFERENCE SEQUENCES FOR INSERTIONS
#'
#' Ensures ref seq for insertions has the same bases removed
#' from the ends as alt seq.
#'
#'
#' @param 145bp reference sequence, string containing bases that were inserted.
#' @return ref sequence with appropriate end bases removed
#' @export
refIns = function(reference_sequence, insertion_bases)
{
  insertion_length <- nchar(insertion_bases)
  to_remove <- ceiling(insertion_length/2)
  reference_sequence <- substr(reference_sequence, to_remove+1, 145-to_remove)
  return(reference_sequence)
}

#' CREATE ALTERNATIVE SEQUENCES FOR SNPS
#'
#' Change base 73 to alt base.
#'
#'
#' @param 145bp reference sequence, string containing base to swap in.
#' @return alt sequence with specified base swapped in
#' @export
altSnps = function(reference_sequence, new_base)
{
  str_sub(reference_sequence, 73, 73) <- new_base
  return(reference_sequence)
}

#' CHECK FOR ALL THREE DIGESTION SITES
#'
#' Check all sequences in a dataframe for digestion sites
#' (SfiI, KpnI, XbaI)
#' (column with sequences tested should be named short_seq)
#'
#'
#' @param dataframe with sequences to test 
#' @return new dataframe containing sequences which had a digestion site
#' @export
checkDigest = function(df_to_test)
{
  ez1 <- "GGTACC"
  ez2 <- "TCTAGA"
  ez3 <- "GGCC.....GGCC"
  df_to_test$enz1 <- grepl(ez1, df_to_test$short_seq)
  df_to_test$enz2 <- grepl(ez2, df_to_test$short_seq)
  df_to_test$enz3 <- grepl(ez3, df_to_test$short_seq)
  passed_w_digestion <- tidyverse::filter(df_to_test, enz1=='TRUE' | enz2=='TRUE' | enz3=='TRUE')
  return(passed_w_digestion)
}

####
#' GET BASES TO BE DELETED
#'
#' Check which bases the program will delete. Can be used to compare to
#' the bases that should be deleted per variant input file.
#'
#'
#' @param reference sequence, string containing bases to delete 
#' @return string containing bases that the program will delete.
#' @export
delBases = function(reference_sequence, deletion_bases)
{
  num_bases_to_delete = nchar(deletion_bases)
  bases_deleted <- substr(reference_sequence, 73, 73+num_bases_to_delete-1)
  return(bases_deleted)
}

#' GET BASE TO BE CHANGED
#'
#' Check which base the program will change. Can be used to compare to
#' the bases that should be changed per variant input file.
#'
#'
#' @param reference sequence 
#' @return string containing base that the program will change.
#' @export
snpBases = function(reference_sequence)
{
  old_base <- substr(reference_sequence, 73, 73)
  return(old_base)
}

# GET BASE AFTER WHICH INSERTION WILL OCCUR?
##