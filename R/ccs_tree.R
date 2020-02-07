# --------------------------------------------------------------------------------- #
# --------------- Code for creating tree of outcomes for CCS codes ---------------- #
# --------------------------------------------------------------------------------- #

# group = character string specifying that only codes beginning with group 
# should be considered
# For example, group = "7" means only codes starting with 7 will be included on the tree
# This is the subtree corresponding to cardiovascular disease
# The default is NULL (returns full tree of CCS codes)

ccs_tree <- function(group = NULL) {
  # Get data.frame showing mapping from ICD9 to multilevel CCS
  ccs_icd <- data.table(icd = unlist(icd9_map_multi_ccs[[1]]))
  for (i in 1:4) {
    ccs_list <- icd9_map_multi_ccs[[i]]
    ccs_df <- data.table(icd = unlist(ccs_list))
    ccs_df[ , ccs := ccs_list %>% 
      names %>% # names of the list entries are the CCS codes
      sapply(FUN = function(nm) rep(nm, length(ccs_list[[nm]]))) %>%
      unlist]
    setnames(ccs_df, "ccs", paste0("l", i))
    ccs_icd <- merge(ccs_icd, ccs_df, by = "icd", all.x = T, all.y = F)
  }
  
  # CCS codes only
  ccs <- ccs_icd
  ccs$icd <- NULL
  ccs <- ccs[!duplicated(ccs), ]
  
  # Order CCS codes appropriately
  ccs_levels <- matrix(nrow = nrow(ccs), ncol = 4)
  for (i in 1:nrow(ccs_levels)) {
    splt <- as.integer(str_split(ccs$l4[i], "\\.")[[1]])
    if (length(splt) != 4) {
      splt <- as.integer(str_split(ccs$l3[i], "\\.")[[1]])
      if (length(splt) != 3) {
        splt <- as.integer(str_split(ccs$l2[i], "\\.")[[1]])
        if (length(splt) != 2) {
          splt <- as.integer(ccs$l1[i])
        }
      }
      splt <- c(splt, rep(0, 4 - length(splt)))
    }
    ccs_levels[i, ] <- splt
  }
  ccs_levels <- as.data.table(ccs_levels)
  names(ccs_levels) <- c("l1", "l2", "l3", "l4")
  ccs <- ccs[order(ccs_levels$l1, ccs_levels$l2, ccs_levels$l3, ccs_levels$l4), ]
  
  # Keep only diseases starting with group
  if (!is.null(group)) {
    group_length <- length(str_split(group, "\\.")[[1]])
    group_lvl <- paste0("l", group_length)
    ccs <- ccs[get(group_lvl) == group]
  }
  if (group_length > 1) {
    ccs <- ccs[ , .SD, 
          .SDcols = sapply(group_length:4, function(x) paste0("l", x))]
  }
  
  # Make tree
  ccs_mat <- as.matrix(ccs)
  edges <- ccs_mat[ , c(1, 2)]
  if (ncol(ccs_mat) > 2) {
    for (i in 3:ncol(ccs_mat)) {
      edges <- rbind(edges, ccs_mat[ , c(i - 1, i)])
    }
  }
  edges <- edges[edges[ , 2] != " ", ]
  edges <- edges[!duplicated(edges), ]
  tr <- igraph::graph_from_edgelist(e = edges, directed = T)
  leaves <- igraph::V(tr)[igraph::degree(tr, mode = "out") == 0]
  igraph::V(tr)$leaf <- FALSE
  igraph::V(tr)$leaf[V(tr) %in% leaves] <- TRUE
  
  # ICD mapping
  ccs$keep <- T
  ccs_icd <- merge(ccs_icd, ccs, by = names(ccs)[-ncol(ccs)],
                   all.x = T, sort = FALSE)
  ccs_icd <- ccs_icd[(keep)]
  ccs_ics <- ccs_icd[ , keep := NULL]
  for (i in 2:4) {
    nm <- paste0("l", i)
    nm_p <- paste0("l", i - 1)
    ccs_icd[get(nm) == " ", (nm) := .SD, .SDcols = nm_p]
  }
  ccs_icd <- ccs_ics[ , c("icd", "l4")]
  setnames(ccs_icd, "l4", "ccs")

  return(list(tr = tr, ccs_icd_mapping = ccs_icd))
}