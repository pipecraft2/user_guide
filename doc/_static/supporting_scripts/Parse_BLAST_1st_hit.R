#!/usr/bin/env Rscript
### Parse BLAST results and discard non-target OTUs (if any)

### Specify input file
 # BLAST 1st hit output from PipeCraft
blast_1st_hit_file = "BLAST_1st_best_hit.txt" 

### Specify target group(s) (if any)
 # Target taxonomic group(s) to keep
target = c("Fungi")
 # Taxonomic level to filter on: 
 #           Kingdom | Phylum | Class | Order | Family | Genus | Species
tax_level = "Kingdom"

### Specify sim_score thresholds for taxonomic levels
 # Minimum sim_score required for reliable assignment at each level
sim_score_thresholds = list(
  Class = 75,
  Order = 80,
  Family = 85,
  Genus = 90,
  Species = 97
)

### Specify e-value thresholds
 # e-value < e-50: reliable for kingdom assignment
 # e-value > e-20: mark as "unknown"
evalue_reliable = 1e-50  # Reliable threshold for kingdom assignment
evalue_unknown = 1e-20   # Unknown threshold, mark OTU as "unknown"  
#--------------------------------------#
#--------------------------------------#

# Load blast_1st_hit_file
blast_1st_hit = read.table("BLAST_1st_best_hit.txt", 
                           header = TRUE, sep = "+", fill = TRUE)

###################################
### Parse 1st best hit taxonomy ###
###################################
### Format: Accession;k__Kingdom;p__Phylum;
#           c__Class;o__Order;f__Family;g__Genus;s__Species
parse_taxonomy = function(taxonomy_string) {
  # Initialize result with NAs
  result = data.frame(
    Accession = NA, Kingdom = NA, Phylum = NA, Class = NA, 
    Order = NA, Family = NA, Genus = NA, Species = NA,
    stringsAsFactors = FALSE
  )
  
  # Handle NA/empty cases
  if (is.na(taxonomy_string) || 
      taxonomy_string == "" || 
      taxonomy_string == "*") {
    return(result)
  }
  
  # Split by semicolon
  ranks = strsplit(taxonomy_string, ";")[[1]]
  
  # Extract accession number (first field)
  if (length(ranks) > 0) {
    result$Accession = ranks[1]
  }
  
  # Parse each rank (format: prefix__taxon_name)
  for (rank in ranks) {
    if (grepl("^k__", rank)) {
      # Kingdom
      result$Kingdom = sub("^k__", "", rank)
    } else if (grepl("^p__", rank)) {
      # Phylum
      result$Phylum = sub("^p__", "", rank)
    } else if (grepl("^c__", rank)) {
      # Class
      result$Class = sub("^c__", "", rank)
    } else if (grepl("^o__", rank)) {
      # Order
      result$Order = sub("^o__", "", rank)
    } else if (grepl("^f__", rank)) {
      # Family
      result$Family = sub("^f__", "", rank)
    } else if (grepl("^g__", rank)) {
      # Genus
      result$Genus = sub("^g__", "", rank)
    } else if (grepl("^s__", rank)) {
      # Species: combine with genus as Genus_species
      species_epithet = sub("^s__", "", rank)
      if (!is.na(result$Genus) && result$Genus != "") {
        result$Species = paste0(result$Genus, "_", species_epithet)
      } else {
        result$Species = species_epithet
      }
    }
  }
  
  return(result)
}

# Apply parsing function 
taxonomy_parsed = do.call(rbind, 
                          lapply(blast_1st_hit$X1st_hit, parse_taxonomy))

# Compile parsed 1st hit file
blast_taxonomy = cbind(
  qseqid = blast_1st_hit$qseqid,
  query_seq = blast_1st_hit$query_seq,
  taxonomy_parsed,
  qlen = blast_1st_hit$qlen,
  evalue = blast_1st_hit$evalue,
  nident = blast_1st_hit$nident,
  mismatch = blast_1st_hit$mismatch,
  qcovs = blast_1st_hit$qcovs,
  adj_qcov = blast_1st_hit$adj_qcov,
  pident = blast_1st_hit$pident,
  sim_score = blast_1st_hit$sim_score)

##############################################
### Apply e-value and sim_score thresholds ###
##############################################
# Mark as "unknown" if e-value > e-20
unknown_mask = blast_taxonomy$evalue > evalue_unknown
blast_taxonomy$Kingdom[unknown_mask] = "unknown"
blast_taxonomy$Phylum[unknown_mask] = NA
blast_taxonomy$Class[unknown_mask] = NA
blast_taxonomy$Order[unknown_mask] = NA
blast_taxonomy$Family[unknown_mask] = NA
blast_taxonomy$Genus[unknown_mask] = NA
blast_taxonomy$Species[unknown_mask] = NA

# Apply sim_score thresholds for valid hits (e-value < e-20)
valid_mask = blast_taxonomy$evalue < evalue_unknown

# Class threshold
class_below = valid_mask & (!is.na(blast_taxonomy$sim_score) & 
                              blast_taxonomy$sim_score < 
                              sim_score_thresholds$Class)
blast_taxonomy$Class[class_below] = NA
blast_taxonomy$Order[class_below] = NA
blast_taxonomy$Family[class_below] = NA
blast_taxonomy$Genus[class_below] = NA
blast_taxonomy$Species[class_below] = NA

# Order threshold
order_below = valid_mask & (!is.na(blast_taxonomy$sim_score) & 
                              blast_taxonomy$sim_score < 
                              sim_score_thresholds$Order) & 
  !is.na(blast_taxonomy$Class)
blast_taxonomy$Order[order_below] = NA
blast_taxonomy$Family[order_below] = NA
blast_taxonomy$Genus[order_below] = NA
blast_taxonomy$Species[order_below] = NA

# Family threshold
family_below = valid_mask & (!is.na(blast_taxonomy$sim_score) & 
                               blast_taxonomy$sim_score < 
                               sim_score_thresholds$Family) & 
  !is.na(blast_taxonomy$Order)
blast_taxonomy$Family[family_below] = NA
blast_taxonomy$Genus[family_below] = NA
blast_taxonomy$Species[family_below] = NA

# Genus threshold
genus_below = valid_mask & (!is.na(blast_taxonomy$sim_score) & 
                              blast_taxonomy$sim_score < 
                              sim_score_thresholds$Genus) & 
  !is.na(blast_taxonomy$Family)
blast_taxonomy$Genus[genus_below] = NA
blast_taxonomy$Species[genus_below] = NA

# Species threshold
species_below = valid_mask & (!is.na(blast_taxonomy$sim_score) & 
                                blast_taxonomy$sim_score < 
                                sim_score_thresholds$Species) & 
  !is.na(blast_taxonomy$Genus)
blast_taxonomy$Species[species_below] = NA

#######################################
### Filter to target group (if any) ###
#######################################
# Discard OTUs that are not classified as target at specified taxonomic level
if (length(target) > 0 && !all(is.na(target)) && target[1] != "") {
  if (tax_level %in% colnames(blast_taxonomy)) {
    n_before = nrow(blast_taxonomy)
    blast_taxonomy = blast_taxonomy[blast_taxonomy[[tax_level]] %in% 
                                       target, , drop = FALSE]
    n_after = nrow(blast_taxonomy)
    n_excluded = n_before - n_after
    cat("\n Filtered to", tax_level, "level:", 
        paste(target, collapse = ", "), "\n")
    cat(" Excluded", n_excluded, "OTUs not matching target.\n")
    cat(" Remaining OTUs:", n_after, "\n")
  } else {
    warning(paste("Column '", tax_level, 
                  "' not found in taxonomy table; no filtering applied.", 
                  sep = ""))
  }
}

# Write output
write.table(blast_taxonomy, file = "BLAST_1st_hit_parsed.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)



###################################
### Consensus from 10 best hits ###
###################################
# Description of the method:
# Consensus taxonomy is derived from up to 10 BLAST hits per OTU using the following steps:
#
# 1. E-value filtering:
#    - Hits with e-value > e-20 are excluded (marked as "unknown" if all hits exceed this)
#    - For Kingdom/Phylum: only hits with e-value < e-50 are considered (reliable threshold)
#    - For lower ranks (Class, Order, Family, Genus, Species): hits with e-value < e-20 are used
#
# 2. Sim_score quality filtering:
#    - consensus is built from high-quality hits only
#    - Only "good" hits are used: hits with sim_score within 5 points of the best hit
#
# 3. Rank-specific consensus (processed from Species → Genus → Family → Order → Class):
#    - For each rank, only hits meeting the sim_score threshold are considered:
#      * Species: sim_score >= 97
#      * Genus: sim_score >= 90
#      * Family: sim_score >= 85
#      * Order: sim_score >= 80
#      * Class: sim_score >= 75
#    - Frequency counting: counts how many hits assign each taxonomic value
#    - Score comparison: gets the best sim_score for each unique taxonomic value
#    - Consensus decision:
#      * If majority value has best score, assign it
#      * If minority value has score >2 points higher than majority, don't assign this rank
#        (check if all hits agree at a higher taxonomic rank instead)
#    - All ranks that can be assigned are assigned (not just the most specific)
#
# 4. Hierarchical assignment:
#    - If a rank cannot be assigned, all more specific ranks are set to NA
#    - If Species/Genus cannot be assigned, checks if all hits agree at Family/Order/Class level
#
# 5. Output:
#    - Returns consensus taxonomy with all ranks that meet their thresholds


# Function to get consensus taxonomy from multiple hits
get_consensus_taxonomy = function(hit_taxonomies, hit_scores, hit_evalues) {
  # hit_taxonomies: list of parsed taxonomy data.frames
  # hit_scores: vector of sim_score values for each hit
  # hit_evalues: vector of e-value values for each hit
  
  if (length(hit_taxonomies) == 0) {
    return(data.frame(
      Accession = NA, Kingdom = NA, Phylum = NA, Class = NA, 
      Order = NA, Family = NA, Genus = NA, Species = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  # Check e-values: if all hits have e-value > e-20, mark as "unknown"
  if (all(hit_evalues > evalue_unknown, na.rm = TRUE)) {
    return(data.frame(
      Accession = NA, Kingdom = "unknown", Phylum = NA, Class = NA, 
      Order = NA, Family = NA, Genus = NA, Species = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  # Filter hits: only use hits with e-value < e-20 for consensus
  valid_hits = hit_evalues < evalue_unknown
  if (sum(valid_hits, na.rm = TRUE) == 0) {
    return(data.frame(
      Accession = NA, Kingdom = "unknown", Phylum = NA, Class = NA, 
      Order = NA, Family = NA, Genus = NA, Species = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  # Separate hits by e-value reliability
  # e-value < e-50: reliable for kingdom/phylum assignment
  reliable_hits = hit_evalues < evalue_reliable
  # e-value < e-20: can be used for lower ranks
  valid_hits_for_lower = hit_evalues < evalue_unknown
  
  # Use only valid hits
  hit_taxonomies_valid = hit_taxonomies[valid_hits_for_lower]
  hit_scores_valid = hit_scores[valid_hits_for_lower]
  hit_evalues_valid = hit_evalues[valid_hits_for_lower]
  reliable_mask = reliable_hits[valid_hits_for_lower]
  
  # Combine all taxonomies into one data frame
  all_tax = do.call(rbind, hit_taxonomies_valid)
  
  # Get consensus for each rank (most common, tie-break by best score)
  consensus = data.frame(
    Accession = NA, Kingdom = NA, Phylum = NA, Class = NA, 
    Order = NA, Family = NA, Genus = NA, Species = NA,
    stringsAsFactors = FALSE
  )
  
  ranks = c("Accession", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rank_thresholds = list(Class = 75, Order = 80, Family = 85, Genus = 90, Species = 97)
  
  # First, handle Kingdom and Phylum (using reliable hits only)
  for (rank in c("Kingdom", "Phylum")) {
    if (sum(reliable_mask, na.rm = TRUE) == 0) {
      consensus[[rank]] = NA
      if (rank == "Kingdom") {
        consensus$Kingdom = "unknown"
        return(consensus)
      }
      next
    }
    
    rank_values = all_tax[[rank]][reliable_mask]
    rank_scores = hit_scores_valid[reliable_mask]
    rank_values = rank_values[!is.na(rank_values) & rank_values != ""]
    
    if (length(rank_values) == 0) {
      consensus[[rank]] = NA
      if (rank == "Kingdom") {
        consensus$Kingdom = "unknown"
        return(consensus)
      }
      next
    }
    
    # Get consensus (most common, tie-break by best score)
    freq_table = table(rank_values)
    most_common = names(freq_table)[which.max(freq_table)]
    matching_indices = which(all_tax[[rank]] == most_common & reliable_mask)
    best_idx = matching_indices[which.max(hit_scores_valid[matching_indices])]
    
    consensus[[rank]] = most_common
    if (rank == "Kingdom" && hit_evalues_valid[best_idx] >= evalue_reliable) {
      consensus$Kingdom = "unknown"
      return(consensus)
    }
  }
  
  # Process ranks from Species down to Class
  # When there's disagreement, check if all hits agree at a higher rank
  taxonomic_ranks = c("Species", "Genus", "Family", "Order", "Class")
  
  # First, filter out hits that are >5 points below the best hit
  if (length(hit_scores_valid) > 0) {
    best_hit_score = max(hit_scores_valid, na.rm = TRUE)
    good_hits_mask = !is.na(hit_scores_valid) & hit_scores_valid >= (best_hit_score - 5)
    
    # Update valid hits to only include "good" hits
    hit_taxonomies_valid = hit_taxonomies_valid[good_hits_mask]
    hit_scores_valid = hit_scores_valid[good_hits_mask]
    hit_evalues_valid = hit_evalues_valid[good_hits_mask]
    reliable_mask = reliable_mask[good_hits_mask]
    
    # Rebuild all_tax with filtered hits
    all_tax = do.call(rbind, hit_taxonomies_valid)
  }
  
  # Process each rank and assign all that can be assigned
  for (rank in taxonomic_ranks) {
    threshold = rank_thresholds[[rank]]
    
    # Filter hits by sim_score threshold for this rank
    # Only use hits that meet the threshold for this rank
    score_threshold_mask = !is.na(hit_scores_valid) & hit_scores_valid >= threshold
    
    if (sum(score_threshold_mask, na.rm = TRUE) == 0) {
      # No hits meet the threshold for this rank
      consensus[[rank]] = NA
      # If this rank can't be assigned, all more specific ranks should also be NA
      current_idx = which(taxonomic_ranks == rank)
      if (current_idx > 1) {
        for (lower_idx in 1:(current_idx - 1)) {
          consensus[[taxonomic_ranks[lower_idx]]] = NA
        }
      }
      next
    }
    
    # Use only hits that meet the sim_score threshold
    rank_values = all_tax[[rank]][score_threshold_mask]
    rank_scores = hit_scores_valid[score_threshold_mask]
    
    # Filter out NA and empty values
    valid_idx = !is.na(rank_values) & rank_values != ""
    if (sum(valid_idx) == 0) {
      consensus[[rank]] = NA
      # If this rank can't be assigned, all more specific ranks should also be NA
      current_idx = which(taxonomic_ranks == rank)
      if (current_idx > 1) {
        for (lower_idx in 1:(current_idx - 1)) {
          consensus[[taxonomic_ranks[lower_idx]]] = NA
        }
      }
      next
    }
    
    rank_values_valid = rank_values[valid_idx]
    rank_scores_valid_subset = rank_scores[valid_idx]
    
    # Count frequencies and get best scores for each unique value
    unique_values = unique(rank_values_valid)
    value_stats = data.frame(
      value = unique_values,
      frequency = sapply(unique_values, function(v) sum(rank_values_valid == v)),
      best_score = sapply(unique_values, function(v) {
        idx = which(rank_values_valid == v)
        if (length(idx) > 0) max(rank_scores_valid_subset[idx], na.rm = TRUE) else NA
      }),
      stringsAsFactors = FALSE
    )
    
    # Sort by frequency (descending), then by best_score (descending)
    value_stats = value_stats[order(-value_stats$frequency, -value_stats$best_score), ]
    
    top_value = value_stats$value[1]
    top_freq = value_stats$frequency[1]
    top_score = value_stats$best_score[1]
    
    # Check if there's disagreement and minority has better score
    should_assign = TRUE
    if (nrow(value_stats) > 1) {
      second_best_score = value_stats$best_score[2]
      # If minority has a significantly better score, don't assign this rank
      if (!is.na(second_best_score) && !is.na(top_score) && 
          second_best_score > top_score + 2) {
        should_assign = FALSE
      }
    }
    
    if (should_assign) {
      # Find which hit has the top value with best score (from filtered hits)
      matching_indices_filtered = which(rank_values_valid == top_value)
      if (length(matching_indices_filtered) > 0) {
        # Get the best score from the filtered subset
        best_score_idx = matching_indices_filtered[which.max(rank_scores_valid_subset[matching_indices_filtered])]
        consensus_score = rank_scores_valid_subset[best_score_idx]
        
        # Threshold already applied in filtering, so consensus_score should meet it
        if (!is.na(consensus_score) && consensus_score >= threshold) {
          consensus[[rank]] = top_value
        } else {
          # Below threshold - don't assign this rank
          consensus[[rank]] = NA
          # If this rank can't be assigned, all more specific ranks should also be NA
          current_idx = which(taxonomic_ranks == rank)
          if (current_idx > 1) {
            for (lower_idx in 1:(current_idx - 1)) {
              consensus[[taxonomic_ranks[lower_idx]]] = NA
            }
          }
        }
      } else {
        consensus[[rank]] = NA
        # If this rank can't be assigned, all more specific ranks should also be NA
        current_idx = which(taxonomic_ranks == rank)
        if (current_idx > 1) {
          for (lower_idx in 1:(current_idx - 1)) {
            consensus[[taxonomic_ranks[lower_idx]]] = NA
          }
        }
      }
    } else {
      # Don't assign this rank due to disagreement
      consensus[[rank]] = NA
      # If this rank can't be assigned, all more specific ranks should also be NA
      current_idx = which(taxonomic_ranks == rank)
      if (current_idx > 1) {
        for (lower_idx in 1:(current_idx - 1)) {
          consensus[[taxonomic_ranks[lower_idx]]] = NA
        }
      }
    }
  }
  
  # If no rank was assigned at species/genus level, check if all hits agree at higher ranks
  if (is.na(consensus$Species) && is.na(consensus$Genus)) {
    for (rank in c("Family", "Order", "Class")) {
      if (!is.na(consensus[[rank]])) {
        # Already assigned, skip
        next
      }
      
      threshold = rank_thresholds[[rank]]
      
      # Filter hits by sim_score threshold for this rank
      score_threshold_mask = !is.na(hit_scores_valid) & hit_scores_valid >= threshold
      
      if (sum(score_threshold_mask, na.rm = TRUE) == 0) {
        next
      }
      
      # Use only hits that meet the sim_score threshold
      rank_values = all_tax[[rank]][score_threshold_mask]
      rank_values = rank_values[!is.na(rank_values) & rank_values != ""]
      
      if (length(unique(rank_values)) == 1 && length(rank_values) > 0) {
        # All hits (that meet threshold) agree at this rank - assign it
        matching_indices = which(all_tax[[rank]] == unique(rank_values)[1] & score_threshold_mask)
        if (length(matching_indices) > 0) {
          best_idx = matching_indices[which.max(hit_scores_valid[matching_indices])]
          consensus_score = hit_scores_valid[best_idx]
          
          if (!is.na(consensus_score) && consensus_score >= threshold) {
            consensus[[rank]] = unique(rank_values)[1]
          }
        }
      }
    }
  }
  
  return(consensus)
}

# Extract hits for each OTU from blast_10_hits
# Each row contains all 10 hits for one OTU
# Column pattern: qseqid, 1st_hit, qlen, slen, ..., qseqid, 2nd_hit, qlen, slen, ...

consensus_list = list()

# Process each row (each row is one OTU with up to 10 hits)
for (row_idx in 1:nrow(blast_10_hits)) {
  otu_row = blast_10_hits[row_idx, ]
  otu_id = otu_row$qseqid
  
  # Extract taxonomy strings, scores, and e-values for all hits (1st through 10th)
  hit_taxonomies = list()
  hit_scores = numeric()
  hit_evalues = numeric()
  
  # Column names: X1st_hit, X2nd_hit, ..., X10th_hit
  # Score columns: sim_score, sim_score.1, sim_score.2, ..., sim_score.9
  # E-value columns: evalue, evalue.1, evalue.2, ..., evalue.9
  hit_cols = c("X1st_hit", "X2nd_hit", "X3rd_hit", "X4th_hit", "X5th_hit",
               "X6th_hit", "X7th_hit", "X8th_hit", "X9th_hit", "X10th_hit")
  
  for (i in 1:10) {
    hit_col = hit_cols[i]
    
    # Check if this hit column exists and has data
    if (hit_col %in% colnames(otu_row)) {
      tax_string = otu_row[[hit_col]]
      
      # Skip if empty/NA
      if (!is.na(tax_string) && tax_string != "" && tax_string != "*") {
        # Parse taxonomy
        parsed = parse_taxonomy(tax_string)
        hit_taxonomies[[length(hit_taxonomies) + 1]] = parsed
        
        # Get sim_score for this hit
        if (i == 1) {
          score_col = "sim_score"
          evalue_col = "evalue"
        } else {
          score_col = paste0("sim_score.", i - 1)
          evalue_col = paste0("evalue.", i - 1)
        }
        
        # Try to get score, default to 100 if not found
        if (score_col %in% colnames(otu_row)) {
          hit_scores[length(hit_scores) + 1] = as.numeric(otu_row[[score_col]])
        } else {
          hit_scores[length(hit_scores) + 1] = 100
        }
        
        # Get e-value for this hit
        if (evalue_col %in% colnames(otu_row)) {
          hit_evalues[length(hit_evalues) + 1] = as.numeric(otu_row[[evalue_col]])
        } else {
          hit_evalues[length(hit_evalues) + 1] = 1e-100  # Very good default
        }
      }
    }
  }
  
  # Get consensus for this OTU
  consensus = get_consensus_taxonomy(hit_taxonomies, hit_scores, hit_evalues)
  consensus_list[[otu_id]] = consensus
}

# Convert to data frame
consensus_taxonomy = do.call(rbind, consensus_list)
rownames(consensus_taxonomy) = names(consensus_list)

# Combine with OTU identifiers and other info
blast_consensus = cbind(
  qseqid = names(consensus_list),
  consensus_taxonomy
)

# Remove Accession column (cannot be parsed to output)
if ("Accession" %in% colnames(blast_consensus)) {
  blast_consensus = blast_consensus[, !colnames(blast_consensus) %in% "Accession", drop = FALSE]
}

# Filter to target if specified
if (length(target) > 0 && !all(is.na(target)) && target[1] != "") {
  if (tax_level %in% colnames(blast_consensus)) {
    n_before = nrow(blast_consensus)
    blast_consensus = blast_consensus[blast_consensus[[tax_level]] %in% target, , drop = FALSE]
    n_after = nrow(blast_consensus)
    n_excluded = n_before - n_after
    cat("\n Consensus filtered to", tax_level, "level:", 
        paste(target, collapse = ", "), "\n")
    cat(" Excluded", n_excluded, "OTUs not matching target.\n")
    cat(" Remaining OTUs:", n_after, "\n")
  }
}

# Write consensus output
write.table(blast_consensus, file = "BLAST_consensus_parsed.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)



