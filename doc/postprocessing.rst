.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 50
  :alt: Alternative text
  :target: https://github.com/pipecraft2/user_guide

.. |output_icon| image:: _static/output_icon.png
  :width: 50
  :alt: Alternative text

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red

.. raw:: html

    <style> .green {color:#00f03d; font-weight:bold; font-size:16px} </style>

.. role:: green
  

.. meta::
    :description lang=en:
        PipeCraft manual. tutorial

.. _postprocessingtools:

=======================================
Post-processing tools |PipeCraft2_logo|
=======================================

.. note:: 

  All post-processing tools accessible under **ADD STEP** -> **POSTPROCESSING**

.. contents:: Contents
   :depth: 2

.. important::

  Note that the input ASV/OTU table must contain "Sequence" column for some postprocessing modules. Read the tooltips too see which ones. 


____________________________________________________

Filter tag-jumps
----------------

Filter out putative tag-jumps with `UNCROSS2 <https://www.drive5.com/usearch/manual/uncross2_algo.html>`_. 

| Input data is tab delimited **OTU table** and corresponding **fasta file** (representative sequences of ASV/OTUs).

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR``, but the file formats do not matter here (just click 'Next').

+--------------------+-------------------------------------------------------------------------------+
| Outputs            |                                                                               |
+====================+===============================================================================+
| *_TagJumpFilt.txt  | output table where tag-jumps have been filtered out                           |
+--------------------+-------------------------------------------------------------------------------+
| TagJump_plot.pdf   | illustration about the presence of tag-jumps based on the selected parameters |
+--------------------+-------------------------------------------------------------------------------+
|| TagJump_stats.txt || tag-jump statistics (Total_reads, Number_of_TagJump_Events,                  |
||                   || TagJump_reads, ReadPercent_removed)                                          |
+--------------------+-------------------------------------------------------------------------------+

+-----------------+------------------------------------------------------------------------------------------+
| Setting         | Tooltip                                                                                  |
+=================+==========================================================================================+
|| ``table``      || select tab-delimited OTU/ASV table, where the 1st column is the OTU/ASV IDs and the     |
||                || following columns represent samples; 2nd column may be Sequence column, with the        |
||                || colName 'Sequence' [**file must be in the SELECT WORKDIR directory**]                   |
+-----------------+------------------------------------------------------------------------------------------+
|| ``fasta file`` || select corresponding fasta file for OTU/ASV table [**fasta file must be in the SELECT** |
||                || **WORKDIR directory**]                                                                  |
+-----------------+------------------------------------------------------------------------------------------+
|| ``f value``    || f-parameter of UNCROSS2, which defines the expected tag-jumps rate. Default is 0.03     |
||                || (equivalent to 3%). A higher value enforces stricter filtering                          |
+-----------------+------------------------------------------------------------------------------------------+
|| ``p value``    || p-parameter, which controls the severity of tag-jump removal. It adjusts the exponent   |
||                || in the UNCROSS formula. Default is 1. Opt for 0.5 or 0.3 to steepen the curve           |
+-----------------+------------------------------------------------------------------------------------------+

___________________________________________________

ASV to OTU
----------

Cluster ASVs (zOTUs) to OTUs using vsearch. 

| **Input data** is tab delimited **ASV table** and **ASV sequences** in fasta format.
| 2nd column of **ASV table** MUST BE 'Sequences' (1st column is ASV IDs; default pipecraft output table).
| For clustering, the ASV size annotation is obtained from the ASV table. 


+---------------------------+--------------------------------------------------------+
| Outputs ``ASVs2OTUs_out`` |                                                        |
+===========================+========================================================+
| OTUs.fasta                | FASTA formated representative OTU sequences            |
+---------------------------+--------------------------------------------------------+
| OTU_table.txt             | OTU distribution table per sample (tab delimited file) |
+---------------------------+--------------------------------------------------------+
| OTUs.uc                   | uclust-like formatted clustering results for OTUs      |
+---------------------------+--------------------------------------------------------+
| ASVs.size.fasta           | size annotated input sequences                         |
+---------------------------+--------------------------------------------------------+


.. _postclustering_lulu:

___________________________________________________

LULU post-clustering
---------------------

Perform OTU post-clustering with `LULU <https://github.com/tobiasgf/lulu>`_ to merge co-occurring 'daughter' OTUs.

LULU description from the `LULU repository <https://github.com/tobiasgf/lulu>`_: the purpose of LULU is to reduce the number of 
erroneous OTUs in OTU tables to achieve more realistic biodiversity metrics. 
By evaluating the co-occurence patterns of OTUs among samples LULU identifies OTUs that consistently satisfy some user selected 
criteria for being errors of more abundant OTUs and merges these. It has been shown that curation with LULU consistently result 
in more realistic diversity metrics. 

Additional information:
 - `LULU repository <https://github.com/tobiasgf/lulu>`_
 - `LULU paper <https://doi.org/10.1038/s41467-017-01312-x>`_
  
| Input data is tab delimited **OTU table** (``table``) and **OTU sequences** (``rep_seqs``) in fasta format (see input examples below). 
| `EXAMPLE table here <https://github.com/tobiasgf/lulu/blob/master/Example_data/otutable_test.txt>`_ *(from LULU repository)*
| `EXAMPLE fasta here <https://github.com/tobiasgf/lulu/blob/master/Example_data/centroids_test.txt>`_ *(from LULU repository)*

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR``, but the file formats do not matter here (just click 'Next').


+------------------------+----------------------------------------------------------------------------+
| Outputs ``lulu_out``   |                                                                            |
+========================+============================================================================+
| lulu_out_table.txt     | curated table in tab delimited txt format                                  |
+------------------------+----------------------------------------------------------------------------+
| lulu_out_RepSeqs.fasta | fasta file for the molecular units (OTUs or ASVs) in the curated table     |
+------------------------+----------------------------------------------------------------------------+
| match_list.lulu        | match list file that was used by LULU to merge 'daughter' molecular units  |
+------------------------+----------------------------------------------------------------------------+
|| discarded_units.lulu  || molecular units (OTUs or ASVs) that were merged with other units based on |
||                       || specified thresholds                                                      |
+------------------------+----------------------------------------------------------------------------+

=============================================== =========================
`Setting <https://github.com/tobiasgf/lulu>`_   Tooltip
=============================================== =========================
``table``                                       | select OTU/ASV table. If no file is selected, then PipeCraft will 
                                                | look OTU_table.txt or ASV_table.txt in the working directory.
                                                | `EXAMPLE table here <https://github.com/tobiasgf/lulu/blob/master/Example_data/otutable_test.txt>`_
``fasta_file``                                  | select fasta formatted sequence file containing your OTU/ASV reads.
                                                | `EXAMPLE file here <https://github.com/tobiasgf/lulu/blob/master/Example_data/centroids_test.txt>`_
``min_ratio_type``                              | sets whether a potential error must have lower abundance than the parent 
                                                | in all samples 'min' (default), or if an error just needs to have lower 
                                                | abundance on average 'avg'
``min_ratio``                                   | set the minimim abundance ratio between a potential error and a 
                                                | potential parent to be identified as an error
``min_match``                                   | specify minimum threshold of sequence similarity for considering 
                                                | any OTU as an error of another
``min_rel_cooccurence``                         | minimum co-occurrence rate. Default = 0.95 (meaning that 1 in 20 samples 
                                                | are allowed to have no parent presence)
``match_list_soft``                             | use either 'blastn' or 'vsearch' to generate match list for LULU. 
                                                | Default is 'vsearch' (much faster)
``vsearch_similarity_type``                     | applies only when 'vsearch' is used as 'match_list_soft'. 
                                                | Pairwise sequence identity definition (--iddef)
``perc_identity``                               | percent identity cutoff for match list. Excluding pairwise comparisons 
                                                | with lower sequence identity percentage than specified threshold
``coverage_perc``                               | percent query coverage per hit. Excluding pairwise comparisons with 
                                                | lower sequence coverage than specified threshold
``strands``                                     | query strand to search against database. Both = search also reverse complement
``cores``                                       | number of cores to use for generating match list for LULU
=============================================== ========================= 


.. _postclustering_dada2_table_filtering:

____________________________________________________

DADA2 collapse ASVs
-------------------

DADA2 `collapseNoMismatch <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function collapses identical ASVs with no internal mismatches (~greedy 100% clustering with end-gapping ignored).
Representative sequence of a collapsed ASV will be the most abundant one. 
and ASVs filtering based on minimum accepted sequence length (custom R functions). 

To **START**, specify working directory under ``SELECT WORKDIR``, but the file formats do not matter here (just click 'Next').


+---------------------------------+-------------------------------------------------------------------------------------+
| Outputs ``filtered_table``      |                                                                                     |
+=================================+=====================================================================================+
| ASVs_table_collapsed.txt        | ASV table after collapsing identical ASVs                                           |
+---------------------------------+-------------------------------------------------------------------------------------+
| ASVs_collapsed.fasta            | ASV sequences after collapsing identical ASVs                                       |
+---------------------------------+-------------------------------------------------------------------------------------+
| ASV_table_collapsed.rds         | ASV table in RDS format after collapsing identical ASVs                             |
+---------------------------------+-------------------------------------------------------------------------------------+
| If length filtering was applied |                                                                                     |
+---------------------------------+-------------------------------------------------------------------------------------+
| ASV_table_lenFilt.tx            | ASV table after filtering out ASVs with shorther than specified sequence length     |
+---------------------------------+-------------------------------------------------------------------------------------+
| ASVs_lenFilt.fasta              | ASV sequences after filtering out ASVs with shorther than specified sequence length |
+---------------------------------+-------------------------------------------------------------------------------------+

========================== ============
Setting                    Tooltip
========================== ============
``DADA2 table``            | select the RDS file (ASV table), output from DADA2 workflow; 
                           | usually in ASVs_out.dada2/ASVs_table.denoised-merged.rds
``collapseNoMismatch``     | collapses ASVs that are identical up to shifts or 
                           | length variation, i.e. that have no mismatches or internal indels
``by_length``              | discard ASVs from the ASV table that are shorter than specified 
                           | value (in base pairs). Value 0 means OFF, no filtering by length
``minOverlap``             | collapseNoMismatch setting. Default = 20. The minimum overlap of 
                           | base pairs between ASV sequences required to collapse them together
``vec``                    | collapseNoMismatch setting. Default = TRUE. Use the vectorized 
                           | aligner. Should be turned off if sequences exceed 2kb in length
========================== ============

__________________________________________________

metaMATE
--------

Determine and filter out putative NUMTs (from mitochondrial coding amplicon genes) and and other erroneous sequences based on relative read abundance thresholds within libraries, phylogenetic clades and/or taxonomic groupings.
 
Additional information:
 - `metaMATE repository <https://github.com/tjcreedy/metamate>`_
 - `metaMATE paper <https://doi.org/10.1111/1755-0998.13337>`_
  

___________________________________________________

ORF-Finder
----------

Filter out putative pseudogenes (NUMTs) from protein coding amplicon dataset (such as COI, rbcL) using NCBI's ORFfinder `(Sayers et al 2022) <https://doi.org/10.1093/nar/gkab1112>`_.
This process translates sequences to open reading frames (ORFs) and retaines the longest ORF per sequence 
if the length of the ORF is between the specified range of ``min length`` and ``max length``.


.. _assign_taxonomy:

____________________________________________________


.. _postprocessing_deicode: 

`DEICODE <https://github.com/biocore/DEICODE>`_ 
-----------------------------------------------

DEICODE (`Martino et al., 2019 <https://doi.org/10.1128/mSystems.00016-19>`_) is used to perform beta diversity analysis 
by applying robust Aitchison PCA on the OTU/ASV table. To consider the compositional nature of data, 
it preprocesses data with rCLR transformation (centered log-ratio on only non-zero values, without adding pseudo count). 
As a second step, it performs dimensionality reduction of the data using robust PCA (also applied only to the non-zero values of the data), 
where sparse data are handled through matrix completion.

Additional information:
 - `DEICODE repository <https://github.com/biocore/DEICODE>`_
 - `DEICODE paper <https://journals.asm.org/doi/10.1128/mSystems.00016-19>`_



| Input data is tab delimited **OTU table** and optionally **subset of OTU ids** to generate results also for the selected subset (see input examples below). 

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR``, but the file formats do not matter here (just click 'Next').

+-------------------------------------------------------------------+------------------------------------------------------------------------+
| Output directory |output_icon| ``DEICODE_out``                                                                                             |
+===================================================================+========================================================================+
| otutab.biom                                                       | full OTU table in BIOM format                                          |
+-------------------------------------------------------------------+------------------------------------------------------------------------+
| rclr_subset.tsv                                                   | rCLR-transformed subset of OTU table \*                                |
+-------------------------------------------------------------------+------------------------------------------------------------------------+
| ``full``/distance-matrix.tsv                                      | distance matrix between the samples, based on full OTU table           |
+-------------------------------------------------------------------+------------------------------------------------------------------------+
| ``full``/ordination.txt                                           | ordination scores for samples and OTUs, based on full OTU table        |
+-------------------------------------------------------------------+------------------------------------------------------------------------+
| ``full``/rclr.tsv                                                 | rCLR-transformed OTU table                                             |
+-------------------------------------------------------------------+------------------------------------------------------------------------+
| ``subs``/distance-matrix.tsv                                      | distance matrix between the samples, based on a subset of OTU table \* |
+-------------------------------------------------------------------+------------------------------------------------------------------------+
| ``subs``/ordination.txt                                           | ordination scores for samples and OTUs, based a subset of OTU table \* |
+-------------------------------------------------------------------+------------------------------------------------------------------------+
| \* files are present only if 'subset_IDs' variable was specified  |                                                                        |
+-------------------------------------------------------------------+------------------------------------------------------------------------+

=============================================== =========================
Setting                                         Tooltip
=============================================== =========================
``table``                                       | select OTU/ASV table. If no file is selected, then PipeCraft will 
                                                | look OTU_table.txt or ASV_table.txt in the working directory.
                                                | See OTU table example below
``subset_IDs``                                  | select list of OTU/ASV IDs for analysing a subset from the full table
                                                | see subset_IDs file example below
``min_otu_reads``                               | cutoff for reads per OTU/ASV. OTUs/ASVs with lower reads then specified 
                                                | cutoff will be excluded from the analysis
``min_sample_reads``                            | cutoff for reads per sample. Samples with lower reads then 
                                                | specified cutoff will be excluded from the analysis
=============================================== =========================


Example of input ``table`` (tab delimited text file):

================== ============== ============== ============== ==============
OTU_id             sample1        sample2        sample3        sample4
================== ============== ============== ============== ==============
00fc1569196587dde  106            271            584            20
02d84ed0175c2c79e  81             44             88             14
0407ee3bd15ca7206  3              4              3              0
042e5f0b5e38dff09  20             83             131            4
07411b848fcda497f  1              0              2              0
07e7806a732c67ef0  18             22             83             7
0836d270877aed22c  1              1              0              0
0aa6e7da5819c1197  1              4              5              0
0c1c219a4756bb729  18             17             40             7
================== ============== ============== ============== ==============

Example of input ``subset_IDs``:

.. code-block::

  07411b848fcda497f
  042e5f0b5e38dff09
  0836d270877aed22c
  0c1c219a4756bb729
  ...

| 


**PERMANOVA and PERMDISP example using the robust Aitchison distance**

.. code-block::

      library(vegan)

      ## Load distance matrix
      dd <- read.table(file = "distance-matrix.tsv")

      ## You will also need to load the sample metadata
      ## However, for this example we will create a dummy data
      meta <- data.frame(
        SampleID = rownames(dd),
        TestData = rep(c("A", "B", "C"), each = ceiling(nrow(dd)/3))[1:nrow(dd)])

      ## NB! Ensure that samples in distance matrix and metadata are in the same order
      meta <- meta[ match(x = meta$SampleID, table = rownames(dd)), ]

      ## Convert distance matrix into 'dist' class
      dd <- as.dist(dd)

      ## Run PERMANOVA
      adon <- adonis2(formula = dd ~ TestData, data = meta, permutations = 1000)
      adon

      ## Run PERMDISP
      permdisp <- betadisper(dd, meta$TestData)
      plot(permdisp)

Example of plotting the ordination scores

.. code-block::

      library(ggplot2)

      ## Load ordination scores
      ord <- readLines("ordination.txt")

      ## Skip PCA summary
      ord <- ord[ 8:length(ord) ]

      ## Break the data into sample and species scores
      breaks <- which(! nzchar(ord))
      ord <- ord[1:(breaks[2]-1)]               # Skip biplot scores
      ord_sp <- ord[1:(breaks[1]-1)]            # species scores
      ord_sm <- ord[(breaks[1]+2):length(ord)]  # sample scores

      ## Convert scores to data.frames 
      ord_sp <- as.data.frame( do.call(rbind, strsplit(x = ord_sp, split = "\t")) )
      colnames(ord_sp) <- c("OTU_ID", paste0("PC", 1:(ncol(ord_sp)-1)))

      ord_sm <- as.data.frame( do.call(rbind, strsplit(x = ord_sm, split = "\t")) )
      colnames(ord_sm) <- c("Sample_ID", paste0("PC", 1:(ncol(ord_sm)-1)))

      ## Convert PCA to numbers
      ord_sp[colnames(ord_sp)[-1]] <- sapply(ord_sp[colnames(ord_sp)[-1]], as.numeric)
      ord_sm[colnames(ord_sm)[-1]] <- sapply(ord_sm[colnames(ord_sm)[-1]], as.numeric)

      ## At this step, sample and OTU metadata could be added to the data.frame

      ## Example plot
      ggplot(data = ord_sm, aes(x = PC1, y = PC2)) + geom_point()

___________________________________________________

BlasCh
----------

False positive chimera detection and recovery module. BlasCh processes BLAST XML results to identify, classify, and recover sequences that were incorrectly flagged as chimeric during initial chimera detection.

| Input data is **chimera sequences** in fasta format (`.chimeras.fasta` files) and optionally a **reference database** (FASTA file or existing BLAST database).


.. note::

  To **START**, specify working directory under ``SELECT WORKDIR``, but the file formats do not matter here (just click 'Next').

+--------------------------------------+-------------------------------------------------------------------------+
| Output directory |output_icon| ``BlasCh_out``                                                                       |
+======================================+=========================================================================+
| *_non_chimeric.fasta                 | recovered non-chimeric sequences (rescued)                             |
+--------------------------------------+-------------------------------------------------------------------------+
| *_borderline.fasta                   | borderline sequences (also rescued as non-chimeric)                    |
+--------------------------------------+-------------------------------------------------------------------------+
| chimera_recovery_report.txt          | summary statistics and classification results                           |
+--------------------------------------+-------------------------------------------------------------------------+
| ``analysis``/*_chimeric.fasta        | confirmed chimeric sequences                                            |
+--------------------------------------+-------------------------------------------------------------------------+
| ``analysis``/*_multiple_alignments.fasta | sequences with multiple HSPs and low coverage                     |
+--------------------------------------+-------------------------------------------------------------------------+
| ``analysis``/*_sequence_details.csv  | detailed classification results for each sequence                      |
+--------------------------------------+-------------------------------------------------------------------------+
| ``databases``/                       | self-databases created from sample FASTA files                         |
+--------------------------------------+-------------------------------------------------------------------------+
| ``reference_db``/                    | reference database files (if created from FASTA input)                 |
+--------------------------------------+-------------------------------------------------------------------------+
| ``xml``/*_blast_results.xml          | BLAST XML output files for each sample                                 |
+--------------------------------------+-------------------------------------------------------------------------+

=============================================== =========================
Setting                                         Tooltip
=============================================== =========================
``reference_db``                                | path to reference database (FASTA file or existing BLAST database). 
                                                | Optional - if not provided, uses only self-databases
``threads``                                     | number of CPU threads for BLAST analysis (default: 8)
``high_identity_threshold``                     | identity threshold for high-quality matches (default: 99.0%)
``high_coverage_threshold``                     | coverage threshold for high-quality matches (default: 99.0%)
``borderline_identity_threshold``               | identity threshold for borderline recovery (default: 80.0%)
``borderline_coverage_threshold``               | coverage threshold for borderline recovery (default: 89.0%)
=============================================== =========================

**Classification Logic:**

The module uses a multi-tier classification system:

1. **Non-chimeric**: High identity (≥99%) and coverage (≥99%) matches to reference database
2. **Borderline**: Moderate identity (≥80%) and coverage (≥89%) - recovered as non-chimeric
3. **Chimeric**: Multiple taxonomies or only self-hits without reference matches
4. **Multiple alignments**: Sequences with multiple HSPs and coverage ≤85%

**Recovered Sequences:**

Both non-chimeric and borderline sequences are saved as "rescued" sequences in the main output directory, while confirmed chimeric and multiple alignment sequences are stored in the ``analysis`` subdirectory.

.. note::

  BlasCh automatically detects .chimeras.fasta files in the working directory and creates self-databases from available sample FASTA files. 
  Original sample files are prioritized over .chimeras files for database creation.