.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text

.. |otu_main| image:: _static/otu_main.png
  :width: 1500
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

|PipeCraft2_logo|
  `github <https://github.com/pipecraft2/pipecraft>`_


.. _postprocessingtools:

=====================
Post-processing tools
=====================

.. note:: 

  All post-processing tools accessible under **ADD STEP** -> **POSTPROCESSING**

.. contents:: Contents
   :depth: 2


____________________________________________________

.. _postprocessing_lulu:

`LULU <https://github.com/tobiasgf/lulu>`_ 
-------------------------------------------

LULU description from the `LULU repository <https://github.com/tobiasgf/lulu>`_: the purpose of LULU is to reduce the number of 
erroneous OTUs in OTU tables to achieve more realistic biodiversity metrics. 
By evaluating the co-occurence patterns of OTUs among samples LULU identifies OTUs that consistently satisfy some user selected 
criteria for being errors of more abundant OTUs and merges these. It has been shown that curation with LULU consistently result 
in more realistic diversity metrics. 

| This is implemented also under POSTCLUSTERING panel, :ref:`see here <postclustering_lulu>` 

____________________________________________________

.. _postprocessing_deicode: 

`DEICODE <https://github.com/biocore/DEICODE>`_ 
-------------------------------------------------

DEICODE (`Martino et al., 2019 <https://doi.org/10.1128/mSystems.00016-19>`_) is used to perform beta diversity analysis 
by applying robust Aitchison PCA on the OTU/ASV table. To consider the compositional nature of data, 
it preprocesses data with rCLR transformation (centered log-ratio on only non-zero values, without adding pseudo count). 
As a second step, it performs dimensionality reduction of the data using robust PCA (also applied only to the non-zero values of the data), 
where sparse data are handled through matrix completion.

Additional information:
 - `DEICODE tutorial <https://library.qiime2.org/plugins/deicode/19/>`_
 - `DEICODE repository <https://github.com/biocore/DEICODE>`_
 - `DEICODE paper <https://journals.asm.org/doi/10.1128/mSystems.00016-19>`_



| Input data is tab delimited **OTU table** and optionally **subset of OTU ids** to generate results also for the selected subset (see input examples below). 

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR``, but the file formats do not matter here (just click 'Next').

| **Output** files in ``DEICODE_out`` directory:
| #   - otutab.biom          =  full OTU table in BIOM format
| #   - rclr_subset.tsv      =  rCLR-transformed subset of OTU table *
| # DEICODE_out/full/
| #   - distance-matrix.tsv  =  distance matrix between the samples, based on full OTU table
| #   - ordination.txt       =  ordination scores for samples and OTUs, based on full OTU table
| #   - rclr.tsv             =  rCLR-transformed OTU table
| # DEICODE_out/subs/
| #   - distance-matrix.tsv  =  distance matrix between the samples, based on a subset of OTU table *
| #   - ordination.txt       =  ordination scores for samples and OTUs, based a subset of OTU table *
| # \*, files are present only if 'subset_IDs' variable was specified


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


