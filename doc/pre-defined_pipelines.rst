.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text

.. meta::
    :description lang=en:
        PipeCraft2 manual. User guide for PipeCraft2

|PipeCraft2_logo|
  `github <https://github.com/pipecraft2/pipecraft>`_

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red

=====================
Pre-defined pipelines
=====================

.. _predefinedpipelines: 

vsearch OTUs
============

.. _otupipe:

OTUs workflow panel
-------------------

.. |otuoff| image:: _static/OTU_off.png
  :width: 50
  :alt: Alternative text

.. |otuon| image:: _static/OTU_on.png
  :width: 50
  :alt: Alternative text

.. note::
  This OTU workflow works with paired-end (e.g. Illumina, MGI-Tech) as well as single-end reads (e.g. PacBio, assembled Illumina reads)


This automated workflow is mostly based on `vsearch <https://github.com/torognes/vsearch>`_ (`Rognes et. al 2016 <https://peerj.com/articles/2584/>`_) [`manual <_static/vsearch_manual_2.22.1.pdf>`_]
 | Note that ``demultiplexing``, ``reorient`` and ``remove primers`` steps are optional. Nevertheless, it is advisable to :ref:`reorient <reorient>` your reads (to 5'-3') and :ref:`remove primers <remove_primers>` before proceeding.

 
.. _otupipe_defaults:

| **Default options:**
| *click on analyses step for more info*

==================================================================== =========================
Analyses step                                                        Default setting
==================================================================== =========================
:ref:`DEMULTIPLEX <demux>` (optional)                                 --
:ref:`REORIENT <reorient>` (optional)                                 --
:ref:`REMOVE PRIMERS <remove_primers>` (optional)                     --
:ref:`MERGE READS <merge_vsearch>`                                   | ``read_R1`` = \\.R1
                                                                     | ``min_overlap`` = 12
                                                                     | ``min_length`` = 32
                                                                     | ``allow_merge_stagger`` = TRUE 
                                                                     | ``include only R1`` = FALSE 
                                                                     | ``max_diffs`` = 20
                                                                     | ``max_Ns`` = 0
                                                                     | ``max_len`` = 600
                                                                     | ``keep_disjoined`` = FALSE 
                                                                     | ``fastq_qmax`` = 41
:ref:`QUALITY FILTERING with vsearch <qfilt_vsearch>`                | ``maxEE`` = 1
                                                                     | ``maxN`` = 0
                                                                     | ``minLen`` = 32
                                                                     | ``max_length`` = undefined
                                                                     | ``qmax`` = 41
                                                                     | ``qmin`` = 0
                                                                     | ``maxee_rate`` = undefined
:ref:`CHIMERA FILTERING with uchime_denovo <chimFilt_vsearch>`       | ``pre_cluster`` = 0.98
                                                                     | ``min_unique_size`` = 1
                                                                     | ``denovo`` = TRUE 
                                                                     | ``reference_based`` = undefined
                                                                     | ``abundance_skew`` = 2
                                                                     | ``min_h`` = 0.28
:ref:`ITS Extractor <itsextractor>` (optional)                       | ``organisms`` = all 
                                                                     | ``regions`` = all
                                                                     | ``partial`` = 50
                                                                     | ``region_for_clustering`` = ITS2
                                                                     | ``cluster_full_and_partial`` = TRUE
                                                                     | ``e_value`` = 1e-2
                                                                     | ``scores`` = 0
                                                                     | ``domains`` = 2
                                                                     | ``complement`` = TRUE 
                                                                     | ``only_full`` = FALSE
                                                                     | ``truncate`` = TRUE 
:ref:`CLUSTERING with vsearch <clustering_vsearch>`                  | ``OTU_type`` = centroid
                                                                     | ``similarity_threshold`` = 0.97
                                                                     | ``strands`` = both
                                                                     | ``remove_singletons`` = false
                                                                     | ``similarity_type`` = 2
                                                                     | ``sequence_sorting`` = cluster_size
                                                                     | ``centroid_type`` = similarity
                                                                     | ``max_hits`` = 1
                                                                     | ``mask`` = dust
                                                                     | ``dbmask`` = dust
:ref:`ASSIGN TAXONOMY with BLAST <assign_taxonomy_blast>` (optional) | ``database_file`` = select a database
                                                                     | ``task`` = blastn
                                                                     | ``strands`` = both
==================================================================== =========================

DADA2 ASVs
=============


.. _asvpipe:

ASVs workflow panel (with `DADA2 <https://benjjneb.github.io/dada2/index.html>`_)
---------------------------------------------------------------------------------

.. note::
  Working directory must contain **at least 2 samples** for DADA2 pipeline.


This automated workflow is based on the `DADA2 tutorial <https://benjjneb.github.io/dada2/tutorial.html>`_ 
 | Note that ``demultiplexing``, ``reorienting``, and ``primer removal`` steps are optional and do not represent parts from the DADA2 tutorial. Nevertheless, it is advisable to :ref:`remove primers <remove_primers>` before proceeding with ASV generation with DADA2.

| The official DADA2 manual is available `here <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_
 
.. _dada2_defaults:

**Default options:**

=========================================================== =========================
Analyses step                                               Default setting
=========================================================== =========================
:ref:`DEMULTIPLEX <demux>` (optional)                       | --
:ref:`REORIENT <reorient>` (optional)                       | --
:ref:`REMOVE PRIMERS <remove_primers>` (optional)           | --
:ref:`QUALITY FILTERING <dada2_qual_filt>`                  | ``read_R1`` = \\.R1
                                                            | ``read_R2`` = \\.R2
                                                            | ``maxEE`` = 2
                                                            | ``maxN`` = 0
                                                            | ``minLen`` = 20
                                                            | ``truncQ`` = 2
                                                            | ``truncLen`` = 0
                                                            | ``maxLen`` = 9999
                                                            | ``minQ`` = 2
                                                            | ``matchIDs`` = TRUE
:ref:`DENOISE <dada2_denoise>`                              | ``pool`` = FALSE
                                                            | ``selfConsist`` = FASLE
                                                            | ``qualityType`` = Auto
:ref:`MERGE PAIRED-END READS <dada2_merge_pairs>`           | ``minOverlap`` = 12
                                                            | ``maxMismatch`` = 0
                                                            | ``trimOverhang`` = FALSE
                                                            | ``justConcatenate`` = FALSE
:ref:`CHIMERA FILTERING <dada2_chimeras>`                   | ``method`` = consensus
:ref:`Filter ASV table <dada2_table_filtering>` (optional)  | ``collapseNoMismatch`` = TRUE
                                                            | ``by_length`` = 250
                                                            | ``minOverlap`` = 20
                                                            | ``vec`` = TRUE
:ref:`ASSIGN TAXONOMY <dada2_taxonomy>` (optional)          | ``minBoot`` = 50
                                                            | ``tryRC`` = FALSE
                                                            | ``dada2 database`` = select a database
=========================================================== =========================

____________________________________________________

.. _dada2_qual_filt:

QUALITY FILTERING [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `filterAndTrim <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function performs quality filtering on input FASTQ files based on user-selected criteria. Outputs include filtered FASTQ files located in the ``qualFiltered_out`` directory.

Quality profiles may be examined using the :ref:`QualityCheck module <interface>`.

================================ =========================
Setting                          Tooltip
================================ =========================
``read_R1``                      | applies only for **paired-end** data. 
                                 | Identifyer string that is common for all R1 reads 
                                 | (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. 
                                 | Note that backslash is only needed to escape dot regex; e.g. 
                                 | when all R1 files have '_R1' string, then enter '_R1'.). 
``read_R2``                      | applies only for **paired-end** data. 
                                 | Identifyer string that is common for all R2 reads 
                                 | (e.g. when all R2 files have '.R2' string, then enter '\\.R2'. 
                                 | Note that backslash is only needed to escape dot regex; e.g. 
                                 | when all R2 files have '_R1' string, then enter '_R2'.).
``maxEE``                        | discard sequences with more than the specified number of expected errors
``maxN``                         | discard sequences with more than the specified number of N’s (ambiguous bases)
``minLen``                       | remove reads with length less than minLen. minLen is enforced 
                                 | after all other trimming and truncation
``truncQ``                       | truncate reads at the first instance of a quality score less than or equal to truncQ
``truncLen``                     | truncate reads after truncLen bases 
                                 | (applies to **R1 reads** when working with **paired-end** data). 
                                 | Reads shorter than this are discarded. 
                                 | Explore quality profiles (with QualityCheck module) and 
                                 | see whether poor quality ends needs to be truncated
``truncLen_R2``                  | applies only for **paired-end** data. 
                                 | Truncate **R2 reads** after truncLen bases. 
                                 | Reads shorter than this are discarded. 
                                 | Explore quality profiles (with QualityCheck module) and 
                                 | see whether poor quality ends needs to truncated
``maxLen``                       | remove reads with length greater than maxLen. 
                                 | maxLen is enforced on the raw reads. 
                                 | In dada2, the default = Inf, but here set as 9999
``minQ``                         | after truncation, reads contain a quality score below minQ will be discarded
``matchIDs``                     | applies only for **paired-end** data. 
                                 | If TRUE, then double-checking (with seqkit pair) that only paired reads 
                                 | that share ids are outputted.
                                 | :red:`Note that 'seqkit' will be used for this process`, because when 
                                 | using e.g. SRA fastq files where original fastq headers have been 
                                 | replaced, dada2 does not recognize those fastq id strings
================================ =========================

see :ref:`default settings <dada2_defaults>`

____________________________________________________

.. _dada2_denoise:

DENOISING [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `dada <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to remove sequencing errors.
Outputs filtered fasta files into ``denoised_assembled.dada2`` directory.

==================== ============
Setting              Tooltip
==================== ============
``pool``             | if TRUE, the algorithm will pool together all samples prior to sample inference. 
                     | Pooling improves the detection of rare variants, but is computationally more expensive. 
                     | If pool = 'pseudo', the algorithm will perform pseudo-pooling between individually 
                     | processed samples.
``selfConsist``      | if TRUE, the algorithm will alternate between sample inference and error rate estimation 
                     | until convergence
``qualityType``      | 'Auto' means to attempt to auto-detect the fastq quality encoding. 
                     | This may fail for PacBio files with uniformly high quality scores, 
                     | in which case use 'FastqQuality'
==================== ============

see :ref:`default settings <dada2_defaults>`

____________________________________________________

.. _dada2_merge_pairs:

MERGE PAIRS [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `mergePairs <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to merge paired-end reads. 
Outputs merged fasta files into ``denoised_assembled.dada2`` directory.

==================== ============
Setting               Tooltip
==================== ============
``minOverlap``       | the minimum length of the overlap required for merging the forward and reverse reads
``maxMismatch``      | the maximum mismatches allowed in the overlap region
``trimOverhang``     | if TRUE, overhangs in the alignment between the forwards and reverse read are  
                     | trimmed off. Overhangs are when the reverse read extends past the start of 
                     | the forward read, and vice-versa, as can happen when reads are longer than the 
                     | amplicon and read into the other-direction primer region
``justConcatenate``  | if TRUE, the forward and reverse-complemented reverse read are concatenated  
                     | rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them
==================== ============

see :ref:`default settings <dada2_defaults>`

.. _dada2_chimeras:

____________________________________________________

CHIMERA FILTERING [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `removeBimeraDenovo <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to remove chimeras. 
Outputs filtered fasta files into ``chimeraFiltered_out.dada2`` and final ASVs to ``ASVs_out.dada2`` directory.

==================== ============
Setting               Tooltip
==================== ============
``method``           | 'consensus' - the samples are independently checked for chimeras, and a consensus 
                     | decision on each sequence variant is made. 
                     | If 'pooled', the samples are all pooled together for chimera identification. 
                     | If 'per-sample', the samples are independently checked for chimeras
==================== ============

see :ref:`default settings <dada2_defaults>`

.. _dada2_table_filtering:

____________________________________________________

filter ASV table [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `collapseNoMismatch <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to collapse identical ASVs; 
and ASVs filtering based on minimum accepted sequence length (custom R functions). 
Outputs filtered ASV table and fasta files into ``ASVs_out.dada2/filtered`` directory.

========================== ============
Setting                    Tooltip
========================== ============
``collapseNoMismatch``     | collapses ASVs that are identical up to shifts or 
                           | length variation, i.e. that have no mismatches or internal indels
``by_length``              | discard ASVs from the ASV table that are shorter than specified 
                           | value (in base pairs). Value 0 means OFF, no filtering by length
``minOverlap``             | collapseNoMismatch setting. Default = 20. The minimum overlap of 
                           | base pairs between ASV sequences required to collapse them together
``vec``                    | collapseNoMismatch setting. Default = TRUE. Use the vectorized 
                           | aligner. Should be turned off if sequences exceed 2kb in length
========================== ============

see :ref:`default settings <dada2_defaults>`

____________________________________________________

.. _dada2_taxonomy:

ASSIGN TAXONOMY [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `assignTaxonomy <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to classify ASVs. 
Outputs classified fasta files into ``taxonomy_out.dada2`` directory.

==================== ============
Setting               Tooltip
==================== ============
``minBoot``          | the minimum bootstrap confidence for assigning a taxonomic level
``tryRC``            | the reverse-complement of each sequences will be used for classification 
                     | if it is a better match to the reference sequences than the forward sequence
``dada2 database``   | select a reference database fasta file for taxonomy annotation
                     | `Download DADA2-formatted reference databases here <https://benjjneb.github.io/dada2/training.html>`_
==================== ============

see :ref:`default settings <dada2_defaults>`

____________________________________________________


UNOISE ASVs
===========

UNOISE3 pipeline for making ASVs (zOTUs) + optionally automatic clustering of those ASVs.
Updating this section soon.


.. _nextits: 

NextITS
=======

NextITS pipeline for PacBio ITS sequences. 
Updating this section soon.
