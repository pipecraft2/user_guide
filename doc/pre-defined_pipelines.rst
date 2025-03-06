.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 50
  :alt: Alternative text
  :target: https://github.com/pipecraft2/user_guide
  
.. |NextITS_seq_cluster| image:: _static/nextits_sequence_clustering.png
  :width: 600
  :height: 200
  :alt: Alternative text

.. |NextITS_extraction| image:: _static/nextits_extraction.png
  :width: 300
  :height: 200
  :alt: Alternative text

.. meta::
    :description lang=en:
        PipeCraft2 manual. User guide for PipeCraft2

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red

.. _predefinedpipelines: 

=======================================
Pre-defined pipelines |PipeCraft2_logo|
=======================================

Working with multiple sequencing runs
======================================

Applicable to: :ref:`DADA2 ASVs <asvpipe>`, :ref:`UNOISE ASVs <unoise_asvs>`,
:ref:`vsearch OTUs <vsearchOTUs>` pre-defined pipelines.


When working with multiple sequencing runs, then **pre-defined pipelines** can automatically process each sequencing run separately, and 
then **merge the results** into a single output OTU/ASV table. Processing each sequencing run separately is 
**necessary for appropriate handling of run-specifiec error profiles and tag-jumps filtering**.

Identical sequences from different runs will be recognized as the same ASV, and therefore merged into a single ASV. 

Directory structure
-------------------

.. important:: 

  When aiming to combine samples from multiple sequencing runs, then follow the below directory structure, 
  with :red:`**multiRunDir**` **being the mandatory directory name** (names of the nested sequencing run directories can be changed).
  
  When specifying a **working directory** in PipeCraft2 for processing multiple sequencing runs, 
  then select the **parent directory** of the **multiRunDir** (e.g. **my_sequencing_runs** in the example below).

| ├─── my_sequencing_runs     *# SELECT THIS FOLDER AS WORKING DIRECTORY*
| │   └─── **multiRunDir**
| │       ├─── **Run1**
| │       │   ├───*sample1_R1.fastq*
| │       │   ├───*sample1_R2.fastq*
| │       │   ├───*sample2_R1.fastq*
| │       │   ├───*sample2_R2.fastq*
| │       │   ├───...
| │       ├─── **Run2**
| │       │   ├───*sample10_R1.fastq*
| │       │   ├───*sample10_R2.fastq*
| │       │   ├───*sample11_R1.fastq*
| │       │   ├───*sample11_R2.fastq*
| │       │   ├───...
| │       ├─── **skip_Run3**   *# this dir will be skipped*
| │       │   ├───*sample20_R1.fastq*
| │       │   ├───*sample20_R2.fastq*
| │       │   ├───*sample21_R1.fastq*
| │       │   ├───*sample21_R2.fastq*
| │       │   ├───...
| │       └─── **merged_runs** *# this is the dir where the merged ASV/OTU table will be saved*
| │           ├───*ASVs.fasta*
| │           ├───*ASV_table.txt*
| │           ├───...

| Note that :red:`you can **skip** processing any sequencing run` by adding a **skip_** prefix to the directory name. In this example here, sequencing run ``skip_Run3`` will be skipped.
|
| ``merged_runs`` directory will contain the merged ASV/OTU table; :red:`avoid naming your sequencing run directories as **merged_runs**!`  
|
| Fastq files with the **same name** will be considered as the same sample and will be merged in the final ASV/OTU table.


Merge sequencing runs
---------------------

When working with multiple sequencing runs, then you can merge the results into a single ASV/OTU table
by enabling the **MERGE SEQUENCING RUNS** option in the 
:ref:`DADA2 ASVs <asvpipe>`, :ref:`UNOISE ASVs <unoise_asvs>`,
:ref:`vsearch OTUs <vsearchOTUs>` pre-defined pipelines.

Note that NextITS and OptimOTU pipelines also support merging of sequencing runs, but require 
slightly different directory structure (see here for NextITS: :ref:`nextits_pipeline` and 
for OptimOTU: :ref:`optimotu_pipeline`).


___________________________________________________

.. _asvpipe:

DADA2 ASVs
=============

This pre-defined workflow is based on the `DADA2 tutorial <https://benjjneb.github.io/dada2/tutorial.html>`_ to form **ASVs and an ASV table**.
This input is the directory that contains per-sample fastq files (**demultiplexed data**).

| Note that ``primer removal`` step do not represent parts from the DADA2 tutorial. Nevertheless, it is advisable to :ref:`remove primers <remove_primers>` before proceeding with ASV generation with DADA2.


.. _dada2_modes:

============= =====================
pipeline mode  when do use
============= =====================
**FORWARD**   | for paired-end Illumina data where amplicons 
              | are expected to be in uniform orientation 
              | (e.g., all reads are in 5'-3' orientation).
**MIXED**     | for paired-end Illumina data where amplicons 
              | are expected to be both, in 5'-3' (forward) 
              | and 3'-5' (reverse) oriented. 
              | In that mode, ``CUT PRIMERS`` is mandatory, 
              | and generates separate directories for forward 
              | and reverse oriented sequences, which will pass 
              | DADA2 pipeline individually. After merging the paired ends, 
              | the reverse oriented sequences are reverse complemented 
              | and aggregated with the forward reads for chimera filtering 
              | and ASV table generation. The output ASVs are all 5'-3' oriented. 
              | **Here, make sure that the R1 and R2 identifiers are** 
              | **\\.R1 and \\.R2 in the QUALITY FILTERING step**
**PACBIO**    | for single-end PacBio data. ``CUT PRIMERS`` step for single-end 
              | data will reoriente all reads to 5'-3' (forward) orientation.
============= =====================


.. note::
  Working directory must contain **at least 2 samples** for DADA2 pipeline.

 
.. _dada2_defaults:

**Default options:**

=========================================================== =========================
Analyses step                                               Default setting
=========================================================== =========================
:ref:`CUT PRIMERS <remove_primers>` (optional)              | --
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

QUALITY FILTERING
-----------------

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
``maxN``                         | discard sequences with more than the specified number of N's (ambiguous bases)
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

DENOISING
---------

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

MERGE PAIRS
-----------

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

CHIMERA FILTERING
-----------------

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

filter ASV table
----------------

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

ASSIGN TAXONOMY
---------------

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


.. _unoise_asvs:

UNOISE ASVs
===========

UNOISE3 pipeline for making ASVs (zOTUs). Uses UNOISE3 algorithm in vsearch. 

This automated workflow is mostly based on `vsearch <https://github.com/torognes/vsearch>`_ (`Rognes et. al 2016 <https://peerj.com/articles/2584/>`_)
to form **zOTUs and an zOTU table** (herein also referred as ASVs). 

The input is the directory that contains per-sample fastq files (**demultiplexed data**).

Pipeline final outputs are in the ``clustering_out`` directory; but per process a separate 
output directory is created (e.g. ``primersCut_out``, ``chimeraFiltered_out`` etc.).

*more to come ...*

__________________________________________________

.. _vsearchOTUs:

vsearch OTUs
============


This automated workflow is mostly based on `vsearch <https://github.com/torognes/vsearch>`_ (`Rognes et. al 2016 <https://peerj.com/articles/2584/>`_)
to form **OTUs and an OTU table**. 
The input is the directory that contains per-sample fastq files (**demultiplexed data**).

Pipeline final outputs are in the ``clustering_out`` directory; but per process a separate 
output directory is created (e.g. ``primersCut_out``, ``chimeraFiltered_out`` etc.).

.. _vsearchOTUs_defaults:

| **Default options:**
| *click on analyses step for more info*

==================================================================== =========================
Analyses step                                                        Default setting
==================================================================== =========================
:ref:`CUT PRIMERS <remove_primers>` (optional)                         --
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

__________________________________________________

.. _nextits_pipeline: 

NextITS
=======

`NextITS <https://next-its.github.io>`_ is an automated pipeline for analysing **full-length ITS** reads 
obtained via **PacBio** sequencing. 

| This pipeline implements:
| * primer trimming
| * quality filtering
| * full-length ITS region extraction
| * correction of homopolymer errors
| * chimera filtering (`get database for reference-based chimera filtering here <https://owncloud.ut.ee/owncloud/s/iaQ3i862pjwYgdy>`_)
| * recovery of sequences false-positively annotated as chimeric
| * detection of tag-switching artifacts per sequencing run
| * multiple options for sequence clustering
| * post-clustering with LULU

.. note:: 

  Please see other details here: https://next-its.github.io
  **Please note that NextITS pipeline accepts only a single primer pair**, i.e., one forward and one reverse primer in STEP_1!

.. important:: 

  NextITS requires your data and folders to be structured in a specific way (see below)! 
  Directory ``my_dir_for_NextITS`` contains ``Input`` [hard-coded requirement here] and one or multiple sequencing runs.
  In the below example, the sequencing runs [``RunID``] are named as Run1, Run2 and Run3 (but naming can be different).

  Although native NextITS requires multiplexed data as an input, the PipeCraft2's implementation **requires demultiplexed data**. So, if you have multiplexed data, then first use the DEMULTIPLEX QuickTool.
  
  In PipeCraft2, following the examples below, select ``my_dir_for_NextITS`` as a **WORKDIR**.
  


Single sequencing run
---------------------

| Select ``my_dir_for_NextITS`` as a WORKDIR in PipeCraft2.
| Directory structure for analysing a single sequencing run:

.. code-block::
   :caption: Required directory structure for NextITS

    my_dir_for_NextITS/   # SELECT THIS FOLDER AS WORKING DIRECTORY
    └── Input/
        ├── Run1/      # name here can be anything (without spaces)
        │   ├── sample1_R1.fastq.gz
        │   ├── sample1_R2.fastq.gz 
        │   ├── sample2_R1.fastq.gz
        │   └── sample2_R2.fastq.gz

Input data for this pipeline **must be demultiplexed**, if your data is multiplexed use the demultiplexer 
from **QuickTools** before running the pipeline.

.. admonition:: Sample naming

  Please avoid non-ASCII symbols in ``SampleID``,
  and do not use the period symbol (.), as it represents the wildcard character in regular expressions.
  Also, it is preferable not to start the sample name with a number.

Multiple sequencing runs
------------------------

| Select ``my_dir_for_NextITS`` as a WORKDIR in PipeCraft2.
| Directory structure for analysing multiple sequencing runs:

.. code-block::
   :caption: Required directory structure for NextITS

    my_dir_for_NextITS/   # SELECT THIS FOLDER AS WORKING DIRECTORY
    └── Input/
        ├── Run1/      # name here can be anything (without spaces)
        │   ├── Run1__sample1_R1.fastq.gz
        │   ├── Run1__sample1_R2.fastq.gz 
        │   ├── Run1__sample2_R1.fastq.gz
        │   └── Run1__sample2_R2.fastq.gz
        ├── Run2/      # name here can be anything (without spaces)
        │   ├── Run2__sample3_R1.fastq.gz
        │   ├── Run2__sample3_R2.fastq.gz
        │   ├── Run2__sample4_R1.fastq.gz
        │   └── Run2__sample4_R2.fastq.gz
        └── Run3/      # name here can be anything (without spaces)
            ├── Run3__sample5_R1.fastq.gz
            └── Run3__sample5_R2.fastq.gz

Input data for this pipeline **must be demultiplexed**, if your data is multiplexed use the demultiplexer 
from **QuickTools** before running the pipeline.

.. admonition:: Sample naming

  Please avoid non-ASCII symbols in ``RunID`` and ``SampleID``,
  and do not use the period symbol (.), as it represents the wildcard character in regular expressions.
  Also, it is preferable not to start the sample name with a number.

NextITS uses the ``SequencingRunID__SampleID`` naming convention (please note the double underscore separating ``RunID`` and ``SampleID`` parts). 
This naming scheme allows to easily trace back sequences, especially if the same sample was sequenced several times and is present in multiple sequencing runs. 
In the later steps, extracting the SampleID part and summarizing read counts for such samples is easy.



**Default settings:**

+---------------------------------------------------------------------------------------------------------------+------------------------------------+
| Analyses step                                                                                                 | Default setting                    |
+===============================================================================================================+====================================+
|| STEP 1: `QUALITY CONTROL, ARTEFACT REMOVAL <https://next-its.github.io/assets/NextITS_Workflow_Step1.webp>`_ || ``primer_mismatch`` = 2           |
||                                                                                                              || ``its_region`` = full             |
||                                                                                                              || ``qc_maxhomopolymerlen`` = 25     |
||                                                                                                              || ``qc_maxn`` = 4                   |
||                                                                                                              || ``ITSx_evalue`` = 1e-2            |
||                                                                                                              || ``ITSx_partial`` = 0              |
||                                                                                                              || ``ITSx_tax`` = all                |
||                                                                                                              || ``chimera_rescue_occurrence`` = 2 |
||                                                                                                              || ``tj f`` = 0.01                   |
||                                                                                                              || ``tj p`` = 1                      |
||                                                                                                              || ``hp`` = TRUE                     |
+---------------------------------------------------------------------------------------------------------------+------------------------------------+
|| STEP 2: `DATA AGGREGATION, CLUSTERING <https://next-its.github.io/assets/NextITS_Workflow_Step2.webp>`_      || ``otu_id`` = 0.98                 |
||                                                                                                              || ``swarm_d`` = 1                   |
||                                                                                                              || ``lulu`` = TRUE                   |
||                                                                                                              || ``unoise`` = FALSE                |
||                                                                                                              || ``otu_id_def`` = 2                |
||                                                                                                              || ``otu_qmask`` = dust              |
||                                                                                                              || ``swarm_fastidious`` = TRUE       |
||                                                                                                              || ``unoise_alpha`` = 2              |
||                                                                                                              || ``unoise_minsize`` = 8            |
||                                                                                                              || ``max_MEEP`` = 0.5                |
||                                                                                                              || ``max_chimera_score`` = 0.5       |
||                                                                                                              || ``lulu_match`` = 95               |
||                                                                                                              || ``lulu_ratio`` = 1                |
||                                                                                                              || ``lulu_ratiotype`` = min          |
||                                                                                                              || ``lulu_relcooc`` = 0.95           |
||                                                                                                              || ``lulu_maxhits`` = 0              |
+---------------------------------------------------------------------------------------------------------------+------------------------------------+


Cut primers
-----------

**Please note that NextITS pipeline accepts only a single primer pair**, i.e., one forward and one reverse primer!

================================ =========================
Setting                          Tooltip
================================ =========================
``primer_forward``               | Specify forward primer, IUPAC codes allowed
``primer_reverse``               | Specify reverse primer, IUPAC codes allowed
``primer_mismatch``              | Specify allowed number of mismatches for primers
================================ =========================

Quality filtering
-----------------

Filter sequences based on expected errors per sequence and per base, compress and correct homopolymers.

================================ =========================
Setting                          Tooltip
================================ =========================
``qc_maxee``                     | Maximum number of expected errors
``qc_maxeerate``                 | Maximum number of expected error per base
``qc_maxn``                      | Discard sequences with more than the specified number of ambiguous nucleotides (N's)
``qc_maxhomopolymerlen``         | Threshold for a homopolymer region lenght in a sequence
``hp``                           | Enable or disable homopolymer correction
================================ =========================

ITS extraction
--------------

| When performing ITS metabarcoding, it may be beneficial to trim the flanking 18S and 28S rRNA genes; because:

 - these conserved regions don't offer species-level differentiation.
 - random errors in these areas can disrupt sequence clustering.
 - chimeric breakpoints, which are common in these regions, are hard to detect in short fragments ranging from 10 to 70 bases.

 | NextITS deploys the `ITSx software <https://microbiology.se/software/itsx/>`_ (`Bengtsson-Palme et al. 2013 <https://doi.org/10.1111/2041-210X.12073>`_) for extracting the ITS sequence. 

================================ =========================
Setting                          Tooltip
================================ =========================
``its_region``                   | ITS part selector (ITS1, ITS2 or full)
``ITSx_tax``                     | Taxonomy profile for ITSx can be used to restrict the search to only taxon(s) of interest.
``ITSx_evalue``                  | E-value cutoff threshold for ITSx
``ITSx_partial``                 | Keep partial ITS sequences (specify a minimum length cutoff)
================================ =========================

Chimera filtering
-----------------

| NextITS employs a two-pronged strategy to detect chimeras: de novo and reference-based chimera filtering.
| A **reference database** for chimera filtering from full-length ITS data is `accessible here <https://owncloud.ut.ee/owncloud/s/iaQ3i862pjwYgdy>`_. This database is based on `EUKARYOME database <https://eukaryome.org>`_

Additional step in NextITS is a **"rescue" of sequences** that were initially flagged as chimeric, but are occur at least in 2 samples (which represent independent PCR reactions); 
thus are likely false-positive chimeric sequences. The chimeric sequence occurrence frequency can be edited using the --chimera_rescue_occurrence parameter.

================================ =========================
Setting                          Tooltip
================================ =========================
``chimera_database`` (optional)  | Database for reference based chimera removal (UDB)
``chimera_rescue_occurence``     | A minimum occurence of initially flagged chimeric sequence required to rescue them
================================ =========================

Tag-jump correction
-------------------

 Tag-jumps, sometimes referred to as index-switches or index cross-talk, may represent a significant concern in high-throughput sequencing (HTS) data. 
 They can cause technical cross-contamination between samples, potentially distorting estimates of community composition. 
 Here, tag-jump events are evaluated the UNCROSS2 algorithm (`Edgar 2018 <https://www.biorxiv.org/content/10.1101/400762v1>`_ ) are removed.

================================ =========================
Setting                          Tooltip
================================ =========================
``tj_f``                         | `UNCROSS <https://www.drive5.com/usearch/manual/uncross_algo.html>`_ parameter f for tag-jump filtering
``tj_p``                         | `UNCROSS <https://www.drive5.com/usearch/manual/uncross_algo.html>`_ parameter p for tag-jump filtering
================================ =========================

UNOISE denoising
----------------

 | The UNOISE algorithm (`Edgar 2016 <https://www.biorxiv.org/content/10.1101/081257v1>`_ ) focuses on error-correction (or denoising) of amplicon reads. Essentially, UNOISE operates on the principle that if a sequence with low abundance closely resembles another sequence with high abundance, the former is probably an error. This helps differentiate between true biological variation and sequencing errors. It's important to note that UNOISE was initially designed and optimized for Illumina data. Because of indel errors stemming from inaccuracies in homopolymeric regions, UNOISE might not work well with data that hasn't undergone homopolymer correction.

================================ =========================
Setting                          Tooltip
================================ =========================
``unoise``                       | Enable or disable denoising with `UNOISE <https://www.drive5.com/usearch/manual/unoise_algo.html>`_ algorithm
``unoise_alpha``                 | Alpha parameter for `UNOISE <https://www.drive5.com/usearch/manual/unoise_algo.html>`_
``unoise_minsize``               | Minimum sequence abundance
================================ =========================

Clustering
----------

NextITS supports 3 different clustering methods:

  - vsearch:
    this employs greedy clustering using a fixed sequence similarity threshold with VSEARCH (`Rognes et al., 2016, <https://peerj.com/articles/2584/>`_ );

  - swarm:
    dynamic sequence similarity threshold for clustering with SWARM (`Mahé et al., 2021, <https://academic.oup.com/bioinformatics/article/38/1/267/6318385?login=false>`_ );

  - unoise:
    creates zero-radius OTUs (zOTUs) based on the UNOISE3 algorithm (`Edgar 2016 <https://www.biorxiv.org/content/10.1101/081257v1>`_ );

================================ =========================
Setting                          Tooltip
================================ =========================
``clustering_method``            | Sequence clustering method (choose from: vsearch, swarm, unoise)
``otu_id``                       | Sequence similarity threshold
``otu_iddef``                    | Sequence similarity definition (applied to UNOISE as well)
``otu_qmask``                    | Method to mask low-complexity sequences (applied to UNOISE as well)
``swarm_d``                      | `SWARM <https://github.com/torognes/swarm>`_ clustering resolution (d)
``swarm_fastidious``             | Link nearby low-abundance swarms (fastidious option)
================================ =========================

Post-clustering with LULU
-------------------------

The purpose of LULU is to reduce the number of erroneous OTUs in OTU tables to achieve more realistic biodiversity metrics. 
By evaluating the co-occurence patterns of OTUs among samples LULU identifies OTUs that 
consistently satisfy some user selected criteria for being errors of more abundant OTUs and merges these OTUs.

================================ =========================
Setting                          Tooltip
================================ =========================
``lulu``                         | Enable or disable post-clustering curation with `lulu <https://github.com/tobiasgf/lulu>`_
``lulu_match``                   | Minimum similarity threshold
``lulu_ratio``                   | Minimum abundance ratio
``lulu_ratiotype``               | Abundance ratio type - "min" or "avg
``lulu_relcooc``                 | Relative co-occurrence
``lulu_maxhits``                 | Maximum number of hits (0 = unlimited)
================================ =========================


__________________________________________________

.. _optimotu_pipeline:

OptimOTU
========

| OptimOTU is a full metabarcoding data analysis pipeline for **paired-end Illumina data** (`arXiv:2502.10350 <https://doi.org/10.48550/arXiv.2502.10350>`_).
| OptimOTU uses taxonomically identified reference sequences to 
| determine optimal genetic distance thresholds for clustering ancestor 
| taxa into groups that best match their descendant taxa (**taxonomically aware OTU clustering**).

PipeCraft2's implementation of OptimOTU is **currently restricted to Fungi (ITS3-ITS4 and g/fITS7-ITS4 amplicons) and Metazoa COI amplicons**.

.. important::

  OptimOTU requires a specific directory structure for input data. See below.
  Note than if you are analysing one a singe sequencing run, you still need to follow the directory structure, but just 
  need to have a single directory in "01_raw" (e.g. "Run1", but you can name it as you want).

.. code-block::
   :caption: Required directory structure for OptimOTU

    my_dir_for_optimotu/   # SELECT THIS FOLDER AS WORKING DIRECTORY
    └── sequences/
        └── 01_raw/
            ├── Run1/      # name here can be anything (without spaces)
            │   ├── sample1_R1.fastq.gz
            │   ├── sample1_R2.fastq.gz
            │   ├── sample2_R1.fastq.gz
            │   └── sample2_R2.fastq.gz
            ├── Run2/      # name here can be anything (without spaces)
            │   ├── sample3_R1.fastq.gz
            │   ├── sample3_R2.fastq.gz
            │   ├── sample4_R1.fastq.gz
            │   └── sample4_R2.fastq.gz
            └── Run3/      # name here can be anything (without spaces)
                ├── sample5_R1.fastq.gz
                └── sample5_R2.fastq.gz

Output files will be saved in the ``my_dir_for_optimotu/output`` directory.
Intermediate files will be saved in the ``my_dir_for_optimotu/sequences/02_trim`` etc directories.

Target taxa and sequence orientation
------------------------------------

Specify if target taxa is fungi or metazoa, and if provided sequences are are expected to be forward, reverse or mixed orientation.

| "fwd" = all sequences are expected to be in 5'-3' orientation.
| "rev" = all sequences are expected to be in 3'-5' orientation.
| "mixed" = the orientation of seqs is expected to be mixed (5'-3' and 3'-5)
| "custom" = the orientation of different files is given in a custom sample table (see :ref:`custom_sample_table`)
| if seqs are "mixed", but using "fwd" setting, then some valid seqs (or samples) will be lost.
| if seqs are "fwd", but using "mixed", then ERROR.

+---------------------+---------------------------------------------------------------------------------+
| Setting             | Tooltip                                                                         |
+=====================+=================================================================================+
| ``target taxa``     | specify if target taxa is fungi or metazoa                                      |
+---------------------+---------------------------------------------------------------------------------+
| ``seq orientation`` | specify if provided sequences are forward (fwd), reverse (rev) or mixed (mixed) |
+---------------------+---------------------------------------------------------------------------------+



Control sequences
-----------------

Two types of control sequences are supported:

1) spike-in sequences: sequences that are added to the samples before PCR
   These sequences are expected to be present in every sample, even
   most types of negative control.
2) positive control sequences: sequences that are added to only a few specific
   positive control samples.  These sequences are expected to be present only
   in the positive control samples, and their presence in other samples is
   indicative of cross-contamination. (Either in the lab or "tag-switching").

In practice both types are treated the same by the pipeline, they are just
reported separately.

The sequences should be in a fasta file.
Specifying either or both type of control sequences is optional.

+----------------------+----------------------------------------------------------+
| Setting              | Tooltip                                                  |
+======================+==========================================================+
| ``spike in``         | (optional) specigy a file with spike-in sequences        |
+----------------------+----------------------------------------------------------+
| ``positive control`` | (optional) specify a file with positive control sequence |
+----------------------+----------------------------------------------------------+

Cut primers and trim reads
--------------------------

Cut primers and trim reads according to the specified parameters (using cutadapt).

+--------------------------+-----------------------------------------------------------------------+
| Setting                  | Tooltip                                                               |
+==========================+=======================================================================+
| ``forward primer``       | specify forward primer sequence (supports only single primer)         |
+--------------------------+-----------------------------------------------------------------------+
| ``reverse primer``       | specify reverse primer sequence (supports only single primer)         |
+--------------------------+-----------------------------------------------------------------------+
| ``max error rate``       | (maximum allowed error rate in the primer search)                     |
+--------------------------+-----------------------------------------------------------------------+
| ``truncQ_R1``            | truncate ends (3') of R1 at first base with quality score <= N        |
+--------------------------+-----------------------------------------------------------------------+
| ``truncQ_R2``            | truncate ends (3') of R2 at first base with quality score <= N        |
+--------------------------+-----------------------------------------------------------------------+
| ``min_length``           | minimum length of the trimmed sequence                                |
+--------------------------+-----------------------------------------------------------------------+
| ``cut_R1``               | remove N bases from start of R1                                       |
+--------------------------+-----------------------------------------------------------------------+
| ``cut_R2``               | remove N bases from start of R2                                       |
+--------------------------+-----------------------------------------------------------------------+
|| ``action``              || trim = trim the primers from the reads;                              |
||                         || retain = retain the primers after primer has been founds             |
+--------------------------+-----------------------------------------------------------------------+
|| ``custom_sample_table`` || custom primer trimming parameters per sample can be given as columns |
||                         || in the sample table. See example below.                              |
+--------------------------+-----------------------------------------------------------------------+


.. _custom_sample_table:

custom sample table
~~~~~~~~~~~~~~~~~~~

Example of custom primer trimming parameters per sample (**tab-delimited**): 

+--------+---------+------------------+------------------+--------+
| seqrun | samples | fastq_R1         | fastq_R2         | orient |
+--------+---------+------------------+------------------+--------+
| run1   | sample1 | sample1_R1.fq.gz | sample1_R2.fq.gz | fwd    |
+--------+---------+------------------+------------------+--------+
| run1   | sample2 | sample2_R1.fq.gz | sample2_R2.fq.gz | fwd    |
+--------+---------+------------------+------------------+--------+
| run2   | sample3 | sample3_R1.fq.gz | sample3_R2.fq.gz | rev    |
+--------+---------+------------------+------------------+--------+
| run2   | sample4 | sample4_R1.fq.gz | sample4_R2.fq.gz | rev    |
+--------+---------+------------------+------------------+--------+
| run3   | sample5 | sample5_R1.fq.gz | sample5_R2.fq.gz | mixed  |
+--------+---------+------------------+------------------+--------+


Quality filtering
-----------------

Quality filtering settings; performed using DADA2. Sequences with ambiguous nucleotides (N's) are discarded. 

+--------------+--------------------------------------------------------------------------------------+
| Setting      | Tooltip                                                                              |
+==============+======================================================================================+
| ``maxEE_R1`` | discard sequences with more than the specified number of expected errors in R1 reads |
+--------------+--------------------------------------------------------------------------------------+
| ``maxEE_R2`` | discard sequences with more than the specified number of expected errors in R2 reads |
+--------------+--------------------------------------------------------------------------------------+


Denoising and merging paired-end reads 
--------------------------------------

There are no adjustable setting for denoising. 
The denoising steps are performed using the DADA2 package (Callahan et al. 2016). 
Error profiles are then learned separately for each
sequencing run, read, and orientation using the learnErrors() function. Sequences with binned quality
scores, as produced by newer Illumina sequencers, are automatically detected, and the error model is adjusted
accordingly. Denoising is then performed using the dada() function, and read pairs are merged using the
mergePairs() function.


Chimera filtering
-----------------

Chimera filtering is performed using the consensus algorithm implemented in DADA2's isBimeraDenovoTable() function.
Additional database provided in the PROTAX CLASSIFICATION step (``with_outgroup`` file) is used for reference-based chimera filtering (vsearch --uchime_ref).


Filter tag-jumps
----------------

Filter potential cases of tag-switching with UNCROSS2 algorithm (Edgar 2018).

+--------------+----------------------------------------------------------------------------------------+
| Setting      | Tooltip                                                                                |
+==============+========================================================================================+
|| ``f value`` || f-parameter of UNCROSS2, which defines the expected tag-jumps rate. Default is 0.03   |
||             || (equivalent to 3%). A higher value enforces stricter filtering                        |
+--------------+----------------------------------------------------------------------------------------+
|| ``p value`` || p-parameter, which controls the severity of tag-jump removal. It adjusts the exponent |
||             || in the UNCROSS formula. Default is 1. Opt for 0.5 or 0.3 to steepen the curve         |
+--------------+----------------------------------------------------------------------------------------+


Amplicon model setting
----------------------

+---------------------+---------+
| Setting             | Tooltip |
+=====================+=========+
| ``model type``      |         |
+---------------------+---------+
| ``model file``      |         |
+---------------------+---------+
| ``numt filter``     |         |
+---------------------+---------+
| ``max model start`` |         |
+---------------------+---------+
| ``min model end``   |         |
+---------------------+---------+
| ``min model score`` |         |
+---------------------+---------+


ProTAX classification
---------------------

+--------------------+--------------------------------------------------------------------------+
| Setting            | Tooltip                                                                  |
+====================+==========================================================================+
|| ``location``      || directory where protax is located. For fungi, default is protaxFungi    |
||                   || and for protaxAnimal for metazoa (included in the PipeCraft2 container) |
+--------------------+--------------------------------------------------------------------------+
|| ``with outgroup`` || additional database which contains also outgroup (non-target)           |
||                   || sequences from the same locus                                           |
+--------------------+--------------------------------------------------------------------------+

Clustering
----------

+-------------------------+--------------------------------------------------------------------+
| Setting                 | Tooltip                                                            |
+=========================+====================================================================+
|| ``cluster thresholds`` || select file with clustering thresholds. Default is pre-calculated |
||                        || thresholds for Fungi (UNITE v x.x)                                |
+-------------------------+--------------------------------------------------------------------+


