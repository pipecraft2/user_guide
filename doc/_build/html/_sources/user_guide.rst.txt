.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text

.. |main_interface| image:: _static/main_interface.png
  :width: 2000
  :alt: Alternative text

.. |asv_main| image:: _static/asv_main.png
  :width: 1500
  :alt: Alternative text

.. |otu_main| image:: _static/otu_main.png
  :width: 1500
  :alt: Alternative text

.. |console| image:: _static/console.png
  :width: 1500
  :alt: Alternative text

.. meta::
    :description lang=en:
        PipeCraft manual. User guide for PipeCraft

|PipeCraft2_logo|
  `github <https://github.com/SuvalineVana/pipecraft>`_

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red

==========
User guide
==========

.. contents:: Contents
   :depth: 3

____________________________________________________

.. _interface: 

The interface
==============

The startup panel:

*(click on the image for enlargement)*
|main_interface|

____________________________________________________

Docker images 
==============

.. |pulling_image| image:: _static/pulling_image.png
  :width: 280
  :alt: Alternative text


Initial PipeCraft installation does not contain any software for sequence data processing. 
All the processes are run through `docker <https://www.docker.com/>`_, where the PipeCraft's simply GUI mediates the 
information exchange. Therefore, whenever a process is initiated for the **first time**, 
a relevant Docker image (contains required software for the analyses step) will be pulled from `Docker Hub <https://hub.docker.com/u/pipecraft>`_.

Example: running DEMULTIPLEXING for the first time |pulling_image|

Thus working **Internet connection** is initially required. Once the Docker images are pulled, PipeCraft can work without an Internet connection. 

:ref:`Docker images <containers>` vary in size, and the speed of the first process is extended by the docker image download time.
 
____________________________________________________

Save workflow
==============

Once the workflow settings are selected, save the workflow by pressin ``SAVE WORKFLOW`` button on the :ref:`right-ribbon <interface>`.
For saving, working directory ( ``SELECT WORKDIR`` ) does not have to be selected. 

.. important::

 When **saiving workflow** settings in **Linux**, specify the file extension as **JSON** (e.g. my_16S_ASVs_pipe.JSON).
 When trying to load the workflow, only .JSON files will be permitted as input. *Windows and Mac OS automatically extend files as JSON (so you may just save "my_16S_ASVs_pipe").*

____________________________________________________

Load workflow
==============

.. note ::

 Prior loading the workflow, make sure that the saved workflow configuration has a .JSON extension. 

Press the ``LOAD WORKFLOW`` button on the :ref:`right-ribbon <interface>` and select appropriate JSON file.
The configuration will be loaded; ``SELECT WORKDIR`` and run analyses.

____________________________________________________

.. _qualitycheck:

Quality and basic statistics screening of the data
==================================================

.. |multiQC_main| image:: _static/multiQC_main.png
  :width: 1000
  :alt: Alternative text

.. |multiQC_1-3| image:: _static/multiQC_1-3.png
  :width: 550
  :alt: Alternative text

.. |multiQC_view_report| image:: _static/multiQC_view_report.png
  :width: 550
  :alt: Alternative text


Quality and basic statistics screening of the data can be done via ``QualityCheck`` panel. 
QualityCheck panel implements `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ and `MultiQC <https://multiqc.info/>`_ to screen the input **fastq** files. 

|multiQC_main|

| 

**To start:** 

 1. **Select folder** (a working directory) which contains only **fastq** (fastq/fq) files that you aim to inspect. 
 2. Press ``CREATE REPORT`` to start MultiQC 
 3. "LOADING ..." will be displayed while the report is being generated

|multiQC_1-3|

|multiQC_view_report|

 4. Click ``VIEW REPORT``. A html file (multiqc_report.html) will open in your default web browser.
    
    *If the summary does not open, check your working floder for the presence of* **multiqc_report.html** *and try to open with some other web browser.*
    *Something went wrong if the file multiqc_report.html* **does not exist** *(may fail when maximum number of fastq files in the folder is extremely large, >10 000).*

 5. Check out  `"using MultiQC reports" <https://multiqc.info/docs/#using-multiqc-reports>`_ in MultiQC web page.
   
.. note::

 Note that '_fastqc.zip' and '_fastqc.html' are generated for each fastq file in the **'quality_check'** directory. These are summarized in **multiqc_report.html**, 
 so you may delete all individual '_fastqc.zip' and '_fastqc.html' files.
 
| 

____________________________________________________

Select workdir and run analyses
===============================

1. Open your working directory by pressing the ``SELECT WORKDIR`` button. E.g., if working with **FASTQ** files,
then be sure that the working directory contains **only relevant FASTQ files** because the selected process will be 
applied to all FASTQ files in the working directory!

.. note::

 When using Windows OS, the selection window might not display the files while browsing through the directories. 

After selecting a working directory, PipeCraft needs you to specify if 
the working directory consists of 

 * multiplexed or demultiplexed data
 * the data is paired-end or single-end
 * and the extension of the data (fastq or fasta)

| ``multiplexed`` --> only one file (or a pair of files, R1 and R2) per sequencing data (library)
| ``demultiplexed`` --> multiple per-sample sequencing files per library
| ``paired-end data`` --> such as data from Illumina or MGI-Tech platforms (R1 and R2 files). :red:`Be sure to have **R1** and **R2** strings in the paired-end files (not simply _1 and _2)`
| ``single-end data`` --> such as data from PacBio, or assembled paired-end data (single file per library or per sample)

2. Select :ref:`ASV <asvpipe>` or :ref:`OTU <otupipe>` workflow panel or press ``ADD STEP`` button
to select relevant :ref:`step <panels>` [or **load the PipeCraft settings file**]; 
edit settings if needed (**SAVE the settings for later use**) and **start
running the analyses** by pressing the ``RUN WORKFLOW`` button.


.. note::

 **Step-by-step analyses**: after ``RUN WORKFLOW`` is finished, then press ``SELECT WORKDIR`` to specify inputs for the next process

.. note::

 The **output files will be overwritten** if running the same 
 analysis step **multiple times in the same working directory**

3. Each process creates a separate output directory with the processed files
inside the selected working directory. 
**README** file about the process and **sequence count summary** statistics are included in the output directory.

____________________________________________________

FULL PIPELINE PANELS
====================

.. |asvoff| image:: _static/ASV_off.png
  :width: 50
  :alt: Alternative text

.. |asvon| image:: _static/ASV_on.png
  :width: 50
  :alt: Alternative text

.. _asvpipe:

ASVs workflow panel (with `DADA2 <https://benjjneb.github.io/dada2/index.html>`_)
----------------------------------------------------------------------------------

.. note::
  Current ASVs workflow supports only **PAIRED-END** reads! 
  Working directory must contain paired-end reads for at **least 2 samples**.

|asv_main|

ASV workflow is active (green icon) |asvon|
; ASV workflow is off |asvoff| 

This automated workflow is based on the `DADA2 tutorial <https://benjjneb.github.io/dada2/tutorial.html>`_ 
 | Note that ``demultiplexing``, ``reorienting``, and ``primer removal`` steps are optional and do not represent parts from the DADA2 tutorial. Nevertheless, it is advisable to :ref:`reorient <reorient>` your reads (to 5'-3') and :ref:`remove primers <remove_primers>` before proceeding with ASV generation with DADA2.

| The official DADA2 manual is available `here <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_
 
.. _dada2_defaults:

**Default options:**

================================================== =========================
Analyses step                                      Default setting
================================================== =========================
:ref:`DEMULTIPLEX <demux>` (optional)              --
:ref:`REORIENT <reorient>` (optional)              --
:ref:`REMOVE PRIMERS <remove_primers>` (optional)  --
:ref:`QUALITY FILTERING <dada2_qual_filt>`         | ``read_R1`` = \\.R1
                                                   | ``read_R2`` = \\.R2
                                                   | ``samp_ID`` = \\.
                                                   | ``maxEE`` = 2
                                                   | ``maxN`` = 0
                                                   | ``minLen`` = 20
                                                   | ``truncQ`` = 2
                                                   | ``truncLen`` = 0
                                                   | ``maxLen`` = 9999
                                                   | ``minQ`` = 2
                                                   | ``matchIDs`` = TRUE
:ref:`DENOISE <dada2_denoise>`                     | ``pool`` = FALSE
                                                   | ``selfConsist`` = FASLE
                                                   | ``qualityType`` = Auto
:ref:`MERGE PAIRED-END READS <dada2_merge_pairs>`  | ``minOverlap`` = 12
                                                   | ``maxMismatch`` = 0
                                                   | ``trimOverhang`` = FALSE
                                                   | ``justConcatenate`` = FALSE
:ref:`CHIMERA FILTERING <dada2_chimeras>`          | ``method`` = consensus
:ref:`ASSIGN TAXONOMY <dada2_taxonomy>` (optional) | ``minBoot`` = 50
                                                   | ``tryRC`` = FALSE
                                                   | ``dada2 database`` = select a database
================================================== =========================

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
``samp_ID``                      | applies only for **paired-end** data. 
                                 | Identifyer string that separates the sample name from redundant 
                                 | charachters (e.g. file name = sample1.R1.fastq, then 
                                 | underscore '\\.' would be the 'identifier string' (sample name = sampl84)); 
                                 | note that backslash is only needed to 
                                 | escape dot regex (e.g. when file name = sample1_R1.fastq then specify as '_')
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _otupipe:

OTUs workflow panel
--------------------

.. |otuoff| image:: _static/OTU_off.png
  :width: 50
  :alt: Alternative text

.. |otuon| image:: _static/OTU_on.png
  :width: 50
  :alt: Alternative text

.. note::
  This OTU workflow works with paired-end (e.g. Illumina, MGI-Tech) as well as single-end reads (e.g. PacBio, assembled Illumina reads)

|otu_main|

OTU workflow is active (green icon) |otuon|
; OTU workflow is off |otuoff| 

This automated workflow is mostly based on `vsearch <https://github.com/torognes/vsearch>`_ (`Rognes et. al 2016 <https://peerj.com/articles/2584/>`_) [`manual <_static/vsearch_2.18.0_manual.pdf>`_]
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
:ref:`MERGE READS <vsearch_merge>`                                   | ``min_overlap`` = 12
                                                                     | ``min_length`` = 32
                                                                     | ``allow_merge_stagger`` = TRUE 
                                                                     | ``include only R1`` = FALSE 
                                                                     | ``max_diffs`` = 20
                                                                     | ``max_Ns`` = 0
                                                                     | ``max_len`` = 600
                                                                     | ``keep_disjoined`` = FALSE 
                                                                     | ``fastq_qmax`` = 41
:ref:`QUALITY FILTERING with vsearch <vsearch_qfilt>`                | ``maxEE`` = 1
                                                                     | ``maxN`` = 0
                                                                     | ``minLen`` = 32
                                                                     | ``max_length`` = undefined
                                                                     | ``qmax`` = 41
                                                                     | ``qmin`` = 0
                                                                     | ``maxee_rate`` = undefined
                                                                     | ``minsize`` = 1
:ref:`CHIMERA FILTERING with vsearch <chimFilt_vsearch>`             | ``pre_cluster`` = 0.98
                                                                     | ``min_unique_size`` = 1
                                                                     | ``denovo`` = TRUE 
                                                                     | ``reference_based`` = undefined
                                                                     | ``abundance_skew`` = 2
                                                                     | ``min_h`` = 0.28
:ref:`ITS Extractor <itsextractor>` (optional)                       | ``organisms`` = Fungi 
                                                                     | ``regions`` = all
                                                                     | ``partial`` = 50
                                                                     | ``e_value`` = 1e-5
                                                                     | ``scores`` = 0
                                                                     | ``domains`` = 2
                                                                     | ``complement`` = TRUE 
                                                                     | ``only_full`` = FALSE
                                                                     | ``truncate`` = TRUE 
:ref:`CLUSTERING with vsearch <clustering_vsearch>`                  | ``OTU_type`` = centroid
                                                                     | ``similarity_threshold`` = 0.97
                                                                     | ``strands`` = both
                                                                     | ``min_OTU_size`` = 2
                                                                     | ``similarity_type`` = 2
                                                                     | ``sequence_sorting`` = cluster_size
                                                                     | ``centroid_type`` = similarity
                                                                     | ``max_hits`` = 1
                                                                     | ``relabel`` = sha1
                                                                     | ``mask`` = dust
                                                                     | ``dbmask`` = dust
                                                                     | ``output_UC`` = FALSE
:ref:`ASSIGN TAXONOMY with BLAST <assign_taxonomy_blast>` (optional) | ``database_file`` = select a database
                                                                     | ``task`` = blastn
                                                                     | ``strands`` = both
==================================================================== =========================

____________________________________________________

.. _panels:

ANALYSES PANELS
===============

.. _demux:

DEMULTIPLEX
------------

If data is **multiplexed, the first step would be demultiplexing** (using `cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_ (`Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_)).
This is done based on the user specified :ref:`indexes file <indexes>`, which includes molecular identifier sequences (so called indexes/tags/barcodes) per sample. 
Note that reverse complementary matches will also be searched. 

| **Fastq/fasta** formatted paired-end and single-end data are supported.
| **Outputs** are fastq/fasta files per sample in ``demultiplexed_out`` directory. Indexes are **truncated** from the sequences. 
| Samples get ``.R1`` and ``.R2`` read identifiers (relevant for :ref:`DADA2 QUALITY FILTERING <dada2_qual_filt>` ).
| **unknown.fastq** file(s) contain sequences where specified index combinations were not found. 

.. note:: 

  If found, sequences with any index combination will be outputted **when using paired indexes**. 
  That means, if, for example, your sample_1 is indexed with *indexFwd_1-indexRev_1* and 
  sample_2 with *indexFwd_2-indexRev_2*, then files with *indexFwd_1-indexRev_2* and *indexFwd_2-indexRev_1*
  are also written (although latter index combinations were not used in the lab to index any sample [i.e. represent tag-switches]). 
  Simply remove those files if not needed or use to estimate tag-switching error if relevant. 

.. _demux_settings:

================================ =========================
Setting                          Tooltip
================================ =========================
``index file``                   | select your fasta formatted indexes file for demultiplexing (:ref:`see guide here <indexes>`), 
                                 | where fasta headers are sample names, and sequences are sample 
                                 | specific index or index combination 
``index mismatch``               | allowed mismatches during the index search
``overlap``                      | number of overlap bases with the index
                                 | Recommended overlap is the maximum length of the index for 
                                 | confident sequence assignments to samples
``min seq length``               | minimum length of the output sequence
``no indels``                    | do not allow insertions or deletions is primer search. 
                                 | Mismatches are the only type of errors accounted in the error rate parameter
================================ =========================


.. note::

 Heterogenity spacers or any redundant base pairs attached to index sequences do not affect demultiplexing. Indexes are trimmed from the best matching position.

.. _indexes:

Indexes file example (fasta formatted)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note::
  Only **IUPAC codes** are allowed.

1. **Demultiplexing using single indexes**:

 | >sample1
 | AGCTGCACCTAA
 | >sample2
 | AGCTGTCAAGCT
 | >sample3
 | AGCTTCGACAGT
 | >sample4
 | AGGCTCCATGTA
 | >sample5
 | AGGCTTACGTGT
 | >sample6
 | AGGTACGCAATT

2. **Demultiplexing using dual (paired) indexes:**

.. note::
 **IMPORTANT!** reverse indexes will be automatically oriented to 5'-3' (for the search); so you can simply copy-paste the indexes from your lab protocol.


| >sample1
| AGCTGCACCTAA...AGCTGCACCTAA
| >sample2
| AGCTGTCAAGCT...AGCTGTCAAGCT
| >sample3
| AGCTTCGACAGT...AGCTTCGACAGT
| >sample4
| AGGCTCCATGTA...AGGCTCCATGTA
| >sample5
| AGGCTTACGTGT...AGGCTTACGTGT
| >sample6
| AGGTACGCAATT...AGGTACGCAATT

.. note::
 Anchored indexes (https://cutadapt.readthedocs.io/en/stable/guide.html#anchored-5adapters) with ^ symbol are **not supported** in PipeCraft demultiplex GUI panel. 

 DO NOT USE, e.g. 

 | >sample1
 | ^AGCTGCACCTAA
 | 
 | >sample1
 | ^AGCTGCACCTAA...AGCTGCACCTAA

|

How to compose indexes.fasta 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In Excel (or any alternative program); 
first column represents sample names,
second (and third) column represent indexes (or index combinations) per sample:

Exaples::

     sample1	AGCTGCACCTAA
     sample2	AGCTGTCAAGCT
     sample3	AGCTTCGACAGT 
     sample4	AGGCTCCATGTA
     sample5	AGGCTTACGTGT
     sample6	AGGTACGCAATT

or ::

     sample1	AGCTGCACCTAA	AGCTGCACCTAA
     sample2	AGCTGTCAAGCT	AGCTGTCAAGCT
     sample3	AGCTTCGACAGT	AGCTTCGACAGT
     sample4	AGGCTCCATGTA	AGGCTCCATGTA
     sample5	AGGCTTACGTGT	AGGCTTACGTGT
     sample6	AGGTACGCAATT	AGGTACGCAATT

Copy those two (or three) columns to text editor that support regular expressions, such as NotePad++ or Sublime Text.
If using **PAIRED** indexes (three columns), proceed to bullet no. 5

* single-end indexes:

  #. Open 'find & replace'
     Find ^   (which denotes the beginning of each line).
     Replace with >  (and DELETE THE LAST > in the beginning of empty row).

  #. Find \\t   (which denotes tab).
     Replace with \\n   (which denotes the new line).

     **FASTA FORMATTED (single-end indexes) indexes.fasta file is ready; SAVE the file.**


* Only for paired-indexes:

  #. Open 'find & replace':
     Find ^   (denotes the beginning of each line);
     replace with >  (and DELETE THE LAST > in the beginning of empty row).

  #. Find .*\\K\\t (which captures the second tab);
     replace with ... (to mark the linked paired-indexes). 

  #. Find \\t (denotes the tab);
     replace with \\n (denotes the new line).

     **FASTA FORMATTED (paired indexes) indexes.fasta file is ready; SAVE the file.**

____________________________________________________

.. _reorient:

REORIENT
--------

Sequences are often (if not always) in both, 5'-3' and 3'-5', orientations in the raw sequencing data sets. 
If the data still contains PCR primers that were used to generate amplicons, 
then by specifying these PCR primers, this panel will perform sequence reorientation 
of all sequences. 

For reorienting, 
first the forward primer will be searched (using `fqgrep <https://github.com/indraniel/fqgrep>`_)  
and if detected then the read is considered as forward complementary (5'-3'). 
Then the reverse primer will be searched (using `fqgrep <https://github.com/indraniel/fqgrep>`_) 
from the same input data and if detected, then the read is considered to be in 
reverse complementary orientation (3'-5'). Latter reads will be transformed to 5'-3' 
orientation and merged with other 5'-3' reads. 
Note that for paired-end data, R1 files will be reoriented to 5'-3' 
but R2 reads will be reoriented to 3'-5' in order to merge paired-end reads.

At least one of the PCR primers must be found in the sequence. 
For example, read will be recorded if forward primer was found even 
though reverse primer was not found (and vice versa). 
**Sequence is discarded if none of the PCR primers are found.** 

Sequences that contain **multiple forward or reverse primers (multi-primer artefacts) 
are discarded** as it is highly likely that these are chimeric sequences. 
Reorienting sequences **will not remove** primer strings from the sequences. 

.. note::

 For single-end data, sequences will be reoriented also during 
 the 'cut primers' process (see below); therefore this step may be skipped
 when working with single-end data (such as data from PacBio machines OR already assembled paired-end data).

Reorienting reads may be relevant for generating ASVs with DADA2 
as reverse complement sequences will represent separate ASVs. 
In the clustering step of an OTU pipeline, both strands of the sequences can be compared prior 
forming OTUs; thus this step may be skipped in the OTU pipeline. 

Supported file formats for paired-end input data are only **fastq**,
but also **fasta** for single-end data.
**Outputs** are fastq/fasta files in ``reoriented_out`` directory. 
Primers are **not truncated** from the sequences; this can be done using :ref:`CUT PRIMER panel <remove_primers>`

================================ =========================
Setting                          Tooltip
================================ =========================
``mismatches``                   | allowed mismatches in the primer search
``forward_primers``              | specify forward primer **(5'-3')**; IUPAC codes allowed; 
                                 | add up to 13 primers
``reverse_primers``              | specify reverse primer **(3'-5')**; IUPAC codes allowed; 
                                 | add up to 13 primers
================================ =========================

____________________________________________________

.. _remove_primers:

CUT PRIMERS
-----------

If the input data contains PCR primers (or e.g. adapters), these can be removed in the ``CUT PRIMERS`` panel.
CUT PRIMERS processes mostly relies on `cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_ (`Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_). 

For generating OTUs or ASVs, it is recommended to truncate the primers from the reads 
(unless ITS Extractor is used later to remove flanking primer binding regions from ITS1/ITS2/full ITS). 
Sequences where PCR primer strings were not detected are discarded by default (but stored in 'untrimmed' directory). 
Reverse complementary search of the primers in the sequences is also performed. 
Thus, primers are clipped from both 5'-3' and 3'-5' oriented reads. However, note that **paired-end reads will not be reoriented** to 5'-3' during this process, 
but **single-end reads will be reoriented** to 5'-3' (thus no extra reorient step needed for single-end data).

.. note::

 For paired-end data, the **seqs_to_keep option should be left as default ('keep_all')**. This will output sequences where at least one primer has been clipped. 
 'keep_only_linked' option outputs only sequences where both the forward and reverse primers are found (i.e. 5'-forward…reverse-3'). 
 'keep_only_linked' may be used for single-end data to keep only **full-length amplicons**.

| **Fastq**/**fasta** formatted paired-end and single-end data are supported.
| **Outputs** are fastq/fasta files in ``primersCut_out`` directory. Primers are **truncated** from the sequences. 

================================ =========================
Setting                          Tooltip
================================ =========================
``forward primers``              | specify forward primer **(5'-3')**; IUPAC codes allowed; 
                                 | add up to 13 primers
``reverse primers``              | specify reverse primer **(3'-5')**; IUPAC codes allowed; 
                                 | add up to 13 primers
``mismatches``                   | allowed mismatches in the primer search
``min overlap``                  | number of overlap bases with the primer sequence. 
                                 | Partial matches are allowed, but short matches may occur by chance, 
                                 | leading to erroneously clipped bases. 
                                 | Specifying higher overlap than the length of primer sequnce 
                                 | will still clip the primer (e.g. primer length is 22 bp, 
                                 | but overlap is specified as 25 - this does not affect the 
                                 | identification and clipping of the primer as long as the match is 
                                 | in the specified mismatch error range)
``seqs to keep``                 | keep sequences where at least one primer was found (fwd or rev); 
                                 | recommended when cutting primers from paired-end data (unassembled), 
                                 | when individual R1 or R2 read lenghts are shorther than the expected 
                                 | amplicon length. 'keep_only_linked' = keep sequences if primers are found 
                                 | in both ends (fwd…rev); discards the read if both primers were not found 
                                 | in this read
``pair filter``                  | **applies only for paired-end data.**
                                 | 'both', means that a read is discarded only if both, corresponding R1 and R2,
                                 | reads  do not contain primer strings (i.e. a read is kept if R1 contains 
                                 | primer string, but no primer string found in R2 read). Option 'any' discards 
                                 | the read if primers are not found in both, R1 and R2 reads
``min seq length``               | minimum length of the output sequence
``no indels``                    | do not allow insertions or deletions is primer search. Mismatches are the 
                                 | only type of errprs accounted in the error rate parameter
================================ =========================

____________________________________________________

|

.. _qual_filt:

QUALITY FILTERING
------------------

Quality filter and trim sequences.

| **Fastq** formatted paired-end and single-end data are supported.
| **Outputs** are fastq files in ``qualFiltered_out`` directory.

.. _vsearch_qfilt:

`vsearch <https://github.com/torognes/vsearch>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

================================ =========================
**vsearch** setting              Tooltip
================================ =========================
``maxEE``                        | maximum number of expected errors per sequence (`see here <https://drive5.com/usearch/manual/exp_errs.html>`_). 
                                 | Sequences with higher error rates will be discarded
``maxN``                         | discard sequences with more than the specified number of Ns
``minLen``                       | minimum length of the filtered output sequence
``max_length``                   | discard sequences with more than the specified number of bases
``qmax``                         | specify the maximum quality score accepted when reading FASTQ files. 
                                 | The default is 41, which is usual for recent Sanger/Illumina 1.8+ files. 
                                 | **For PacBio data use 93**
``qmin``                         | the minimum quality score accepted for FASTQ files. The default is 0, which is 
                                 | usual for recent Sanger/Illumina 1.8+ files. 
                                 | Older formats may use scores between -5 and 2
``maxee_rate``                   | discard sequences with more than the specified number of expected errors per base
``minsize``                      | discard sequences with an abundance lower than the specified value
================================ =========================

| 

.. _trimmomatic_qfilt:

`trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

================================ =========================
**trimmomatic** setting          Tooltip
================================ =========================
``window_size``                  | the number of bases to average base qualities
                                 | Starts scanning at the 5'-end of a sequence and trimms the read once the 
                                 | average required quality (required_qual) within the window size falls 
                                 | below the threshold
``required_quality``             | the average quality required for selected window size
``min_length``                   | minimum length of the filtered output sequence
``leading_qual_threshold``       | quality score threshold to remove low quality bases from the beginning of the read. 
                                 | As long as a base has a value below this threshold the base is removed and 
                                 | the next base will be investigated
``trailing_qual_threshold``      | quality score threshold to remove low quality bases from the end of the read. 
                                 | As long as a base has a value below this threshold the base is removed and 
                                 | the next base will be investigated
``phred``                        | phred quality scored encoding. 
                                 | Use phred64 if working with data from older Illumina (Solexa) machines
================================ =========================


| 

.. _fastp_qfilt:

`fastp <https://github.com/OpenGene/fastp>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

================================ =========================
**fastp** setting                Tooltip
================================ =========================
``window_size``                  | the window size for calculating mean quality
``required_qual``                | the mean quality requirement per sliding window (window_size)
``min_qual``                     | the quality value that a base is qualified. Default 15 means 
                                 | phred quality >=Q15 is qualified
``min_qual_thresh``              | how many percents of bases are allowed to be unqualified (0-100)
``maxNs``                        | discard sequences with more than the specified number of Ns
``min_length``                   | minimum length of the filtered output sequence. Shorter sequences are discarded
``max_length``                   | reads longer than 'max length' will be discarded, default 0 means no limitation
``trunc_length``                 | truncate sequences to specified length. Shorter sequences are discarded; 
                                 | thus check that 'min length' setting is lower than 'trunc length'
``aver_qual``                    | if one read's average quality score <'aver_qual', then this read/pair is discarded. 
                                 | Default 0 means no requirement
``low_complexity_filter``        | enables low complexity filter and specify the threshold for low complexity filter. 
                                 | The complexity is defined as the percentage of base that is different from its 
                                 | next base (base[i] != base[i+1]). 
                                 | E.g. vaule 30 means then 30% complexity is required. 
                                 | Not specified = filter not applied
``cores``                        | number of cores to use
================================ =========================

| 

.. _dada2_qfilt:

`DADA2 <https://github.com/benjjneb/dada2>`_ ('filterAndTrim' function)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

================================ =========================
**DADA2** setting                Tooltip
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
``samp_ID``                      | applies only for **paired-end** data. 
                                 | Identifyer string that separates the sample name from redundant 
                                 | charachters (e.g. file name = sample1.R1.fastq, then 
                                 | underscore '\\.' would be the 'identifier string' (sample name = sampl84)); 
                                 | note that backslash is only needed to 
                                 | escape dot regex (e.g. when file name = sample1_R1.fastq then specify as '_')
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
                                 | after truncation, reads contain a quality score below minQ will be discarded
================================ =========================

____________________________________________________

| 

.. _merge_pairs:

ASSEMBLE PAIRED-END reads 
--------------------------

Assemble paired-end sequences (such as those from Illumina or MGI-Tech platforms). 

``include_only_R1`` represents additional in-built module. If TRUE, 
unassembled R1 reads will be included to the set of assembled reads per sample. 
This may be relevant when working with e.g. ITS2 sequences, because the ITS2 region in some 
taxa is too long for paired-end assembly using current short-read sequencing technology. 
Therefore longer ITS2 amplicon sequences are discarded completely after the assembly process. 
Thus, including also unassembled R1 reads (``include_only_R1`` = TRUE), partial ITS2 sequences for 
these taxa will be represented in the final output. But when using :ref:`ITSx <itsextractor>`  
, keep ``only_full`` = FALSE and include ``partial`` = 50.

**Fastq** formatted paired-end data is supported.
**Outputs** are fastq files in ``assembled_out`` directory.


.. _vsearch_merge:

`vsearch <https://github.com/torognes/vsearch>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

================================ =========================
Setting                          Tooltip
================================ =========================
``min_overlap``                  | minimum overlap between the merged reads
``min_length``                   | minimum length of the merged sequence
``allow_merge_stagger``          | allow to merge staggered read pairs. Staggered pairs are pairs 
                                 | where the 3' end of the reverse read has an overhang to the left 
                                 | of the 5' end of the forward read. This situation can occur when a 
                                 | very short fragment is sequenced
``include_only_R1``              | include unassembled R1 reads to the set of assembled reads per sample
``max_diffs``                    | the maximum number of non-matching nucleotides allowed in the overlap region
``max_Ns``                       | discard sequences with more than the specified number of Ns
``max_len``                      | maximum length of the merged sequence
``keep_disjoined``               | output reads that were not merged into separate FASTQ files
``fastq_qmax``                   | maximum quality score accepted when reading FASTQ files. 
                                 | The default is 41, which is usual for recent Sanger/Illumina 1.8+ files
================================ =========================

|


.. _dada2_merge:

`DADA2 <https://github.com/benjjneb/dada2>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. important::

  Here, dada2 will perform denoising (function 'dada') before assembling paired-end data. 
  Because of that, input sequences (in **fastq** format) must consist of 
  only A/T/C/Gs. 

================================ =========================
Setting                          Tooltip
================================ =========================
``read_R1``                      | identifyer string that is common for all R1 reads 
                                 | (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. 
                                 | Note that backslash is only needed to escape dot regex; e.g. 
                                 | when all R1 files have '_R1' string, then enter '_R1'.)
``read_R2``                      | identifyer string that is common for all R2 reads 
                                 | (e.g. when all R2 files have '.R2' string, then enter '\\.R2'. 
                                 | Note that backslash is only needed to escape dot regex; e.g. 
                                 | when all R2 files have '_R1' string, then enter '_R2'.)
``samp_ID``                      | identifyer string that separates the sample name from redundant 
                                 | charachters (e.g. file name = sample1.R1.fastq, then 
                                 | underscore '\\.' would be the 'identifier string' (sample name = sampl84)); 
                                 | note that backslash is only needed to escape dot regex
                                 | (e.g. when file name = sample1_R1.fastq then specify as '_')
``minOverlap``                   | the minimum length of the overlap required for merging the forward and 
                                 | reverse reads
``maxMismatch``                  | the maximum mismatches allowed in the overlap region
``trimOverhang``                 | if TRUE, overhangs in the alignment between the forwards and reverse read are  
                                 | trimmed off. Overhangs are when the reverse read extends past the start of 
                                 | the forward read, and vice-versa, as can happen when reads are longer than the 
                                 | amplicon and read into the other-direction primer region
``justConcatenate``              | if TRUE, the forward and reverse-complemented reverse read are concatenated  
                                 | rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them
``pool``                         | denoising setting. If TRUE, the algorithm will pool together all samples 
                                 | prior to sample inference. Pooling improves the detection of rare variants, 
                                 | but is computationally more expensive. 
                                 | If pool = 'pseudo', the algorithm will perform pseudo-pooling between  
                                 | individually processed samples.
``selfConsist``                  | denoising setting. If TRUE, the algorithm will alternate between sample 
                                 | inference and error rate estimation until convergence
``qualityType``                  | 'Auto' means to attempt to auto-detect the fastq quality encoding. 
                                 | This may fail for PacBio files with uniformly high quality scores, 
                                 | in which case use 'FastqQuality'
================================ =========================


.. _chimFilt:

____________________________________________________

|

CHIMERA FILTERING
-----------------

Perform de-novo and or reference database based chimera filtering. 

Chimera filtering is performed by **sample-wise approach** (i.e. each sample (input file) is treated separately). 

| **Fastq/fasta** formatted single-end data is supported [fastq inputs will be converted to fasta].
| **Outputs** are fasta files in ``chimera_Filtered_out`` directory.

.. _chimFilt_vsearch:

`vsearch <https://github.com/torognes/vsearch>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

================================ =========================
Setting                          Tooltip
================================ =========================
``pre_cluster``                  | identity percentage when performing 'pre-clustering' with --cluster_size 
                                 | for denovo chimera filtering with --uchime_denovo
``min_unique_size``              | minimum amount of a unique sequences in a fasta file. If value = 1, then 
                                 | no sequences are discarded after dereplication; if value = 2, then sequences,
                                 | which are represented only once in a given file are discarded; and so on
``denovo``                       | if TRUE, then perform denovo chimera filtering with --uchime_denovo
``reference_based``              | perform reference database based chimera filtering with --uchime_ref. 
                                 | Select fasta formatted reference database (e.g. `UNITE for ITS reads <https://unite.ut.ee/sh_files/uchime_reference_dataset_28.06.2017.zip>`_). 
                                 | If denovo = TRUE, then reference based chimera filtering will be performed 
                                 | after denovo. 
``abundance_skew``               | the abundance skew is used to distinguish in a threeway alignment which 
                                 | sequence is the chimera and which are the parents. The assumption is that 
                                 | chimeras appear later in the PCR amplification process and are therefore 
                                 | less abundant than their parents. The default value is 2.0, which means that 
                                 | the parents should be at least 2 times more abundant than their chimera. 
                                 | Any positive value equal or greater than 1.0 can be used
``min_h``                        | minimum score (h). Increasing this value tends to reduce the number of false 
                                 | positives and to decrease sensitivity. Values ranging from 0.0 to 1.0 included 
                                 | are accepted
================================ =========================

.. _itsextractor:

____________________________________________________

|

`ITS Extractor <https://microbiology.se/software/itsx/>`_
-----------------------------------------------------------

When working with ITS amplicons, then 
extract ITS regions with `ITS Extractor <https://microbiology.se/software/itsx/>`_ (`Bengtsson-Palme et al. 2013 <https://doi.org/10.1111/2041-210X.12073>`_)

| **Fastq/fasta** formatted single-end data is supported [fastq inputs will be converted to fasta].
| **Outputs** are fasta files in ``ITSx_out`` directory.

================================ =========================
Setting                          Tooltip
================================ =========================
``organisms``                    | set of profiles to use for the search. Can be used to restrict the search to 
                                 | only a few organism groups types to save time, if one or more of the origins 
                                 | are not relevant to the dataset under study
``regions``                      | ITS regions to output (note that 'all' will output also full ITS region [ITS1-5.8S-ITS2])
``partial``                      | if larger than 0, ITSx will save additional FASTA-files for full and partial ITS sequences 
                                 | longer than the specified cutoff value. If his setting is left to 0 (zero), 
                                 | it means OFF
``e-value``                      | domain e-value cutoff a sequence must obtain in the HMMER-based step to be 
                                 | included in the output
``scores``                       | domain score cutoff that a sequence must obtain in the HMMER-based step to 
                                 | be included in the output
``domains``                      | the minimum number of domains (different HMM gene profiles) that must match 
                                 | a sequence for it to be included in the output (detected as an ITS sequence). 
                                 | Setting the value lower than two will increase the number of false positives, 
                                 | while increasing it above two will decrease ITSx detection abilities
                                 | on fragmentary data
``complement``                   | if TRUE, ITSx checks both DNA strands for matches to HMM-profiles
``only full``                    | If TRUE, the output is limited to full-length ITS1 and ITS2 regions only
``truncate``                     | removes ends of ITS sequences if they are outside of the ITS region. 
                                 | If FALSE, the whole input sequence is saved
================================ =========================

____________________________________________________

|

.. _clustering:

CLUSTERING
----------

Cluster sequences, form OTUs.

| Supported file format for the input data is **fasta**.
| **Outputs** are **OTUs.fasta** and **OTU_table.txt** files in ``clustering_out`` directory.

.. note::

 OTU table filed separator is 'tab'.

.. _clustering_vsearch:

`vsearch <https://github.com/torognes/vsearch>`_ 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=============================================== =========================
`Setting <_static/vsearch_2.18.0_manual.pdf>`_  Tooltip
=============================================== =========================
``OTU_type``                                    | centroid" = output centroid sequences; "consensus" = output 
                                                | consensus sequences
``similarity_threshold``                        | define OTUs based on the sequence similarity threshold; 0.97 = 97% 
                                                | similarity threshold
``strands``                                     | when comparing sequences with the cluster seed, check both strands 
                                                | (forward and reverse complementary) or the plus strand only
``min_OTU_size``                                | minimum read count per output OTU (e.g., if value = 2, then 
                                                | singleton OTUs will be discarded [OTUs with only one sequence])
``similarity_type``                             | pairwise sequence identity definition `--iddef <_static/vsearch_2.18.0_manual.pdf>`_
``sequence_sorting``                            | size = sort the sequences by decreasing abundance; 
                                                | "length" = sort the sequences by decreasing length (--cluster_fast); 
                                                | "no" = do not sort sequences (--cluster_smallmem --usersort)
``centroid_type``                               | "similarity" = assign representative sequence to the closest (most similar) 
                                                | centroid (distance-based greedy clustering); 
                                                | "abundance" = assign representative sequence to the most abundant centroid 
                                                | (abundance-based greedy clustering; --sizeorder), ``max_hits`` should be > 1
``max_hits``                                    | maximum number of hits to accept before stopping the search 
                                                | (should be > 1 for abundance-based selection of centroids [centroid type])
``relabel``                                     | relabel sequence identifiers (none = do not relabel)
``mask``                                        | mask regions in sequences using the "dust" method, or do not mask ("none")
``dbmask``                                      | prior the OTU table creation, mask regions in sequences using the 
                                                | "dust" method, or do not mask ("none")
``output_UC``                                   | output clustering results in tab-separated UCLAST-like format
=============================================== =========================

.. _assign_taxonomy:

____________________________________________________

|


ASSIGN TAXONOMY
---------------

Implemented tools for taxonomy annotation:

.. _assign_taxonomy_blast:

`BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_ (`Camacho et al. 2009 <https://doi.org/10.1186/1471-2105-10-421>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| BLAST search sequences againt selected database. 

.. important::

 **BLAST database needs to be an unzipped fasta file in a separate folder** (fasta will be automatically converted to BLAST database files). 
 If converted BLAST database files (.ndb, .nhr, .nin, .not, .nsq, .ntf, .nto) already exist, then just SELECT **one** of those files as BLAST database in 
 'ASSIGN TAXONOMY' panel.

| Supported file format for the input data is **fasta**.
| **Outputs** are **BLAST_1st_best_hit.txt** and **BLAST_10_best_hits.txt** files in ``taxonomy_out`` directory.

.. note::

 BLAST values filed separator is '+'. When pasting the taxonomy results to e.g. Excel, then first denote '+' as 
 as filed separator to align the columns.

================================ =========================
Setting                          Tooltip
================================ =========================
 ``database_file``               | select a database file in fasta format.
                                 | Fasta format will be automatically converted to BLAST database
``task``                         | BLAST search settings according to blastn or megablast
``strands``                      | query strand to search against database. Both = search also reverse complement
``e_value``                      | a parameter that describes the number of hits one can expect to see 
                                 | by chance when searching a database of a particular size. 
                                 | The lower the e-value the more 'significant' the match is
``word_size``                    | the size of the initial word that must be matched between the database 
                                 | and the query sequence
``reward``                       | reward for a match
``penalty``                      | penalty for a mismatch
``gap_open``                     | cost to open a gap
``gap_extend``                   | cost to extend a gap
================================ =========================

.. _databases:

A list of public databases available for taxonomy annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

==================================================================== ======== ===================================================================================================================================================================
Database                                                             Version  Description (click to download) 
==================================================================== ======== ===================================================================================================================================================================
`UNITE <https://unite.ut.ee/>`_                                      | 8.3    | `ITS region, all Eukaryotes <https://plutof.ut.ee/#/doi/10.15156/BIO/1281567>`_
`SILVA <https://www.arb-silva.de/>`_                                 | 138.1  | `16S/18S (SSU), Bacteria, Archaea and Eukarya <https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz>`_
`SILVA <https://www.arb-silva.de/projects/ssu-ref-nr/>`_ 99%         | 138.1  | `16S/18S (SSU), Bacteria, Archaea and Eukarya <https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz>`_
`MIDORI <http://www.reference-midori.info/>`_                        | 246    | `Eukaryota mitochondrial genes <http://www.reference-midori.info/download.php#>`_
`CO1 Classifier <https://github.com/terrimporter/CO1Classifier>`_    | 4      | `Metazoa COI <https://github.com/terrimporter/CO1Classifier/releases/tag/v4-ref>`_
==================================================================== ======== ===================================================================================================================================================================

____________________________________________________

.. _expert_mode:

Expert-mode (PipeCraft console)
===============================

All bioinformatic tools used by PipeCraft are available as docker images and can be used via dockers CLI.

Specify the working directory under the -v flag, 
or enter the working directory via the terminal to specify the working directory as $pwd 

docker run --interactive --tty -v $pwd:shared pipecraft/dada2:1.20 

docker run -v $pwd\:/shared -w /shared -i -t pipecraft/dada2:1.20 bash


Quit (exit from the container)

:: 

 > exit

|console|



Docker images (latest)
----------------------

====================================  =============================================================== 
Image                                 Software                                                         
====================================  ===============================================================
ewels/multiqc                         mutliqc
staphb/fastqc                         fastqc               
pipecraft/cutadapt                    cutadapt, seqkit                                      
pipecraft/dada2                       dada2                                                   
pipecraft/reorient                    fqgrep, seqkit                                                         
pipecraft/trimmomatic                 trimmomatic, seqkit                             
pipecraft/vsearch                     vsearch, seqkit           
pipecraft/itx                         ITSx, biopython, seqkit, mothur                                       
pipecraft/blast                       BLAST                           
pipecraft/deicode                     DEICODE, qiime2 
pipecraft/fastp                       fastp                            
====================================  ===============================================================

