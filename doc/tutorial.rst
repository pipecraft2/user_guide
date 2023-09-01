.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
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

=====================
Example data analyses 
=====================

Some example analyses pipelines.

.. note::

  When samples of interest are distributed between different sequencing libraries, then first demultiplex (if needed)
  libraries separately and place samples of interest into separate working directory.

.. contents:: Contents
   :depth: 4

____________________________________________________


Inspect quality profiles
========================

Examine the quality profiles and basic statistics of the your data set using **QualityCheck module**. :ref:`See here <qualitycheck>`.

____________________________________________________

|

Paired-end Illumina (or MGI-Tech) data
======================================
Example analyses of paired-end data. Starting with raw paired-end fastq files, finishing with ASV/OTU table and taxonomy table. 



.. _demultiplexed_PE_ASV:

Demultiplexed paired-end data; ASV workflow with DADA2
------------------------------------------------------

.. note::

 This tutorial follows `DADA2 Pipeline Tutorial <https://benjjneb.github.io/dada2/tutorial.html>`_.

Here, we perform example analyses of paired-end data using `mothur MiSeq SOP example data set <https://mothur.org/wiki/miseq_sop/>`_.
`Download example data set here <https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip>`_ (35.1 Mb) and unzip it. 
This data set represents demultiplexed set (per-sample fastq files) of 16S rRNA gene V4 amplicon sequences where sample indexes and primers have already been removed. 

* If working with multiplexed data, :ref:`see here <multiplexed_PE>`.
* If you need to trim the primers/adapters, :ref:`see here <remove_primers>`.

.. warning::
 
 Be sure that all sequences have **same orientation** (5'-3' or 3'-5') in your input data set(s)! If sequences are in **mixed orientation** 
 (i.e. some sequences are recorded as 5'-3' and some as 3'-5'; as usually in the raw data), 
 then exactly the same ASV may be reported twice, where one is just the reverse complementary ASV: 1) ASV with sequence orientation of 5'-3'; and 2) ASV with sequence orientation of 3'-5'. **Reorient sequences** based on primer sequene using ``REORIENT`` panel; :ref:`see here <reorient>`.
 
.. important::

  When working with your own data, then please check that the paired-end data file names contain "R1" and "R2" strings 
  (to correctly identify the paired-end reads by PipeCraft). 

  | *Example:*
  | *F3D0_S188_L001_R1_001.fastq*
  | *F3D0_S188_L001_R2_001.fastq*


:: 

 1. Select working directory by pressing the 'SELECT WORKDIR' button.

| Secify
| ``sequencing data format`` as **demultiplexed**;
| ``sequence files extension`` as **\*.fastq**;  
| ``sequencing read types`` as **paired-end**.

:: 

 2. Select 'ASVs workflow' panel (right-ribbon) and check that docker is running (green icon);

.. |ASVs_tut1| image:: _static/ASVs_tut1.png
  :width: 1000
  :alt: Alternative text

* Here, working with demultiplexed data, where primers have already been removed; so **do not tick** ``DEMULTIPLEX``, ``REORIENT``, ``CUT PRIMERS`` 
  (:ref:`see here <multiplexed_PE>` to analyse multiplexed data, and :ref:`here <remove_primers>` if you need to cut primers/adapters).
  
|ASVs_tut1|

| 
|

.. _rest_of_PE_ASV_workflow:

:: 

 3. 'QUALITY FILTERING'
   
.. |DADA2_quality_filt_expand_dadaTut| image:: _static/DADA2_quality_filt_expand_dadaTut.png
  :width: 550
  :alt: Alternative text

.. |fastqc_per_base_sequence_quality_plot| image:: _static/fastqc_per_base_sequence_quality_plot.png
  :width: 850
  :alt: Alternative text

Before adjusting quality filtering settings, let's have a look on the **quality profile** of our example data set. Below quality profile plot was generated using ``QualityCheck`` panel (:ref:`see here <qualitycheck>`).

|fastqc_per_base_sequence_quality_plot|

In this case, all **R1 files are represented with green lines**, indicating good average quality per file. However, all **R2 files are either yellow or red**, indicating a drop in quality scores. 
Lower qualities of R2 reads are characteristic for Illumina sequencing data, and is not too alarming. DADA2 algoritm is robust to lower quality sequences, but removing the low quality read parts
will improve the DADA2 sensitivity to rare sequence variants.


* **Click on** ``QUALITY FILTERING`` **to expand the panel**
* specify identifier strings for ``read R1`` and ``read R2``. Here, fastq file names = F3D0_S188_L001_R1_001.fastq, F3D0_S188_L001_R2_001.fastq etc...; **_R1** and **_R2** are common identifiers for all files.
* specify ``samp ID`` (sample identifier). Here **_** (underscore), which denotes that sample name is a string before the first **_** in the fastq file name.
* trim reads to specified length to remove low quality ends. Set ``truncLen`` to 240 for trimming R1 reads and ``truncLen R2`` to 160 to trim R2 reads. Latter positions represent the approximate positions where sequence quality drps notably
  (quality profile figure above). Be sure to consider the amplicon length before applying ``truncLen`` options, so that R1 and R2 reads would still overlap for the ``MERGE PAIRS`` process. 
* other settings as default.

*(click on the image for enlargement)*
|DADA2_quality_filt_expand_dadaTut|

| *This step performs quality filtering.* 
| *Quality filtering settings* :ref:`here <dada2_qual_filt>`
| 
| **Output** directory = ``qualFiltered_out``:
| \*_filt.fastq          = quality filtered sequences per sample in FASTQ format
| seq_count_summary.txt = summary of sequence counts per sample
| FASTA/\*_filt.fasta    = quality filtered sequences per sample in FASTA format

:: 

 4. Here, we use default 'DENOISE' and 'MERGE PAIRS' settings 


| *This step performs denoising and merging of paired-end sequences.* 
| *Denoise settings* : :ref:`here <dada2_denoise>`, *merge pairs settings* :ref:`here <dada2_merge_pairs>`)
|
| **Output** directory = ``denoised_assembled.dada2``. 
| \*.merged_ASVs.fasta   = denoised and assembled ASVs per sample. 'Size' denotes the abundance of the ASV sequence
| Error_rates_R1.pdf    = plots for estimated R1 error rates
| Error_rates_R2.pdf    = plots for estimated R2 error rates
| seq_count_summary.txt = summary of sequence and ASV counts per sample

:: 

 5. Default settings for 'CHIMERA FILTERING'

(method = consensus)

| *This step performs chimera filtering on denoised and merged reads.*
| *ASV table is generated during this step* 
| *Chimera filtering settings* :ref:`here <dada2_chimeras>`
|
| **Output** directories -> 
| ``chimeraFiltered_out.dada2``: 
| \*.chimFilt_ASVs.fasta = chimera filtered ASVs per sample. 'Size' denotes the abundance of the ASV sequence.  
| seq_count_summary.txt = summary of sequence and ASV counts per sample
| \*.chimeras.fasta      = ASVs per sample that were flagged as chimeras (and thus discarded)

| ``ASVs_out.dada2``: 
| ASVs_table.txt = ASV distribution table per sample (tab delimited file)
| ASVs.fasta     = FASTA formatted representative ASV sequences (this file is used for taxonomy assignment)

|

:: 

 6. 'ASSIGN TAXONOMY'

* Click on 'ASSIGN TAXONOMY' to expand the panel  
* press ``DOWNLOAD DATABASES`` which direct you to the DADA2-formatted reference databases `web page <https://benjjneb.github.io/dada2/training.html>`_.
* download SILVA (silva_nr99_v138.1_wSpecies_train_set.fa.gz) database for assigning taxonomy to our 16S ASVs. `Download link here <https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1>`_
* specify the location of your downloaded DADA2 database by pressing ``SELECT FASTA``
* since primers were already removed from this data set, we could not :ref:`reorient all sequences to uniform orientation as based on primers <reorient>`. Therefore, **swithc ON** ``tryRC`` 
  to include reverse-complement search. 
  
.. |DADA2_assign_taxRC| image:: _static/DADA2_assign_taxRC.png
  :width: 550
  :alt: Alternative text

|DADA2_assign_taxRC|

| *This step assigns taxonomy to ASVs using DADA2* `assignTaxonomy <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function.
| *Assign taxonomy settings* :ref:`here <dada2_taxonomy>`
|
| **Output** directory = ``taxonomy_out.dada2``:
| taxonomy.txt = classifier results with bootstrap values


:: 

 6.1. Save the workflow by pressing ``SAVE WORKFLOW`` button on the right-ribbon.

::

 7.  Press** 'START' **to start the analyses.

.. |workflow_finished| image:: _static/workflow_finished.png
  :width: 300
  :alt: Alternative text

.. note ::

  When running the panel for the first time, a docker image will be pulled first to start the process.

:: 

 When done, 'workflow finished' window will be displayed.

|workflow_finished|

.. |stop_workflow| image:: _static/stop_workflow.png
  :width: 200
  :alt: Alternative text

.. note::
 
 Press ``STOP WORKFLOW`` to stop. 
   |stop_workflow|

|

->

Examine the outputs
~~~~~~~~~~~~~~~~~~~

Several process-specific output folders are generated:

| ``qualFiltered_out`` -> quality filtered paired-end **fastq** files per sample
| ``denoised_assembled.dada2`` -> denoised and assembled **fasta** files per sample (and error rate plots)
| ``chimeraFiltered_out.dada2`` --> chimera filtered **fasta** files per sample
| ``ASVs_out.dada2`` --> **ASVs table** (ASVs_table.txt), and ASV sequences (ASVs.fasta) file
| ``taxonomy_out.dada2``--> ASVs **taxonomy table** (taxonomy.txt)

.. _seq_count_summary:

Each folder (except ASVs_out.dada2 and taxonomy_out.dada2) contain 
**summary of the sequence counts** (seq_count_summary.txt). 
Examine those to track the read counts throughout the pipeline. 

 For example, merging the seq_count_summary.txt file in ``qualFiltered_out`` with the seq_count_summary.txt file from ``chimeraFiltered_out.dada2`` forms a table for examining sequence counts throughout the 
 pipeline and number of ASVs per sample. 

======= ===== ============ ====== ================ ==========
sample  input qualFiltered merged chimeraFiltered  no.of ASVs
======= ===== ============ ====== ================ ==========
F3D0    7793  7113         6540   6528             106
F3D141  5958  5463         4986   4863             74
F3D142  3183  2914         2595   2521             48
F3D143  3178  2941         2552   2518             56
F3D144  4827  4312         3627   3488             47
F3D145  7377  6741         6079   5820             72
F3D146  5021  4560         3968   3879             84
F3D147  17070 15637        14231  13006            103
F3D148  12405 11413        10529  9935             97
F3D149  13083 12017        11154  10653            112
F3D150  5509  5032         4349   4240             78
F3D1    5869  5299         5028   5017             100
F3D2    19620 18075        17431  16835            134
F3D3    6758  6250         5853   5491             68
F3D5    4448  4052         3716   3716             86
F3D6    7989  7369         6865   6679             90
F3D7    5129  4765         4428   4217             61
F3D8    5294  4871         4576   4547             99
F3D9    7070  6504         6092   6015             106
Mock    4779  4314         4269   4269             20
======= ===== ============ ====== ================ ==========

|

``ASVs_out.dada2`` directory contains **ASVs table** (ASVs_table.txt), where the **1st column** represents ASV identifiers, 
**2nd column** representative sequences of ASVs,
and all following columns represent samples (number of sequences per ASV in a sample). This is tab delimited text file.

*ASVs_table.txt; first 4 samples*

===== ==============  ===== ======  ======  ======
ASV   Sequence        F3D0  F3D141  F3D142  F3D143
===== ==============  ===== ======  ======  ======
ASV_1 TACGGAGGATG...  579   444     289     228
ASV_2 TACGGAGGATG...  345   362     304     176
ASV_3 TACGGAGGATG...  449   345     158     204
ASV_4 TACGGAGGATG...  430   502     164     231
ASV_5 TACGGAGGATC...  154   189     180     130
ASV_6 TACGGAGGATG...  470   331     181     244
ASV_7 TACGGAGGATG...  282   243     163     152
ASV_8 TACGGAGGATT...  184   321     89      83
ASV_9 TACGGAGGATG...  45    167     89      109
===== ==============  ===== ======  ======  ======

The **ASV sequences** are representad also in the fasta file (ASVs.fasta) in ``ASVs_out.dada2`` directory. 

Result from the taxonomy annotation process - **taxonomy table** (taxonomy.txt) - is located at the ``taxonomy_out.dada2`` directory. 
"NA" denotes that the ASV was not assigned to corresponding taxonomic unit.  
Last columns with integers (for 'Kingdom' to 'Species') represent bootstrap values for the correspoinding taxonomic unit. 

*taxonomy.txt; first 10 ASVs*

=======  ========== ======== ============ =========== ===============  ===============  ============================== ========== ======= ====== ===== ===== ====== ===== =======
ASV      Sequence   Kingdom   Phylum      Class       Order            Family           Genus                          Species    Kingdom Phylum Class Order Family Genus Species
=======  ========== ======== ============ =========== ===============  ===============  ============================== ========== ======= ====== ===== ===== ====== ===== =======
ASV_1    TACGGAG... Bacteria Bacteroidota Bacteroidia Bacteroidales    Muribaculaceae   NA                             NA         100     100    100   100   100    100   100
ASV_2    TACGGAG... Bacteria Bacteroidota Bacteroidia Bacteroidales    Muribaculaceae   NA                             NA         100     100    100   100   100    100   100
ASV_3    TACGGAG... Bacteria Bacteroidota Bacteroidia Bacteroidales    Muribaculaceae   NA                             NA         100     100    100   100   100    100   100
ASV_4    TACGGAG... Bacteria Bacteroidota Bacteroidia Bacteroidales    Rikenellaceae    Alistipes                      NA         100     100    100   100   100    100   100
ASV_5    TACGGAG... Bacteria Bacteroidota Bacteroidia Bacteroidales    Muribaculaceae   NA                             NA         100     100    100   100   100    100   100
ASV_6    TACGGAG... Bacteria Bacteroidota Bacteroidia Bacteroidales    Muribaculaceae   NA                             NA         100     100    100   100   100    95    95
ASV_7    TACGTAG... Bacteria Firmicutes   Clostridia  Lachnospirales   Lachnospiraceae  Lachnospiraceae NK4A136 group  NA         100     100    100   100   100    100   99
ASV_8    TACGGAG... Bacteria Bacteroidota Bacteroidia Bacteroidales    Muribaculaceae   NA                             NA         100     100    100   100   100    100   100
ASV_9    TACGGAG... Bacteria Bacteroidota Bacteroidia Bacteroidales    Bacteroidaceae   Bacteroides                    caecimuris 100     100    100   100   100    100   77
ASV_10   TACGGAG... Bacteria Bacteroidota Bacteroidia Bacteroidales    Muribaculaceae   NA                             NA         100     100    100   100   100    99    99
=======  ========== ======== ============ =========== ===============  ===============  ============================== ========== ======= ====== ===== ===== ====== ===== =======


____________________________________________________

|

.. _demultiplexed_PE_OTU:

Demultiplexed paired-end data; OTU workflow
-------------------------------------------

.. note::

 Built-in panel for OTU workflow with (mostly) vsearch.

Here, we perform example analyses of paired-end data using `mothur MiSeq SOP example data set <https://mothur.org/wiki/miseq_sop/>`_.
`Download example data set here <https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip>`_ (35.1 Mb) and unzip it. 
This data set represents demultiplexed set (per-sample fastq files) of 16S rRNA gene V4 amplicon sequences where sample indexes and primers have already been removed. 

* If working with multiplexed data, :ref:`see here <multiplexed_PE>`.
* If you need to trim the primers/adapters, :ref:`see here <remove_primers>`.

.. note::

 When working with your own data, then consider **reorienting** reads; :ref:`see here <reorient>`. Although, in the OTU formation step (clustering), both sequence strands 
 will be compared to generate OTUs, the time for BLAST (taxonomy annotation step) can be reduced when there is no need to search reverse complementary matches. 

.. important::

  When working with your own data, then please check that the paired-end data file names contain "R1" and "R2" strings 
  (to correctly identify the paired-end reads by PipeCraft). 

  | *Example:*
  | *F3D0_S188_L001_R1_001.fastq*
  | *F3D0_S188_L001_R2_001.fastq*


:: 

 1. Select working directory by pressing the 'SELECT WORKDIR' button.

| Secify
| ``sequencing data format`` as **demultiplexed**;
| ``sequence files extension`` as **\*.fastq**;  
| ``sequencing read types`` as **paired-end**.

:: 

 2. Select 'OTU workflow' panel (right-ribbon) and check that docker is running (green icon);

.. |OTUs_workflow| image:: _static/OTUs_workflow.png
  :width: 1000
  :alt: Alternative text

* Here, working with demultiplexed data, where primers have already been removed; so **do not tick** ``DEMULTIPLEX``, ``REORIENT``, ``CUT PRIMERS`` 
  (:ref:`see here <multiplexed_PE>` to analyse multiplexed data, and :ref:`here <remove_primers>` if you need to cut primers/adapters).
  
|OTUs_workflow|

| 
|

Before proceeding, let's have a look on the **quality profile** of our example data set. Below quality profile plot was generated using ``QualityCheck`` panel (:ref:`see here <qualitycheck>`).

|fastqc_per_base_sequence_quality_plot|

In this case, all **R1 files are represented with green lines**, indicating good average quality per file. However, all **R2 files are either yellow or red**, indicating a drop in quality scores. 
Lower qualities of R2 reads are characteristic for Illumina sequencing data, and is not too alarming. Nevertheless, we need to quality filter the data set. 

.. _rest_of_PE_OTU_workflow:

:: 

 3. 'MERGE PAIRS' 

* Here, we use default settings.
  
.. note:: 

 If ``include_only_R1`` option = TRUE, 
 then unassembled R1 reads will be included to the set of assembled reads per sample. 
 This may be useful when working with e.g. ITS2 sequences, because the ITS2 region in some 
 taxa is too long for paired-end assembly using current short-read sequencing technology. 
 Therefore longer ITS2 amplicon sequences are discarded completely after the assembly process. 
 Thus, including also unassembled R1 reads (``include_only_R1`` = TRUE), partial ITS2 sequences for 
 these taxa will be represented in the final output. But when using :ref:`ITSx <itsextractor>`  
 , keep ``only_full`` = FALSE and include ``partial`` = 50.
 | 
 If include only R1 option = TRUE, then other specified options (lenght, max error rate etc.) have not been 
 applied to R1 reads in the 'assembled' file. Thus, additional quality filtering (if this was done before assembling) 
 should be run on the 'assembled' data. But in this built-in OTU workflow, the quality filtering step is anyway performed after merge pairs step. 


| *This step performs merging of paired-end sequences using vsearch --fastq_mergepairs.* 
| *Merge pairs settings* :ref:`here <merge_pairs>`)
|
| **Output** directory = ``assembled_out``. 

:: 

 4. 'QUALITY FILTERING'

.. |vsearch_qfilt| image:: _static/vsearch_qfilt.png
  :width: 600
  :alt: Alternative text   

* **Click on** ``QUALITY FILTERING`` **to expand the panel**
* specify ``maxee`` (maximum number of expected errors per sequence), here we use 1 (`see here what is maxee <https://drive5.com/usearch/manual/exp_errs.html>`_).
* specify ``maxNs`` (maximum number of Ns in the sequences). Here, we will discard any sequence that contains N (ambiguously recorded nucleotide) by setting the value to 0.
* other settings as default.

|vsearch_qfilt|

| *This step performs quality filtering using vsearch.* 
| *vsearch quality filtering settings* :ref:`here <qfilt_vsearch>`
| 
| **Output** directory = ``qualFiltered_out``

|

:: 

 5. 'CHIMERA FILTERING'

* **Click on** ``CHIMERA FILTERING`` **to expand the panel**
* specify ``pre cluster`` threshold as 0.97 (that is 97%; when planning to use 97% sequence similarity threshold also for clustering reads into OTUs).
* here, we perform only ``denovo`` chimera filtering 
* other settings as default.

.. note::

 Tick ``reference based`` if there is appropriate database for reference based chimera filtering 
 (such as e.g. `UNITE for ITS reads <https://unite.ut.ee/sh_files/uchime_reference_dataset_28.06.2017.zip>`_).

.. |vsearch_chimeraFilt| image:: _static/vsearch_chimeraFilt.png
  :width: 600
  :alt: Alternative text   

|vsearch_chimeraFilt|

| *This step performs chimera filtering using vsearch* 
| *Chimera filtering settings* :ref:`here <chimFilt>`
|
| **Output** directory = ``chimeraFiltered_out``

|

:: 

 6. Consideration when working with ITS data

Identify and extract the ITS regions using ITSx; :ref:`see here <itsextractor>`

.. note::

  because ITSx outputs multiple directories for different ITS sub-regions
  ``CLUSTERING`` and ``ASSIGN TAXONOMY`` will be disabled after 'ITS EXTRACTOR'.
  Select appropriate ITSx output folder for CLUSTERING after the process is finished 
  ['ADD STEP' -> ``CLUSTERING`` (vsearch)].

| *This step extracts ITS reads using ITSx* 
| *ITSx settings* :ref:`here <itsextractor>`
|
| **Output** directory = ``ITSx_out`` 
| 

::

 7. 'CLUSTERING' 

* Here, we use default settings by clustering the reads using 97% similarity threshold

| *This step performs clustering using vsearch.* 
| *vsearch clustering settings* :ref:`here <clustering>`
|
| **Output** directory = ``clustering_out`` 

| 

::

  8. 'ASSIGN TAXONOMY'

.. |assign_taxonomy_blast| image:: _static/assign_taxonomy_blast.png
  :width: 600
  :alt: Alternative text   

* Tick ``ASSIGN TAXONOMY`` to perform taxonomy assignment with BLAST
* download SILVA 99% database :ref:`here <databases>` (SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz)
* **unzip** the downloaded database and place this into separete folder (to automatically make blast database from that fasta file)
* specify the location of your downloaded SILVA database by pressing ``SELECT FILE`` under 'database file' option
* since primers were already removed from this data set, we could not :ref:`reorient all sequences to uniform orientation as based on primers <reorient>`. Therefore, **keep ON** the ``strands`` = both to include reverse-complement search. 


|assign_taxonomy_blast|


| *This step assigns taxonomy to OTUs using BLAST* 
| *Assign taxonomy settings* :ref:`here <assign_taxonomy>`
|
| **Output** directory = ``taxonomy_out`` 


:: 

 8.1. Save the workflow by pressing ``SAVE WORKFLOW`` button on the right-ribbon.

::

 1.  Press** 'START' **to start the analyses.

.. note ::

  When running the panel for the first time, a docker image will be pulled first to start the process.


:: 

 When done, 'workflow finished' window will be displayed.

|workflow_finished|

.. note::
 
 Press ``STOP WORKFLOW`` to stop. 
   |stop_workflow|

|

->

Examine the outputs
~~~~~~~~~~~~~~~~~~~

Several process-specific output folders are generated:

| ``assembled_out`` --> assembled **fastq** files per sample
| ``qualFiltered_out`` --> quality filtered **fastq** files per sample
| ``chimeraFiltered_out`` --> chimera filtered **fasta** files per sample
| ``clustering_out`` --> **OTU table** (OTU_table.txt), and OTU sequences (OTUs.fasta) file
| ``taxonomy_out``--> BLAST hits for the OTUs (BLAST_1st_best_hit.txt and BLAST_10_best_hits.txt)


Each folder (except clustering_out and taxonomy_out) contain 
**summary of the sequence counts** (seq_count_summary.txt). 
Examine those to track the read counts throughout the pipeline (:ref:`example here <seq_count_summary>`)


``clustering_out`` directory contains **OTUs table** (OTUs_table.txt), where the **1st column** represents OTU identifiers, 
and all following columns represent samples (number of sequences per OTU in a sample).
The **OTU sequences** are representad in the fasta file (OTUs.fasta) in ``clustering_out`` directory. 

*OTUs_table.txt; first 4 samples*

========================================  ============== ============== ============== ==============
OTU_id                                    F3D0_S188_L001 F3D1_S189_L001 F3D2_S190_L001 F3D3_S191_L001
========================================  ============== ============== ============== ==============
00fc1569196587dde0462c7d806cc05774f61bfa  106            271            584            20
02d84ed0175c2c79e8379a99cffb6dbc7f6a6bd9  81             44             88             14
0407ee3bd15ca7206a75d02bb41732516adaaa88  3              4              3              0
042e5f0b5e38dff09f7ad58b6849fb17ec5503b9  20             83             131            4
07411b848fcda497fd29944d351b8a2ec7dc2bd4  1              0              2              0
07e7806a732c67ef090b6b279b74a87fefad9e8e  18             22             83             7
0836d270877aed22cd247f7e703b9247fb339127  1              1              0              0
0aa6e7da5819c11973f186cb35b3f4f58275fb04  1              4              5              0
0c1c219a4756bb729e5f0ceb7d82d932bbfa0c5e  18             17             40             7
========================================  ============== ============== ============== ==============


Results from the taxonomy annotation process (BLAST) are located at the ``taxonomy_out`` directory (BLAST_1st_best_hit.txt and BLAST_10_best_hits.txt).
**Blast values are separated by** ``+`` and ``tab`` [be sure to specify the delimiter when aligning columns in e.g. LibreOffice or Excel]. 
"NO_BLAST_HIT" denotes that the OTU sequence did not get any match againt the selected database. 


============= =================================================
blast values 
============= =================================================
score          blast score
e-value        blast e-value
query len      query (i.e. OTU/ASV) sequence length
query start    start position of match in the query seq
query end      end position of match in the query seq
target len     target seq length in the database
target start   start position of match in the target seq
target end     end position of match in the target seq
align len      alignment length of query and target
identities     number of identical matches
gaps           number of gaps in the alignment
coverage      | query coverage percentage against the target sequence 
              | (100 percent is full-length match; 
              | low coverage may indicate presence of **chimeric** sequence/OTU/ASV)
id             identity percentage against the target sequence
============= =================================================


____________________________________________________

|

.. _multiplexed_PE:

Multiplexed library
-------------------

Working with paired-end raw multiplexed data. 

:: 

 1. Select working directory by pressing the 'SELECT WORKDIR' button.

| Secify
| ``sequencing data format`` as **multiplexed**;
| ``sequence files extension`` as **may be fastq or fasta** formatted files;  
| ``sequencing read types`` as **paired-end**.

:: 

 2. 'DEMULTIPLEX'

2.1 Press ``ADD STEP`` -> ``DEMULTIPLEX``

.. |add_setp_demux| image:: _static/add_setp_demux.png
  :width: 300
  :alt: Alternative text

|add_setp_demux|

or 

2.2. Select ``ASVs workflow`` or ``OTUs workflow`` panel

* tick ``DEMULTIPLEX``, ``REORIENT`` and ``CUT PRIMERS``;
* check that the docker is running (green icon [red = not running])

.. |ASVs_demux_tut1| image:: _static/ASVs_demux_tut1.png
  :width: 550
  :alt: Alternative text

*(click on the image for enlargement)*
|ASVs_demux_tut1|

|

::

 3. Click on 'DEMULTIPLEX' to expand the panel

.. |demux_expand| image:: _static/demux_expand.png
  :width: 550
  :alt: Alternative text

* select your FASTA formatted **index_file.fasta** (:ref:`general index file guide here <indexes>`)
* adjust ``overlap`` setting to fully match the length (in base pairs) of the indexes in the index_file.fasta. 
  
*(click on the image for enlargement)*
|demux_expand|

| This step distributes sequences to samples according to the information in the index_file.fasta. See :ref:`specifics here <demux_settings>`
| 
| **Output** directory = ``demultiplex_out``:
| * fastq or fasta files per sample (as specified in the :ref:`index file <indexes>`)
| * unknown.fastq/fasta files contain sequences where specified index combinations were not found. 

|

::
  
  1.  'REORIENT'

.. |reorient_expand| image:: _static/reorient_expand.png
  :width: 550
  :alt: Alternative text

* specify allowed ``mismatches`` during the primer search; >2 not recommended.
* specify ``forward primer``: 5'-GTGYCAGCMGCCGCGGTAA-3' (example)
* specify ``reverse primer``: 3'-GGCCGYCAATTYMTTTRAGTTT-5' (example)

*(click on the image for enlargement)*
|reorient_expand|

| *This step reorients sequences to 5'-3' as based on specified forward and reverse primers. See* :ref:`specifics here <reorient>`
| 
| **Output** directory = ``reorient_out``

|

::

 5. Click on 'CUT PRIMERS' to expand the panel

.. |cut_primers_expand| image:: _static/cut_primers_expand.png
  :width: 550
  :alt: Alternative text

* specify ``forward primer``: 5'-GTGYCAGCMGCCGCGGTAA-3' (example)
* specify ``reverse primer``: 3'-GGCCGYCAATTYMTTTRAGTTT-5' (example)
* specify allowed ``mismatches`` during the primer search; >2 not recommended
* for paired-end reads keep ``seqs to keep`` and ``pair filter`` as default (**keep_all** and **both**, respectively)


*(click on the image for enlargement)*
|cut_primers_expand|

| *This step clipps specified primer sequences from the reads (if primers are found). See* :ref:`specifics here <remove_primers>`.
| *Discards the reads where primer sequences are not detected.*
|
| **Output** directory = ``primersCut_out``

| 

**6.** Follow the rest of the :ref:`ASV workflow <rest_of_PE_ASV_workflow>` or :ref:`OTU workflow <rest_of_PE_OTU_workflow>`

____________________________________________________



Single-end (PacBio or assembled paired-end) data
================================================

coming soon ...

.. hide:

  .. _multiplexed_SE_OTU:

  Multiplexed single-end data; OTU workflow
  -----------------------------------------

  .. note::
    ASV workflow currently not available for single-end data

  reorient not needed when cutting primers


  Examine the outputs
  ~~~~~~~~~~~~~~~~~~~~

  ____________________________________________________


  .. _demultiplexed_SE_OTU:

  Demultiplexed single-end data; OTU workflow
  -------------------------------------------

  .. note::
    ASV workflow currently not available for single-end data

  reorient not needed when cutting primers


  Examine the outputs
  ~~~~~~~~~~~~~~~~~~~~

