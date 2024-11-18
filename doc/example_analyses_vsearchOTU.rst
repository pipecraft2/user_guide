.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 50
  :alt: Alternative text
  :target: https://github.com/pipecraft2/user_guide

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red

.. raw:: html

    <style> .green {color:#00f03d; font-weight:bold; font-size:16px} </style>

.. role:: green

.. |fastqc_per_base_sequence_quality_plot| image:: _static/fastqc_per_base_sequence_quality_plot.png
  :width: 850
  :alt: Alternative text

.. |workflow_finished| image:: _static/workflow_finished.png
  :width: 300
  :alt: Alternative text

.. |stop_workflow| image:: _static/stop_workflow.png
  :width: 200
  :alt: Alternative text

.. meta::
    :description lang=en:
        PipeCraft manual. tutorial

|

vsearch OTUs pipeline
---------------------

This example data analyses follows vsearch OTUs workflow as implemented in PipeCraft2's pre-compiled pipelines panel. 


Starting with demultiplexed paired-end Illumina fastq files.




.. note::

 Built-in panel for OTU workflow with (mostly) vsearch.

Here, we perform example analyses of paired-end data using `mothur MiSeq SOP example data set <https://mothur.org/wiki/miseq_sop/>`_.
`Download example data set here <https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip>`_ (35.1 Mb) and unzip it. 
This data set represents demultiplexed set (per-sample fastq files) of 16S rRNA gene V4 amplicon sequences where sample indexes and primers have already been removed. 

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


