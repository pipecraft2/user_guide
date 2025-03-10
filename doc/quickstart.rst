.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 50
  :alt: Alternative text
  :target: https://github.com/pipecraft2/user_guide

.. |main_interface| image:: _static/main_interface.png
  :width: 2000
  :alt: Alternative text

.. |select_pipeline_or_quicktools| image:: _static/select_pipeline_or_quicktools.png
  :width: 1000
  :alt: Alternative text

.. |select_wd| image:: _static/select_wd.png
  :width: 1000
  :alt: Alternative text

.. meta::
    :description lang=en:
        PipeCraft2 manual. User guide for PipeCraft2

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red

.. _quickstart:

============================
QuickStart |PipeCraft2_logo|
============================

Required data formats
=====================

There are few specifc requirements for the input data for PipeCraft. 

- **Paired-end data** 
    * :red:`must contain **R1** and **R2** strings in the paired-end files`,
    * :red:`sample names maynot contain 'R1/R2' strings`
        + **OK file names**: ``my_sample_01_R1_L001.fastq`` and ``my_sample_01_R2_L001.fastq``
        + **NOT-OK** file names: ``my_R1sample_01_R1_L001.fastq`` and ``my_R1sample_01_R2_L001.fastq`` or  ``my_sample_01_1.fastq`` and ``my_sample_01_2.fastq``
    * :red:`Paired-end reads should have common read identifiers (e.g. _R1/_R2) and file extensions (e.g. fastq/fq) for the reads that are processed together`.


- **index/barcodes file** for demultiplexing:
    * :ref:`see formatting requirements here <indexes>` 
  
- **Please avoid spaces and non-ASCII symbols** in sample names.

- **Use at least 2 samples per sequencing run** for the pre-defined pipelines.
  
- specific directory structure for NextITS pipeline 
    * :ref:`see NextITS page here <nextits_pipeline>` 

- specific directory structure for OptimOTU pipeline 
    * :ref:`see OptimOTU page here <optimotu_pipeline>` 

- specific directory structure when aiming to combine multiple sequencing runs with vsearch, unoise, or DADA2 pipelines
    * :ref:`see here <multi_run_dir>` 
  



____________________________________________________

How to START
============

1. To ``START`` any analyses, you must specify the working directory (WORKDIR) by pressing the ``SELECT WORKDIR`` button. E.g., if working with **fastq** files,
then be sure that the working directory contains **only relevant fastq files** because the selected process will be 
applied to all fastq files in the working directory!

.. note::

 When using Windows OS, the selection window might not display the files while browsing through the directories. 

After selecting a working directory, PipeCraft needs you to specify if 

 * if the data is paired-end or single-end
 * and the extension of the data (fastq or fasta)

| ``paired-end data`` --> such as data from Illumina or MGI-Tech platforms (R1 and R2 files). :red:`Be sure to have **R1** and **R2** strings in the paired-end files (not simply _1 and _2; and sample names maynot contain R1/R2 strings)`
| ``single-end data`` --> such as data from PacBio, or assembled paired-end data (single file per library or per sample)

|select_wd|

2. ``SELECT PIPELINE`` or press ``Quick Tools`` button
to select relevant :ref:`step <quicktools>`; 
edit settings if needed and **start
running the analyses** by pressing the ``START`` button.

|select_pipeline_or_quicktools|


.. note::

 **When running 'step-by-step analyses with Quick Tools'**: when one workflow is finished, then press ``SELECT WORKDIR`` to specify inputs for the next process to ensure the correct workflow piping.  


.. note:: 

  Note that a ``pipecraft2_last_run_configuration.json`` file will be generated into the working directory upon starting a workflow.
  As the file name says, it is the workflow configuration file for your last PipeCraft run in this **working directory**. 

  This ``JSON`` file can be loaded into PipeCraft2 to **automatically configure your next runs exactly the same way**.

  
.. warning::

 The **outputs will be overwritten** if running the same 
 analysis step **multiple times in the same working directory**.
 If needed, edit the default output directory name to prevent that.


Each process creates a separate output directory with the processed files. 
The **README** file in the output directory states some of the details about the finished process.

____________________________________________________


 :ref:`Ready-to-run pre-defined pipelines here <predefinedpipelines>`

 :ref:`QuickTools page here <quicktools>`

____________________________________________________


Save workflow
==============

Once the workflow settings are selected, save the workflow by pressin ``SAVE WORKFLOW`` button on the :ref:`right-ribbon <interface>`.

.. note ::

  starting from version 0.1.4, PipeCraft2 will automatically save the settings into selected WORKDIR prior starting the analyses (file name = "**pipecraft2_last_run_configuration.json**")

.. important::

 When **saving workflow** settings in **Linux**, specify the file extension as **json** (e.g. my_16S_ASVs_pipe.json).
 When trying to load the workflow, only .JSON files will be permitted as input. *Windows and Mac OS automatically extend files as json (so you may just save "my_16S_ASVs_pipe").*

____________________________________________________

Load workflow
==============

| Press the ``LOAD WORKFLOW`` button on the :ref:`right-ribbon <interface>` and select appropriate JSON file.
| The configuration will be loaded. 
| **Then you need to** ``SELECT WORKDIR`` and after that may run PipeCraft.

.. note ::

 Prior loading the workflow, make sure that the saved workflow configuration has a .json extension. Note also that **workflows saved in older PipeCraft2 version** might not run in newer version, but anyhow the selected options will be visible.



