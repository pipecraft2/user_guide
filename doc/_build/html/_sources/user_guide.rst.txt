.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text

.. |main_interface| image:: _static/main_interface.png
  :width: 2000
  :alt: Alternative text

.. |asv_main| image:: _static/asv_main.png
  :width: 1500
  :alt: Alternative text

.. |console| image:: _static/console.png
  :width: 1500
  :alt: Alternative text

.. meta::
    :description lang=en:
        PipeCraft2 manual. User guide for PipeCraft2

|PipeCraft2_logo|
  `github <https://github.com/pipecraft2/pipecraft>`_

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red

=================
General overview
=================

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

Glossary
========

List of terms that you may encounter in this user guide.

=========================== ===================================
**working directory**       | the directory (folder) that contains the files for the analyses.
                            | The outputs will be written into this directory
**paired-end data**         | obtained by sequencing two ends of the same DNA fragment, 
                            | which results in read 1 (R1) and read 2 (R2) files per library or per sample
**single-end data**         | only one sequencing file per library or per sample. 
                            | Herein, may mean also assembled paired-end data.
**demultiplexed data**      | sequences are sorted into separate files, representing individual samples 
**multiplexed data**        | file(s) that represent a pool of sequences from different samples
**read/sequence**           | DNA sequence; herein, reads and sequences are used interchangeably 
=========================== ===================================

____________________________________________________

Docker images 
==============

.. |pulling_image| image:: _static/pulling_image.png
  :width: 280
  :alt: Alternative text


Initial PipeCraft2 installation does not contain any software for sequence data processing. 
All the processes are run through `docker <https://www.docker.com/>`_, where the PipeCraft's simply GUI mediates the 
information exchange. Therefore, whenever a process is initiated for the **first time**, 
a relevant Docker image (contains required software for the analyses step) will be pulled from `Docker Hub <https://hub.docker.com/u/pipecraft>`_.

Example: running DEMULTIPLEXING for the first time |pulling_image|

Thus working **Internet connection** is initially required. Once the Docker images are pulled, PipeCraft2 can work without an Internet connection. 

:ref:`Docker images <dockerimages>` vary in size, and the speed of the first process is extended by the docker image download time.
 
____________________________________________________

Save workflow
==============

Once the workflow settings are selected, save the workflow by pressin ``SAVE WORKFLOW`` button on the :ref:`right-ribbon <interface>`.
For saving, working directory ( ``SELECT WORKDIR`` ) does not have to be selected. 

.. note ::

  starting from version 0.1.4, PipeCraft2 will automatically save the settings into selected WORKDIR prior starting the analyses (file name = "pipecraft2_config.json")

.. important::

 When **saiving workflow** settings in **Linux**, specify the file extension as **JSON** (e.g. my_16S_ASVs_pipe.JSON).
 When trying to load the workflow, only .JSON files will be permitted as input. *Windows and Mac OS automatically extend files as JSON (so you may just save "my_16S_ASVs_pipe").*

____________________________________________________

Load workflow
==============

.. note ::

 Prior loading the workflow, make sure that the saved workflow configuration has a .JSON extension. Note also that **workflows saved in older PipeCraft2 version** might not run in newer version, but anyhow the selected options will be visible.

Press the ``LOAD WORKFLOW`` button on the :ref:`right-ribbon <interface>` and select appropriate JSON file.
The configuration will be loaded; ``SELECT WORKDIR`` and run analyses.

____________________________________________________

.. _qualitycheck:

Quality scores and basic statistics screening of the data
==========================================================

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

 * if the data is paired-end or single-end
 * and the extension of the data (fastq or fasta)

| ``paired-end data`` --> such as data from Illumina or MGI-Tech platforms (R1 and R2 files). :red:`Be sure to have **R1** and **R2** strings in the paired-end files (not simply _1 and _2; and sample names maynot contain R1/R2 strings)`
| ``single-end data`` --> such as data from PacBio, or assembled paired-end data (single file per library or per sample)

2. ``SELECT PIPELINE`` or press ``Quick Tools`` button
to select relevant :ref:`step <quicktools>` [or **load the PipeCraft settings file**]; 
edit settings if needed (**SAVE the settings for later use**) and **start
running the analyses** by pressing the ``START`` button.


.. note::

 **Step-by-step analyses**: after ``START`` is finished, then press ``SELECT WORKDIR`` to specify inputs for the next process

.. note::

 The **output files will be overwritten** if running the same 
 analysis step **multiple times in the same working directory**

3. Each process creates a separate output directory with the processed files
inside the selected working directory. 
**README** file about the process and **sequence count summary** statistics are included in the output directory.

____________________________________________________


See ready-to-run :ref:`pre-defined pipelines here <predefinedpipelines>`

____________________________________________________
