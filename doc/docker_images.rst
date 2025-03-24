.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 50
  :alt: Alternative text
  :target: https://github.com/pipecraft2/user_guide

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red


.. _dockerimages:

===============================
Docker images |PipeCraft2_logo|
===============================

Docker images (with the bioinformatic tools) used in PipeCraft2 are stored on `Dockerhub <https://hub.docker.com/u/pipecraft>`_. 
These images can be used to launch any tool with the Docker CLI to utilize the compiled tools.


Images in use
-------------

====================================  ========================================================================== 
Image                                 Software                                                         
====================================  ==========================================================================
pipecraft/vsearch_dada2:3             vsearch v2.29.4, dada2 v1.34, seqkit v2.9.0, lulu v0.1.0, ORFfinder v0.4.3 R v4.4.1
pipecraft/vsearch_dada2_m:3           For Mac M chips, vsearch v2.29.4, dada2 v1.34, seqkit v2.9.0, lulu v0.1.0, ORFfinder v0.4.3 R v4.4.1
ewels/multiqc:1.10                    mutliqc v1.10
staphb/fastqc:0.11.9                  fastqc v0.11.9               
pipecraft/cutadapt:4.4                cutadapt v4.4, seqkit v2.3.0, python3.10.12, biopython v1.81                                       
pipecraft/reorient:1                  fqgrep v0.4.4, seqkit v2.3.0                                                       
pipecraft/trimmomatic:0.39            trimmomatic v0.39, seqkit v2.3.0                             
pipecraft/itsx:1.1.3                  ITSx v1.1.3, seqkit v2.3.0, mothur v1.46.1                                                          
pipecraft/deicode:0.2.4               DEICODE v0.2.4, qiime2-2002.2
pipecraft/fastp:0.23.2                fastp v0.23.2
pipecraft/blast:2.14                  BLAST 2.14.0+, biopython v1.81, python3.10.12, gawk v5.1.0
pipecraft/metamate:1                  metamate v0.4.0, seqkit v2.8.1, python3.10.13, biopython v1.83
pipecraft/metaworks:1.12.0            metaworks v1.12.0, seqkit v2.3.0, ORFfinder v0.4.3 R v4.1.2
pipecraft/optimotu:5                  optimotu v0.9.3, optimotu.pipeline v0.5.2 optimotu_targets v5.0.0 R v4.3 dada2 v1.30 cutadapt v5.0 vsearch v2.30.0
====================================  ==========================================================================

Other images
----------------

====================================  ================================================================================================== 
Image                                 Software                                                         
====================================  ==================================================================================================                                  
pipecraft/dada2:1.20                  dada2 v1.20, seqkit v2.3.0, lulu v0.1.0, R                                                                           
pipecraft/vsearch:2.18                vsearch v2.18, seqkit v2.3.0, GNU parallel                  
pipecraft/vsearch_dada2:1             has issues with MacOS. vsearch v2.22.1, dada2 v1.20, seqkit v2.3.0, lulu v0.1.0, R, GNU parallel
pipecraft/vsearch_dada2:2             vsearch v2.23, dada2 v1.27, seqkit v2.3.0, lulu v0.1.0
====================================  ==================================================================================================
