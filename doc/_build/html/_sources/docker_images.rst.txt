.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text


|PipeCraft2_logo|
  `github <https://github.com/pipecraft2/pipecraft>`_
 

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red


.. _dockerimages:

===============
Docker images
===============

Bioinformatic tools used by PipeCraft2 are stored on `Dockerhub <https://hub.docker.com/u/pipecraft>`_ as Docker images. 
These images can be used to launch any tool with the Docker CLI to utilize the compiled tools.


Images in use
---------------

====================================  ========================================================================== 
Image                                 Software                                                         
====================================  ==========================================================================
pipecraft/vsearch_dada2:1             vsearch v2.22.1, dada2 v1.20, seqkit v2.3.0, lulu v0.1.0, R, GNU parallel
ewels/multiqc:latest                  mutliqc v1.12
staphb/fastqc:0.11.9                  fastqc v0.11.9               
pipecraft/cutadapt:3.5                cutadapt v3.5, seqkit v2.3.0, python3, biopython                                        
pipecraft/reorient:1                  fqgrep v0.4.4, seqkit v2.3.0                                                       
pipecraft/trimmomatic:0.39            trimmomatic 0.39, seqkit v2.3.0                             
pipecraft/itsx:1.1.3                  ITSx v1.1.3, seqkit v2.3.0, mothur v1.46.1                                                          
pipecraft/deicode:0.2.4               DEICODE v0.2.4, qiime2-2002.2
pipecraft/fastp:0.23.2                fastp v0.23.2
pipecraft/blast:2.12                  BLAST 2.12.0+, biopython, python3, gawk                             
====================================  ==========================================================================

Other images
----------------

====================================  ========================================================================== 
Image                                 Software                                                         
====================================  ==========================================================================                                  
pipecraft/dada2:1.20                  dada2 v1.20, seqkit v2.3.0, lulu v0.1.0, R                                                                           
pipecraft/vsearch:2.18                vsearch v2.18, seqkit v2.3.0, GNU parallel                                    
====================================  ==========================================================================



