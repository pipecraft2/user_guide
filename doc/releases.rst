.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text


|PipeCraft2_logo|
  `github <https://github.com/SuvalineVana/pipecraft>`_
 

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red


.. _releases:

=========
Releases
=========

.. contents:: Contents
   :depth: 2

____________________________________________________

 
.. _0.1.2:

0.1.2 (07.06.2022)
==================

* added LULU post-clustering 
* added DEICODE (postprocessing)
* added fastp quality filtering
* added DADA2 assignTaxonomy under genereal 'ADD STEP' -> 'ASSIGN TAXONOMY' panel
* added --fastq_truncqual option for vsearch quality filtering
    
Implemented software:
*(software in red font denote new additions; 'version' in bold denotes version upgrade)*

=======================================================================  ========  =========================================================================================
Software                                                                 version   Reference                                                                                  
=======================================================================  ========  =========================================================================================
`DADA2 <https://benjjneb.github.io/dada2/index.html>`_                   **1.20**  `Callahan et. al 2016 <https://www.nature.com/articles/nmeth.3869>`_                      
`vsearch <https://github.com/torognes/vsearch>`_                         2.18.0    `Rognes et. al 2016 <https://peerj.com/articles/2584/>`_                                  
`trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_          0.39      `Bolger et al. 2014 <https://doi.org/10.1093/bioinformatics/btu170>`_                     
`seqkit <https://bioinf.shenwei.me/seqkit/>`_                            2.0.0     `Shen et al. 2016 <https://doi.org/10.1371/journal.pone.0163962>`_                        
`cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_                 3.5       `Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_                                     
`mothur <https://github.com/mothur/mothur>`_                             1.46.1    `Schloss et al. 2009 <https://doi.org/10.1128/AEM.01541-09>`_                             
`ITS Extractor <https://microbiology.se/software/itsx/>`_                1.1.3     `Bengtsson-Palme et al. 2013 <https://doi.org/10.1111/2041-210X.12073>`_                  
`fqgrep <https://github.com/indraniel/fqgrep>`_                          0.4.4     `Indraniel Das 2011 <https://github.com/indraniel/fqgrep>`_                               
`BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_                      2.11.0+   `Camacho et al. 2009 <https://doi.org/10.1186/1471-2105-10-421>`_                         
`FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_   0.11.9    `Andrews 2019 <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_              
`MultiQC <https://multiqc.info/>`_                                       1.12      `Ewels et al. 2016 <https://doi.org/10.1093/bioinformatics/btw354>`_                      
:red:`LULU` `(link) <https://github.com/tobiasgf/lulu>`_                 0.1.0     `Froslev et al. 2017 <https://doi.org/10.1038/s41467-017-01312-x>`_
:red:`fastp` `(link) <https://github.com/OpenGene/fastp>`_               0.23.2    `Chen et al. 2018 <https://doi.org/10.1093/bioinformatics/bty560>`_
:red:`DEICODE` `(link) <https://github.com/biocore/DEICODE>`_            0.2.4     `Martion et al. 2019 <https://journals.asm.org/doi/10.1128/mSystems.00016-19>`_
=======================================================================  ========  =========================================================================================

____________________________________________________

.. _0.1.1:

0.1.1 (01.04.2022)
==================

Minor cosmetic changes and bug fixes. 
`DOWNLOADS <https://github.com/SuvalineVana/pipecraft/releases/tag/0.1.1>`_

* separate output forlder for unused index combinations in demultiplexing.  
* resolved issues with sample renaiming when using dual combinational indexes for paired-end data 
  (DEMULTIPLEX)
* minBoot option fixed in DADA2 taxonomy annotation
* vsearch quality filtering "minsize" not working (option currently removed).

____________________________________________________

.. _0.1.0:

0.1.0 pre-release (14.12.2021)
==============================

`DOWNLOADS <https://github.com/SuvalineVana/pipecraft/releases/tag/0.1.0>`_

* ASV workflow with DADA2 for paired-end data.
* vsearch based OTU workflow.
* QualityCheck module with MultiQC and FastQC

Implemented software:

=======================================================================  ========  =========================================================================================
Software                                                                 version   Reference                                                                                  
=======================================================================  ========  =========================================================================================
`DADA2 <https://benjjneb.github.io/dada2/index.html>`_                   1.14      `Callahan et. al 2016 <https://www.nature.com/articles/nmeth.3869>`_                      
`vsearch <https://github.com/torognes/vsearch>`_                         2.18.0    `Rognes et. al 2016 <https://peerj.com/articles/2584/>`_                                  
`trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_          0.39      `Bolger et al. 2014 <https://doi.org/10.1093/bioinformatics/btu170>`_                     
`seqkit <https://bioinf.shenwei.me/seqkit/>`_                            2.0.0     `Shen et al. 2016 <https://doi.org/10.1371/journal.pone.0163962>`_                        
`cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_                 3.5       `Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_                                     
`mothur <https://github.com/mothur/mothur>`_                             1.46.1    `Schloss et al. 2009 <https://doi.org/10.1128/AEM.01541-09>`_                             
`ITS Extractor <https://microbiology.se/software/itsx/>`_                1.1.3     `Bengtsson-Palme et al. 2013 <https://doi.org/10.1111/2041-210X.12073>`_                  
`fqgrep <https://github.com/indraniel/fqgrep>`_                          0.4.4     `Indraniel Das 2011 <https://github.com/indraniel/fqgrep>`_                               
`BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_                      2.11.0+   `Camacho et al. 2009 <https://doi.org/10.1186/1471-2105-10-421>`_                         
`FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_   0.11.9    `Andrews 2019 <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_              
`MultiQC <https://multiqc.info/>`_                                       1.12      `Ewels et al. 2016 <https://doi.org/10.1093/bioinformatics/btw354>`_                      
=======================================================================  ========  =========================================================================================
