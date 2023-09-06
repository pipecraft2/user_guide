.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text


|PipeCraft2_logo|
  `github <https://github.com/pipecraft2/pipecraft>`_
 

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red


.. _releases:

========
Releases
========

.. contents:: Contents
   :depth: 2

____________________________________________________

.. hide 

    for next release - BLAST dabasese resource link to GUI
    xxx





.. _1.0.0:

1.0.0 (01.09.2023)
==================

`DOWNLOAD link for v1.0.0 <XXX>`_

* major updates in the front-end; individual tools on the right, pipelines on left.
* added debugging mode and improved log info 
* added NextITS pipeline for PacBio ITS sequences (not available for MacOS release)
* added ORFfinder + HMM bsed pseudogene/off-targets filtering for protein coding genes
* added RDP classifier
* added DADA2 pipeline for PacBio data 
* added DADA2 pipeline for paired-end mixed oriented amplicons (fwd_orient and rev_orient are denoised separately and then merged)
* implemented DADA2 denoising sensitivity editing
* all features will get sha1 ID 
* added ASVs to OTUs module (cluster ASVs into OTUs with vsearch)
* added tag-jumps filtering module (UNCROSS2)
* fixed the vsearch_dada2 container issues for MacOS 
  
Implemented software:
*(software in red font denote new additions; 'version' in bold denotes version upgrade)*

=======================================================================  ==========
Software                                                                 version                                                                                       
=======================================================================  ==========
:red:`NextITS pipeline` `(link) <https://next-its.github.io/>`_          **0.5.0**
:red:`ORFfinder` `(link) <https://www.ncbi.nlm.nih.gov/orffinder/>`_     **v0.4.3**
:red:`RDP classifier`                                                    **v2.13**
`DADA2 <https://benjjneb.github.io/dada2/index.html>`_                   **1.27**
`vsearch <https://github.com/torognes/vsearch>`_                         **2.23**
`trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_          0.39
`seqkit <https://bioinf.shenwei.me/seqkit/>`_                            2.3.0
`cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_                 **4.4**
`mothur <https://github.com/mothur/mothur>`_                             1.46.1
`ITS Extractor <https://microbiology.se/software/itsx/>`_                1.1.3
`fqgrep <https://github.com/indraniel/fqgrep>`_                          0.4.4
`BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_                      **2.12**
`FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_   0.11.9
`MultiQC <https://multiqc.info/>`_                                       1.12
`LULU <https://github.com/tobiasgf/lulu>`_                               0.1.0
`fastp <https://github.com/OpenGene/fastp>`_                             0.23.2
`DEICODE <https://github.com/biocore/DEICODE>`_                          0.2.4
=======================================================================  ==========

____________________________________________________

.. _0.1.4:

0.1.4 (15.12.2022)
==================

`DOWNLOAD link for v0.1.4 <https://github.com/pipecraft2/pipecraft/releases/tag/v0.1.4>`_

* added 2nd round of cut primers to properly remove fwd and rev primers form the paired-end data set
* added UNOISE3 module to generate zOTUs (under clustering)
* added uchime3 chimera filtering (for denoised amplicons)
* edited sequence count statistics process after the process (using seqkit)
* only fasta (fa, fas) format is accepted for clustering
* edited OTU table making strategy for OTU clustering (was --usearch_global before)
* added table filtering options for DADA2 ASV table (collapse mismatch, filter by length)
* added ASV to OTU module (clustering DADA2 ASVs into OTUs)
* select region to cluster after ITSx in OTUs workflow
* automatically saves the PipeCraft workflow settings into loadable JSON file
* outputs log file (in development)
* merged vsearch and dada2 containers (had a lot in common)
  
Implemented software:
*(software version in bold denotes version upgrade)*

=======================================================================  ==========  =========================================================================================
Software                                                                 version     Reference                                                                                  
=======================================================================  ==========  =========================================================================================
`DADA2 <https://benjjneb.github.io/dada2/index.html>`_                   1.20        `Callahan et. al 2016 <https://www.nature.com/articles/nmeth.3869>`_                      
`vsearch <https://github.com/torognes/vsearch>`_                         **2.22.1**  `Rognes et. al 2016 <https://peerj.com/articles/2584/>`_                                  
`trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_          0.39        `Bolger et al. 2014 <https://doi.org/10.1093/bioinformatics/btu170>`_                     
`seqkit <https://bioinf.shenwei.me/seqkit/>`_                            **2.3.0**   `Shen et al. 2016 <https://doi.org/10.1371/journal.pone.0163962>`_                        
`cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_                 3.5         `Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_                                     
`mothur <https://github.com/mothur/mothur>`_                             1.46.1      `Schloss et al. 2009 <https://doi.org/10.1128/AEM.01541-09>`_                             
`ITS Extractor <https://microbiology.se/software/itsx/>`_                1.1.3       `Bengtsson-Palme et al. 2013 <https://doi.org/10.1111/2041-210X.12073>`_                  
`fqgrep <https://github.com/indraniel/fqgrep>`_                          0.4.4       `Indraniel Das 2011 <https://github.com/indraniel/fqgrep>`_                               
`BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_                      2.11.0+     `Camacho et al. 2009 <https://doi.org/10.1186/1471-2105-10-421>`_                         
`FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_   0.11.9      `Andrews 2019 <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_              
`MultiQC <https://multiqc.info/>`_                                       1.12        `Ewels et al. 2016 <https://doi.org/10.1093/bioinformatics/btw354>`_                      
`LULU <https://github.com/tobiasgf/lulu>`_                               0.1.0       `Froslev et al. 2017 <https://doi.org/10.1038/s41467-017-01312-x>`_
`fastp <https://github.com/OpenGene/fastp>`_                             0.23.2      `Chen et al. 2018 <https://doi.org/10.1093/bioinformatics/bty560>`_
`DEICODE <https://github.com/biocore/DEICODE>`_                          0.2.4       `Martion et al. 2019 <https://journals.asm.org/doi/10.1128/mSystems.00016-19>`_
=======================================================================  ==========  =========================================================================================

____________________________________________________


.. _0.1.3:

0.1.3 (28.07.2022)
==================

`DOWNLOAD link for v0.1.3 <https://github.com/SuvalineVana/pipecraft/releases/tag/v0.1.3>`_

* updated BLAST 2.11.0+ to BLAST 2.12.0+ and added biopython to BLAST container (fixed the coverage% calculation)
* fixed the megaBLAST, when gapextend=undefined
* quality Check module edit (does not stop when browsing around)
* fixed ASVs workflow error message when using <2 samples
* added lock panels when starting a process
* few cosmetic front-end adds  

.. _0.1.2:

0.1.2 (07.06.2022)
==================

`DOWNLOAD link for v0.1.2 <https://github.com/SuvalineVana/pipecraft/releases/tag/v0.1.2>`_

* added LULU post-clustering 
* added DEICODE (postprocessing)
* added fastp quality filtering
* added DADA2 quality filtering under 'ADD STEP' -> 'QUALITY FILTERING' panel
* added DADA2 denoise and assemble paired-end data under 'ADD STEP' -> 'ASSEMBLE PAIRED-END' panel
* added DADA2 assignTaxonomy under 'ADD STEP' -> 'ASSIGN TAXONOMY' panel
* added trunc_length option for vsearch quality filtering
* python3 module fix for ITSx for removing empty sequeces 
    
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
`DOWNLOAD link for v0.1.1 <https://github.com/SuvalineVana/pipecraft/releases/tag/0.1.1>`_

* separate output forlder for unused index combinations in demultiplexing.  
* resolved issues with sample renaiming when using dual combinational indexes for paired-end data 
  (DEMULTIPLEX)
* minBoot option fixed in DADA2 taxonomy annotation
* vsearch quality filtering "minsize" not working (option currently removed).

____________________________________________________

.. _0.1.0:

0.1.0 pre-release (14.12.2021)
==============================

`DOWNLOAD link for v0.1.0 <https://github.com/SuvalineVana/pipecraft/releases/tag/0.1.0>`_

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
