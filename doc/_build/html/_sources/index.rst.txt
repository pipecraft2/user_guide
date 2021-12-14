.. |PipeCraft2_logo| image:: _static/PipeCraft2_logo.png
  :width: 150
  :alt: Alternative text

.. |add_step_interface| image:: _static/add_step.png
  :width: 600
  :alt: Alternative text

|PipeCraft2_logo|

==========================================
Manual for PipeCraft 2
==========================================

| **PipeCraft** in a Graphical User Interfece software that implements :ref:`various popular tools <tools>` for **metabarcoding** data analyzes that are linked to generate custom bioinformatics pipeline/workflow. 
| Pre-defined full pipelines for generating :ref:`OTUs <otupipe>` or :ref:`ASVs <asvpipe>` are also implemented.

|add_step_interface|

| :ref:`Panels <panels>` for pipeline processes contain key options for sequence data analyses, but full use of any implemented program may be accessed via :ref:`PipeCraft console (command line) <expert_mode>`. 
| Default settings in the panels represent commonly used options for amplicon sequence data analyses, which may be tailored according to user experience or needs. 
 Custom designed pipeline settings can be saved and thus the exact same pipeline may be easily re-run on other sequencing data (and for reproducibility, may be used as a supplement material in the manuscript). 
 PipeCraft enables generation of the full pipeline (user specifies the input data and output will be e.g. OTU/ASV table with taxonomic annotations of the OTUs/ASVs), 
 but supports also single-step mode where analyses may be performed in a step-by-step manner *(e.g. perform quality filtering, then examine the output and decide whether to adjust the quality filtering options of 
 to proceed with next step, e.g. with chimera filtering step)*.

.. _tools:

Currently implemented software
------------------------------

=======================================================================  =========================================================================================  =============
Software                                                                 Reference                                                                                  Task
=======================================================================  =========================================================================================  =============
`docker <https://www.docker.com/>`_                                      https://www.docker.com                                                                     building, sharing and running applications
`DADA2 <https://benjjneb.github.io/dada2/index.html>`_                   `Callahan et. al 2016 <https://www.nature.com/articles/nmeth.3869>`_                       ASVs workflow (from raw reads to ASV table)
`vsearch <https://github.com/torognes/vsearch>`_                         `Rognes et. al 2016 <https://peerj.com/articles/2584/>`_                                   quality filtering, assemble paired-end reads, chimera filtering, clustering
`trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_          `Bolger et al. 2014 <https://doi.org/10.1093/bioinformatics/btu170>`_                      quality filtering
`seqkit <https://bioinf.shenwei.me/seqkit/>`_                            `Shen et al. 2016 <https://doi.org/10.1371/journal.pone.0163962>`_                         multiple sequence manipulation operations
`cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_                 `Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_                                      demultiplexing, cut primers
`biopython <https://biopython.org/>`_                                    `Cock et al. 2009 <https://academic.oup.com/bioinformatics/article/25/11/1422/330687>`_    multiple sequence manipulation operations
`GNU Parallel <https://doi.org/10.5281/zenodo.4710607>`_                 `Tangle 2021 <https://doi.org/10.5281/zenodo.4710607>`_                                    executing jobs in parallel
`mothur <https://github.com/mothur/mothur>`_                             `Schloss et al. 2009 <https://doi.org/10.1128/AEM.01541-09>`_                              submodule in ITSx to make unique and deunique seqs
`ITS Extractor <https://microbiology.se/software/itsx/>`_                `Bengtsson-Palme et al. 2013 <https://doi.org/10.1111/2041-210X.12073>`_                   extract ITS regions
`fqgrep <https://github.com/indraniel/fqgrep>`_                          `Indraniel Das 2011 <https://github.com/indraniel/fqgrep>`_                                core for reorient reads
`BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_                      `Camacho et al. 2009 <https://doi.org/10.1186/1471-2105-10-421>`_                          assign taxonomy
`FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_   `Andrews 2019 <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_               QualityCheck module
`MultiQC <https://multiqc.info/>`_                                       `Ewels et al. 2016 <https://doi.org/10.1093/bioinformatics/btw354>`_                       QualityCheck module
=======================================================================  =========================================================================================  =============

Let us know if you would like to have a specific software implemeted to PipeCraft; :ref:`contacts <contact>`


Contents
--------

.. toctree::
   :maxdepth: 2

   installation
   user_guide
   licence
   contact
   citation

| 
| 

*Manual may contain some typos! Fixing those on the way.*
