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
  :class: center

.. |stop_workflow| image:: _static/stop_workflow.png
  :width: 200
  :alt: Alternative text

.. |DADA2_PE_FWD_cp| image:: _static/DADA2_PE_FWD_cp.png
  :width: 700
  :alt: Alternative text

.. |cut_primers_expand_example| image:: _static/cut_primers_expand_example.png
  :width: 600
  :alt: Alternative text 

.. |DADA2_quality_filt_expand| image:: _static/DADA2_quality_filt_expand.png
  :width: 600
  :alt: Alternative text

.. |DADA2_denoise_expand| image:: _static/DADA2_denoise_expand.png
  :width: 600
  :alt: Alternative text

.. |DADA2_assign_tax_expand| image:: _static/DADA2_assign_tax_expand.png
  :width: 600
  :alt: Alternative text

.. |DADA2_filter_table_expand_COI| image:: _static/DADA2_filter_table_expand_COI.png
  :width: 600
  :alt: Alternative text

.. |DADA2_2samples_needed| image:: _static/troubleshoot/DADA2_2samples_needed.png
  :width: 300
  :alt: Alternative text

.. |ASVs_to_OTUs_COI| image:: _static/ASVs_to_OTUs.png
  :width: 650
  :alt: Alternative text

.. |LULU| image:: _static/lulu.png
  :width: 650
  :alt: Alternative text

.. |output_icon| image:: _static/output_icon.png
  :width: 50
  :alt: Alternative text

.. |save| image:: _static/save.png
  :width: 50
  :alt: Alternative text

.. |pulling_image| image:: _static/pulling_image.png
  :width: 280
  :alt: Alternative text

.. |COI_example_quality_plot| image:: _static/COI_example_quality_plot.png
  :width: 500
  :alt: Alternative text

.. meta::
    :description lang=en:
        PipeCraft manual. tutorial

|

DADA2 ASVs pipeline, COI |PipeCraft2_logo|
==========================================

This example dataset consists of **COI mtDNA gene amplicon sequences with the target length of 313 bp**:

| `Download example data set here <https://raw.githubusercontent.com/pipecraft2/user_guide/master/data/example_data_COI_313bp.zip>`_ and unzip it. 

For this example data run, we are using a subset of `CO1Classifier <https://github.com/terrimporter/CO1Classifier>`_ database in the taxonomy annotation process, `download it from here <https://raw.githubusercontent.com/pipecraft2/user_guide/master/data/Databases/SINTAX_COIv5.1.0.subset.zip>`_.


____________________________________________________

Starting point 
--------------

This example dataset consists of **COI mtDNA gene amplicon sequences with the target length of 313 bp**:

- **paired-end** Illumina MiSeq data;
- **demultiplexed** set (per-sample fastq files);
- primers are not not **removed**;
- sequences in this set are **5'-3' (fwd) oriented**.


.. admonition:: when working with your own data ...

  ... then please check that the paired-end data file names contain **R1** and **R2** strings *(not just _1 and _2)*, so that 
  PipeCraft can correctly identify the paired-end reads.

  | *Example:*
  | *sample1_COI.R1.fq.gz*
  | *sample1_COI.R2.fq.gz*

**At least 2 samples** (2x R1 + 2x R2 files) are required for this workflow! Otherwise ERROR in the denoising step:

|DADA2_2samples_needed| 

____________________________________________________

| **To select DADA2 pipeline**, press
| ``SELECT PIPELINE`` --> ``DADA" ASVs``.
| 
| **To select input data**, press ``SELECT WORKDIR``
| and specify
| ``sequence files extension`` as **\*.fastq**;  
| ``sequencing read types`` as **paired-end**.

____________________________________________________

Workflow mode
-------------

Because we are working with sequences that are **5'-3' oriented**, we are selecting hte ``PAIRED-END FORWARD`` mode of the pipeline. 

|DADA2_PE_FWD_cp| 

.. admonition:: if sequences are in mixed orientation
 
 If some sequences in your library are in 5'-3' and some as 3'-5' orientation, 
 then with the 'PAIRED-END FORWARD' mode exactly the same ASV may be reported twice, where one ASV is just the reverse complementary of another. 
 To avoid that, select **PAIRED-END MIXED** mode. 
 *Sequences have mixed orientation in libraries where sequenceing adapters have been ligated, rather than attached to amplicons during PCR.*

 **Specifying primers** (for CUT PRIMERS) **is mandatory for the PAIRED-END MIXED** mode. Based on the priemr sequences, the library will be split into two: 
 1) fwd oriented sequences, and 2) rev oriented sequences. Both batches are processed independently to produce ASVs, after which the rev oriented batch ASVs are 
 reverse complemented and merged with the fwd oriented ASVs. Identical ASVs are merged to form a final data set. This is a reccomended workflow for accurate denoising compared with first 
 reorienting all sequences to 5'-3', and then performing a standard 'PAIRED-END FORWARD' workflow.

____________________________________________________

Cut primers
-----------

The example dataset **contains primer sequences**. Generally, we need to remove these to proceed the analyses only with the variable metabarcode of interest.
If there are some additional sequence fragments, from eg. sequencing adapters or poly-G tails, then clipping the primers will remove those fragments as well.

Tick the box for ``CUT PRIMERS`` and specify forward and reverse primers.
For the example data, the **forward primer is GGWACWGGWTGAACWGTWTAYCCYCC** and **reverse primer is TANACYTCNGGRTGNCCRAARAAYCA**.

|cut_primers_expand_example|

The primers are 26 bp - to keep a bit of flexibility in the primer search, we are requesting the ``min overlap`` of **22 bp** and are allowing maximum of 2 ``mismatches`` . 
Note that too low ``min overlap`` may lead to random matches. Check :ref:`other CUT PRIMER options here <remove_primers>`.


__________________________________________________

Quality filtering 
-----------------

Before adjusting quality filtering settings, let's have a look on the **quality profile** of our example data set. 
Below quality profile plot was generated using ``QualityCheck`` panel (:ref:`see here <qualitycheck>`).

|COI_example_quality_plot|

All files are represented with **green lines, indicating good average quality per file** (i.e., sample). 
However, if you see lower qualities of especially towards the end of R2 reads, then it not too alarming, since 
those can be clipped off with ``truncLen R2`` setting. DADA2 algoritm is robust to lower quality sequences, 
but removing the low quality read parts will improve the DADA2 sensitivity to rare sequence variants. 
But herein, we do not need to clip the ends, because the overall quality of the sequences is good enough.

____________________________________________________

**Click on** ``QUALITY FILTERING`` **to expand the panel**

.. |COI_ex_qFilt.png| image:: _static/COI_ex_qFilt.png
  :width: 600
  :alt: Alternative text

|COI_ex_qFilt.png|

Here, we can leave the settings as DEFAULT by discarding sequences with 
**maximum error rate of >2** and with **ambiguous bases of >0**. 

+-----------------------+-------------------------------------------------------+
| Output directory |output_icon|          ``qualFiltered_out``                  |
+=======================+=======================================================+
| \*.fq.gz              | quality filtered sequences per sample in FASTQ format |
+-----------------------+-------------------------------------------------------+
| \*.rds                | R objects for the following DADA2 workflow processes  |
+-----------------------+-------------------------------------------------------+
| seq_count_summary.csv | summary of sequence counts per sample                 |
+-----------------------+-------------------------------------------------------+

____________________________________________________

Denoise and merge pairs
-----------------------

This step performs desiosing (as implemented in DADA2), which first forms ASVs per R1 and R2 files. 
Then during merging/assembling process the paired ASV mates are assembled to output full amplicon length ASV. 

|DADA2_denoise_expand| 

Here, we are working with Illumina data, so let's make sure that the ``errorEstFun`` setting is **loessErrfun**. 
We can leave all settings as DEFAULT. 

+----------------------------------+--------------------------------------------------------+
| Output directory |output_icon|          ``denoised_assembled.dada2``                      |
+==================================+========================================================+
| \*.fasta                         | denoised and assembled ASVs per sample in FASTA format |
+----------------------------------+--------------------------------------------------------+
| \*.rds                           | R objects for the following DADA2 workflow processes   |
+----------------------------------+--------------------------------------------------------+
| Error_rates_R*.pdf               | plots for estimated R1/R2 error rates                  |
+----------------------------------+--------------------------------------------------------+
| seq_count_summary.csv            | summary of sequence counts per sample                  |
+----------------------------------+--------------------------------------------------------+

___________________________________________________

Chimera filtering
-----------------

This step performs chimera filtering according to DADA2 removeBimeraDenovo function. During this step, the **ASV table** is also generated. 

.. important:: 

  make sure that primers have been removed from your amplicons; otherwise many false-positive chimeras may be filtered out from your dataset. 

Here, we filter chimeras using the **consensus** method. 

+----------------------------------------+-------------------------------------------------------------------+
| Output directory                       | ``chimeraFiltered_out.dada2``                                     |
+========================================+===================================================================+
| \*.fasta                               | chimera filtered ASVs per sample                                  |
+----------------------------------------+-------------------------------------------------------------------+
| seq_count_summary.csv                  | summary of sequence counts per sample                             |
+----------------------------------------+-------------------------------------------------------------------+
| 'chimeras' dir                         | ASVs per sample identified as chimeras                            |
+----------------------------------------+-------------------------------------------------------------------+
| Output directory                       | ``ASVs_out.dada2``                                                |
+----------------------------------------+-------------------------------------------------------------------+
| ASVs_table.txt                         | denoised and chimera filtered ASV-by-sample table                 |
+----------------------------------------+-------------------------------------------------------------------+
| ASVs.fasta                             | corresponding FASTA formated ASV Sequences                        |
+----------------------------------------+-------------------------------------------------------------------+
| ASVs per sample identified as chimeras | rds formatted denoised and chimera filtered ASV table (for DADA2) |
+----------------------------------------+-------------------------------------------------------------------+

____________________________________________________


Curate ASV table
----------------

This process first removes putative **tag jumps** and then **collapses the ASVs that are identical** up to shifts or length variation, 
i.e. ASVs that have no internal mismatches; and finally 
filters out ASVs that are shorter/longer than specified length (in base pairs).

Here, we are **enabling this process** by checking the box for ``CURATE ASV TABLE`` in the DADA2 ASV workflow panel. 

|DADA2_filter_table_expand_COI|

The ``f_value`` and ``p_value`` settings are used to filter out putative tag jumps (using UNCROSS2 algorithm). 
Generally, we recommend to use p_value of 1 (default), and **f_value of 0.03** when using combinational indexing strategy; 
f_value of 0.05 when using single-indexes, and f_value of 0.01 when using unique dual-indexes.

The expected amplicon length (without primers) in our example dataset in **313 bp**.
Assuming that shorter sequences are non-target sequences, 
we use 307 in the ``min length`` setting and 319 in ``max length`` setting. 
This will discard ASVs that are shorter than 307 bp or longer than 319 bp.

We are also setting the ``collapseNoMismatch`` to TRUE, to collapse identical ASVs. 
This is basically equivalent to 100% clustering by ignoring the end gaps.

+----------------------------+-------------------------------------------------------------------+
| Output directory |output_icon|            ``ASVs_out.dada2/curated``                           |
+============================+===================================================================+
|| ASVs_collapsed.fasta      || tag-jump filtered and collapsed and size filtered                |
||                           || ASV Sequences                                                    |
+----------------------------+-------------------------------------------------------------------+
| ASVs_table_collapsed.txt   | corresponding ASV-by-sample table (curated)                       |
+----------------------------+-------------------------------------------------------------------+
| ASVs_table_TagJumpFilt.txt | only tag-jump filtered ASV-by-sample table                        |
+----------------------------+-------------------------------------------------------------------+
| ASVs.fasta                 | corresponding ASV Sequences with ASVs_table_TagJumpFilt.txt table |
+----------------------------+-------------------------------------------------------------------+

.. admonition:: If there is nothing to collapse or filter out based on the length
  
  then there are no corresponding files in the ``ASVs_out.dada2/curated`` directory, and only 
  ASVs_table_TagJumpFilt.txt and ASVs.fasta files will be generated 
  (even when there is nothing to tag-jump filter - in which case ASVs_table_TagJumpFilt.txt is the same 
  ASVs_table.txt in the ``ASVs_out.dada2`` directory).


.. note:: 

  The pre-defined pipeline ends here. Outputs COI ASVs should be further filtered to **remove 
  putative pseudogenes (NUMTs)**, and optionally ASVs can be clustered into OTUs. See below. 

____________________________________________________

Save workflow
-------------

Once we have decided about the settings in our workflow, we can save the configuration file by pressing ``save workflow`` button on the right-ribbon
|save|

If you forget the save, then no worries, a ``pipecraft2_last_run_configuration.json`` file will be generated for you upon starting the workflow.
As the file name says, it is the workflow configuration file for your last PipeCraft run in this **working directory**. 

This ``JSON`` file can be loaded into PipeCraft2 to **automatically configure your next runs exactly the same way**.

___________________________________________________

Start the workflow
------------------

Press ``START`` on the left ribbon **to start the analyses**.

.. admonition:: when running the module for the first time ...
  
  ... a docker image will be first pulled to start the process. 

  For example: |pulling_image|


When you need to STOP the workflow, press ``STOP`` button |stop_workflow|


.. admonition:: When the workflow has completed ...

  ... a message window will be displayed.

  |workflow_finished|

___________________________________________________



Additional processes via QuickTools
-----------------------------------

The following are **additional steps** that can be performed, with **removing putative NUMTs as the most important one** for COI data.

__________________________________________________

Cluster ASVs into OTUs
~~~~~~~~~~~~~~~~~~~~~~

If the aim is not to do the haplotype-level analyses, then ASVs can be clustered into OTUs using vsearch.
The ASVs approach may not accurately reflect species composition in the community of as 
COI gene has highly variable levels of intraspecific polymorphism. Thus, 
one species may be represented by many ASVs, whereas other species may be represented by very few ASVs.

Here, we are **clustering ASVs to OTUs** (using vsearch) via **QuickTools** panel (on the right ribbon).

|ASVs_to_OTUs_COI|

| **Input data** are ``ASVs_table_collapsed.txt`` and ``ASVs_collapsed.fasta`` files in the ``ASVs_out.dada2/curated`` directory.

.. note::

  2nd column of **ASV table must be 'Sequences'** (1st column is ASV IDs; default PipeCraft2 output table from the previous step).
  For clustering, the ASV size annotation is obtained from the ASV table. 


+---------------------------+--------------------------------------------------------+
| Outputs ``ASVs2OTUs_out`` |                                                        |
+===========================+========================================================+
| OTUs.fasta                | FASTA formated representative OTU sequences            |
+---------------------------+--------------------------------------------------------+
| OTU_table.txt             | OTU distribution table per sample (tab delimited file) |
+---------------------------+--------------------------------------------------------+
| OTUs.uc                   | uclust-like formatted clustering results for OTUs      |
+---------------------------+--------------------------------------------------------+
| ASVs.size.fasta           | size annotated input sequences                         |
+---------------------------+--------------------------------------------------------+


___________________________________________________

LULU post-clustering
~~~~~~~~~~~~~~~~~~~~

Additionally, we can perform `LULU post-clustering <https://github.com/tobiasgf/lulu>`_ to merge co-occurring 'daughter' OTUs.

LULU description from the `LULU repository <https://github.com/tobiasgf/lulu>`_: the purpose of LULU is to reduce the number of 
erroneous OTUs in OTU tables to achieve more realistic biodiversity metrics. 
By evaluating the co-occurence patterns of OTUs among samples LULU identifies OTUs that consistently satisfy some user selected 
criteria for being errors of more abundant OTUs and merges these. **It has been shown that curation with LULU consistently result 
in more realistic diversity metrics. **

|LULU| 

The **input data** are ``OTU_table.txt`` and ``OTUs.fasta`` files in the ``ASVs2OTUs_out`` directory.

+-----------------------+----------------------------------------------------------------------------+
| Outputs ``lulu_out``  |                                                                            |
+=======================+============================================================================+
| OTU_table.lulu.txt    | curated table in tab delimited txt format                                  |
+-----------------------+----------------------------------------------------------------------------+
| OTUs.lulu.fasta       | fasta file for the molecular units (OTUs or ASVs) in the curated table     |
+-----------------------+----------------------------------------------------------------------------+
| match_list.lulu       | match list file that was used by LULU to merge 'daughter' molecular units  |
+-----------------------+----------------------------------------------------------------------------+
|| discarded_units.lulu || molecular units (OTUs or ASVs) that were merged with other units based on |
||                      || specified thresholds                                                      |
+-----------------------+----------------------------------------------------------------------------+

__________________________________________________

Remone NUMTs
~~~~~~~~~~~~

When working with processing **protein-coding markers** (such as COI),
the amplified sequences of **nuclear mitochondrial pseudogenes (NUMTs)** 
may inflate the richness estimates and thus introduce biases in biodiversity research using metabarcoding.
Therefore, it is important to remove these sequences from the data set.

In PipeCraft2, there are tools such as `**metaMATE** <https://github.com/tjcreedy/metamate>`_ and `**ORFfinder** <https://www.ncbi.nlm.nih.gov/orffinder/>`_ 
that can be used to remove NUMTs from the data set.


__________________________________________________



Assign taxonomy
---------------

Assign taxonomy **is not the part of the full per-defined pipeline**, but can be run as a **separate step in QuickTools**.
Here, we are using the :ref:`DADA2 classifier <assign_taxonomy_dada2>`, the assignTaxonomy function, 
which itself implements the RDP Naive Bayesian Classifier algorithm. 
See :ref:`other assign taxonomy options here <assign_taxonomy>`.

We need to specify the location of the **reference DATABASE** for the taxonomic classification of our ASVs. Click on the header of ``dada2_database`` setting, 
which directs you to the `DADA2-formatted reference databases web page <https://benjjneb.github.io/dada2/training.html>`_.
Here, we are using ``silva_nr99_v138.2_toSpecies_trainset.fa.gz``. 

|DADA2_assign_tax_expand|

Specify the location of your downloaded DADA2 database by pressing ``SELECT FASTA``. 
The default minBoot (minimum bootstrap; ranging from 0-100; ~assignment confidence) is 50, but here we are setting this to **80**. 
This means that taxonomic ranks with at least bootstrap value of 80 will get classification (unclassified for <80). 

``tryRC`` may be OFF, since we expact that all of our ASVs are in 5'-3' orientation. 
:ref:`See DADA2 assign taxonomy settings here <assign_taxonomy_dada2>`

+-------------------------------------------------------------+
| Output directory   |output_icon| ``taxonomy_out.dada2``     |
+==================+==========================================+
| taxonomy.csv     | classifier results with bootstrap values |
+------------------+------------------------------------------+



