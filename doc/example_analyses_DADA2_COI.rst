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

.. |DADA2_PE_FWD| image:: _static/DADA2_PE_FWD.png
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

.. |DADA2_filter_table_expand| image:: _static/DADA2_filter_table_expand.png
  :width: 600
  :alt: Alternative text

.. |DADA2_2samples_needed| image:: _static/troubleshoot/DADA2_2samples_needed.png
  :width: 300
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

.. meta::
    :description lang=en:
        PipeCraft manual. tutorial

|

:red:`THIS IS UNCOMPLETE PAGE`


COI example pipeline |PipeCraft2_logo|
--------------------------------------

This example dataset consists of **COI mtDNA gene amplicon sequences with the target length of 313 bp**:

| `Download example data set here <https://raw.githubusercontent.com/pipecraft2/user_guide/master/data/example_data_COI_313bp.zip>`_ and unzip it. 

For this example data run, we are using a subset of `CO1Classifier <https://github.com/terrimporter/CO1Classifier>`_ database in the taxonomy annotation process, `download it from here <https://raw.githubusercontent.com/pipecraft2/user_guide/master/data/Databases/SINTAX_COIv5.1.0.subset.zip>`_.


____________________________________________________

Starting point 
~~~~~~~~~~~~~~

This example dataset consists of **COI mtDNA gene amplicon sequences with the target length of 313 bp**:

- **paired-end** Illumina MiSeq data;
- **demultiplexed** set (per-sample fastq files);
- primers are not not **removed**;
- sequences in this set are **5'-3' (fwd) oriented**.


.. admonition:: when working with your own data ...

  ... then please check that the paired-end data file names contain **R1** and **R2** strings *(not just _1 and _2)*, so that 
  PipeCraft can correctly identify the paired-end reads.

  | *Example:*
  | *F3D0_S188_L001_R1_001.fastq*
  | *F3D0_S188_L001_R2_001.fastq*

  
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
~~~~~~~~~~~~~

Because we are working with sequences that are **5'-3' oriented**, we are selecting hte ``PAIRED-END FORWARD`` mode of the pipeline. 

|DADA2_PE_FWD| 

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
~~~~~~~~~~~

The example dataset **contains primer sequences**. Generally, we need to remove these to proceed the analyses only with the variable metabarcode of interest.
If there are some additional sequence fragments, from eg. sequencing adapters or poly-G tails, then clipping the primers will remove those fragments as well.

Tick the box for ``CUT PRIMERS`` and specify forward and reverse primers.
For the example data, the **forward primer is GGWACWGGWTGAACWGTWTAYCCYCC** and **reverse primer is TANACYTCNGGRTGNCCRAARAAYCA**.

|cut_primers_expand_example|

The primers are 26 bp - to keep a bit of flexibility in the primer search, we are requesting the ``min overlap`` of **22 bp** and are allowing maximum of 2 ``mismatches`` . 
Note that too low ``min overlap`` may lead to random matches. Check :ref:`other CUT PRIMER options here <remove_primers>`.


__________________________________________________

**THIS IS UNCOMPLETE PAGE**
