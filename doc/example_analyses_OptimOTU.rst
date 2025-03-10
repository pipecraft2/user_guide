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

.. |workflow_finished| image:: _static/workflow_finished.png
  :width: 300
  :alt: Alternative text

.. |stop_workflow| image:: _static/stop_workflow.png
  :width: 200
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


OptimOTU pipeline, ITS2 |PipeCraft2_logo|
-----------------------------------------

This example data analyses follows OptimOTU workflow as implemented in PipeCraft2's pre-compiled pipelines panel. 

| `Download example data set here <https://raw.githubusercontent.com/pipecraft2/user_guide/master/data/example_data_ITS2.zip>`_ (15.1 Mb) and unzip it. 
| This is **ITS2 Illumina MiSeq** dataset. 

____________________________________________________

Starting point 
~~~~~~~~~~~~~~

This example dataset consists of **ITS2 rRNA gene amplicon sequences**; targeting fungi:

- **paired-end** Illumina MiSeq data;
- **demultiplexed** set (per-sample fastq files);
- primers **are not removed**;
- sequences in this set are **5'-3' (fwd) oriented**.


.. code-block::
   :caption: Required directory structure for OptimOTU

    my_dir/   
    └── sequences/         # SELECT THIS FOLDER AS WORKING DIRECTORY (name here can be anything)
        └── 01_raw/
            ├── Run1/      # name here can be anything (without spaces)
            │   ├── sample1_R1.fastq.gz
            │   ├── sample1_R2.fastq.gz
            │   ├── sample2_R1.fastq.gz
            │   └── sample2_R2.fastq.gz
            ├── Run2/      # name here can be anything (without spaces)
            │   ├── sample3_R1.fastq.gz
            │   ├── sample3_R2.fastq.gz
            │   ├── sample4_R1.fastq.gz
            │   └── sample4_R2.fastq.gz
            └── Run3/      # name here can be anything (without spaces)
                ├── sample5_R1.fastq.gz
                └── sample5_R2.fastq.gz

____________________________________________________


| **To select OptimOTU pipeline**, press
| ``SELECT PIPELINE`` --> ``OptimOTU``.
| 
| **To select input data**, press ``SELECT WORKDIR``
| and specify
| ``sequence files extension`` as **\*.fastq.gz**;  
| ``sequencing read types`` as **paired-end**.

___________________________________________________


Target taxa and sequence orientation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we are specifying that target taxa is **fungi**, and sequence orientation is **fwd**.

__________________________________________________

Control sequences
~~~~~~~~~~~~~~~~~

Control sequences are sequences that are not target taxa, but are used to estimate the error rate of the sequencing.

__________________________________________________


Cut primers and trim reads
~~~~~~~~~~~~~~~~~~~~~~~~~~

The example dataset **contains primer sequences**. Generally, we need to remove these to proceed the analyses only with the variable metabarcode of interest.
If there are some additional sequence fragments, from eg. sequencing adapters or poly-G tails, then clipping the primers will remove those fragments as well.

For the example data, the **forward primer is fITS7 GTGARTCATCGAATCTTTG** and **reverse primer is ITS4 TCCTCCGCTTATTGATATGC**.

  
____________________________________________________


Quality filtering 
~~~~~~~~~~~~~~~~~

Quality filtering here removes sequences which does not meet the threshold for the allowed maximum number of expected errors. 
See :ref:`here for more inforamtion about sequence quality <qualitycheck>` 
and `here for the additional information about expected errors <https://drive5.com/usearch/manual/exp_errs.html>`_.


____________________________________________________

Denoising and merging paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

...

__________________________________________________

Chimera filtering
~~~~~~~~~~~~~~~~~

...

__________________________________________________

Filter tag-jumps
~~~~~~~~~~~~~~~~~

...

__________________________________________________

Amplicon model setting
~~~~~~~~~~~~~~~~~~~~~~

``model_type`` = CM
``model_file`` = gITS7_ITS4.cm

__________________________________________________

Protax classification
~~~~~~~~~~~~~~~~~~~~~

``location`` = protaxFungi
``with_outgroup`` = UNITE_SHs

__________________________________________________  


Clustering
~~~~~~~~~~

``cluster thresholds`` = Fungi_GSSP

__________________________________________________

Save workflow
~~~~~~~~~~~~~

Once we have decided about the settings in our workflow, we can save the configuration file by pressing ``save workflow`` button on the right-ribbon
|save|

If you forget the save, then no worries, a ``pipecraft2_last_run_configuration.json`` file will be generated for you upon starting the workflow.
As the file name says, it is the workflow configuration file for your last PipeCraft run in this **working directory**. 

This ``JSON`` file can be loaded into PipeCraft2 to **automatically configure your next runs exactly the same way**.

___________________________________________________

Start the workflow
~~~~~~~~~~~~~~~~~~

Press ``START`` on the left ribbon **to start the analyses**.

.. admonition:: when running the module for the first time ...
  
  ... a docker image will be first pulled to start the process. 

  For example: |pulling_image|


When you need to STOP the workflow, press ``STOP`` button |stop_workflow|


.. admonition:: When the workflow has completed ...

  ... a message window will be displayed.

  |workflow_finished|

___________________________________________________

Examine the outputs
~~~~~~~~~~~~~~~~~~~

Several process-specific output folders are generated |output_icon|

