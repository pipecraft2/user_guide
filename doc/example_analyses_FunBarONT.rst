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

.. |funbaront_workflow| image:: _static/funbaront_workflow.png
  :width: 800
  :alt: FunBarONT workflow overview

.. meta::
    :description lang=en:
        PipeCraft manual. FunBarONT workflow tutorial

|

FunBarONT pipeline, ITS |PipeCraft2_logo|
-------------------------------------------

This example data analyses follows **FunBarONT** workflow as implemented in PipeCraft2's pre-compiled pipelines panel. 

FunBarONT is a specialized pipeline for processing **Oxford Nanopore Technologies (ONT) fungal barcoding data**, 
specifically targeting the **ITS rRNA gene region**. This pipeline is optimized for long-read sequencing data and incorporates 
quality filtering, demultiplexing, sequence polishing, and taxonomic assignment to generate high-confidence fungal identifications.

| `Download example data set here: <https://raw.githubusercontent.com/pipecraft2/user_guide/master/data/example_data_FunBarONT.zip>`_ and **unzip** it.
| This is a sample dataset for **fungal identification and barcoding** using ITS amplicons.

____________________________________________________

Starting point 
~~~~~~~~~~~~~~

This example dataset consists of **ITS rRNA gene amplicon sequences**; targeting fungi:

- **single-end** Oxford Nanopore sequencing data;
- **demultiplexed** set (per-sample fastq files, typically demultiplexed using cutadapt, MinKNOW, or similar);
- barcodes and adapters have already been **removed**;
- sequences are generated using **long-read sequencing technology** (read lengths typically 1-10+ kb).

.. admonition:: when working with your own ONT data ...

  ... then please ensure that the fastq files are properly demultiplexed and contain sample identifiers in the file names.
  For mixed barcode runs that have not been demultiplexed, preprocessing may be required to split reads by barcode before 
  using the FunBarONT pipeline.

  | *Example demultiplexed file naming:*
  | *sample1.fastq*
  | *sample2.fastq*
  | *sample3.fastq*

____________________________________________________

| **To select FunBarONT pipeline**, press
| ``SELECT PIPELINE`` --> ``FunBarONT``.

|funbaront_workflow|
| 
| **To select input data**, press ``SELECT WORKDIR``
| and specify
| ``sequence files extension`` as **\*.fastq** or **\*.fastq.gz**;  
| ``sequencing read types`` as **single-end**.

____________________________________________________

Workflow overview
~~~~~~~~~~~~~~~~~

The FunBarONT pipeline consists of the following processing steps designed to handle the characteristics of Oxford Nanopore long-read sequencing data:

1. **Quality Control (NanoPlot)** - Generates quality reports and statistics for each sample
2. **Quality Filtering (chopper)** - Filters reads based on quality scores and length thresholds
3. **Clustering (VSEARCH)** - Groups similar sequences into clusters/OTUs
4. **Sequence Polishing (racon + medaka)** - Corrects sequencing errors to generate high-accuracy consensus sequences
5. **ITS Extraction (ITSx)** - Extracts the ITS region from fungal sequences (optional)
6. **Taxonomy Assignment (BLAST)** - Assigns taxonomic classification using BLAST against a reference database

.. note::

  The FunBarONT pipeline is specifically designed for Oxford Nanopore fungal barcoding data. 
  All steps run automatically in sequence once the workflow is started.

____________________________________________________

Input data preparation
~~~~~~~~~~~~~~~~~~~~~~

Before starting the FunBarONT workflow in PipeCraft2, ensure your input data is properly prepared:

**File naming convention:**

All fastq files should follow a consistent naming pattern with the sample identifier at the beginning:

.. code-block::
   :caption: Recommended file naming structure

    sample_name.fastq
    sample_name.fastq.gz


**Directory structure:**

.. code-block::
   :caption: Recommended directory structure for FunBarONT

    my_fungal_barcoding/   
    └── sequences/         # SELECT THIS FOLDER AS WORKING DIRECTORY
        ├── sample1.fastq
        ├── sample2.fastq
        ├── sample3.fastq
        └── ...

**Data quality considerations:**

- **Read quality**: Oxford Nanopore reads can contain sequencing errors, particularly towards the ends of reads. 
  The pipeline includes quality filtering steps to handle this.

- **Read length**: Ensure that the expected amplicon length (including ITS and flanking regions) matches your read lengths. 
  The default minimum length filtering is typically set to accommodate full-length ITS amplicons (~500-700 bp for fungi).

- **Mixed samples**: If your data contains mixed fungal species or environmental samples, the clustering steps will group similar sequences together.

____________________________________________________

Quality Control (NanoPlot)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The FunBarONT pipeline uses **NanoPlot** to assess the quality of your Oxford Nanopore sequencing data. 
This step generates comprehensive quality reports and statistics for each sample.

NanoPlot produces:

- **Quality distribution plots** - Visualize the distribution of read quality scores
- **Read length distribution** - Shows the length distribution of your sequencing reads
- **NanoStats.txt** - Summary statistics including read counts, mean quality, and length metrics
- **NanoPlot-report.html** - Interactive HTML report with all quality metrics

Key quality metrics to consider:

- **Mean quality score**: Typically > Q10 is acceptable for ONT basecalling data
- **Read length distribution**: Should center around your expected amplicon length
- **Total reads**: Verify adequate sequencing depth per sample

+-----------------------------------------------+---------------------------------------------------+
| Output directory |output_icon|                | ``01_quality_reports``                            |
+===============================================+===================================================+
| ``<sample>_NanoPlot_results``/                | NanoPlot results folder per sample                |
+-----------------------------------------------+---------------------------------------------------+
| NanoPlot-report.html                          | interactive HTML quality report                   |
+-----------------------------------------------+---------------------------------------------------+
| NanoStats.txt                                 | summary statistics (read counts, quality, length) |
+-----------------------------------------------+---------------------------------------------------+

____________________________________________________

Quality Filtering (chopper)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quality filtering uses **chopper** to remove low-quality reads that do not meet specified thresholds. 
This step is critical for Oxford Nanopore data, which can have variable error rates.

**Configurable parameters:**

- **chopper_quality** (default: 10) - Minimum read quality score. Reads below this threshold are discarded.
- **chopper_min_read_length** (default: 150 bp) - Minimum read length. Shorter reads are removed.
- **chopper_max_read_length** (default: 1000 bp) - Maximum read length. Longer reads are removed.

**For ITS fungal barcoding:**

Adjust the length thresholds based on your expected amplicon size. The ITS region in fungi typically ranges from 400-800 bp, 
so default settings should work well for most applications.

+--------------------------------------------+-------------------------------------------------------+
| Output directory |output_icon|             | ``02_filtered_sequences``                             |
+============================================+=======================================================+
| \*.chopper.fasta.gz                        | quality filtered sequences per sample in FASTA format |
+--------------------------------------------+-------------------------------------------------------+

____________________________________________________

Clustering (VSEARCH)
~~~~~~~~~~~~~~~~~~~~

The clustering step uses **VSEARCH** to group similar sequences into clusters. 
This reduces the impact of sequencing errors and produces representative sequences for downstream analysis.

**Configurable parameters:**

- **vsearch_cluster_id** (default: 0.95) - Clustering identity threshold (0-1). Sequences with similarity above this threshold will be clustered together.
- **vsearch_cluster_strand** (default: "both") - Check both strands or plus strand only during clustering.

**Why clustering for ONT data:**

- Reduces impacts of sequencing errors by grouping similar sequences
- Produces representative centroid sequences for each cluster
- Improves computational efficiency for downstream polishing and taxonomy assignment

+---------------------------------------+-----------------------------------------------+
| Output directory |output_icon|        | ``03_clusters``                               |
+=======================================+===============================================+
| \*.centroids.fasta.gz                 | representative centroid sequences per sample  |
+---------------------------------------+-----------------------------------------------+

____________________________________________________

Sequence Polishing (racon + medaka)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Oxford Nanopore long reads often contain random errors that are corrected using a two-step polishing process:

1. **Racon** - First-pass polishing using multiple sequence alignment
2. **Medaka** - Neural network-based polishing for high-accuracy consensus sequences

**Configurable parameters:**

- **medaka_model** (default: r1041_e82_400bps_hac_variant_v4.3.0) - Select the medaka model based on your flowcell, kit, and basecaller.
- **racon_quality_threshold** (default: 20) - Minimum average base quality for windows used by Racon.
- **racon_window_length** (default: 100) - Window length used by Racon for polishing.

.. note::

  Select the appropriate **medaka model** based on your sequencing setup:
  
  - **r1041_e82** models are for R10.4.1 flowcells with E8.2 chemistry
  - **r941** models are for R9.4.1 flowcells
  - Choose **hac** (high accuracy) or **sup** (super accuracy) based on your basecalling model
  - The model affects consensus accuracy, so matching your setup is important

+---------------------------------------+-----------------------------------------------+
| Output directory |output_icon|        | ``04_polished_sequences``                     |
+=======================================+===============================================+
| \*.racon.fasta                        | Racon-polished sequences per sample           |
+---------------------------------------+-----------------------------------------------+
| \*.medaka.consensus.fasta             | Medaka-polished consensus sequences per sample|
+---------------------------------------+-----------------------------------------------+

____________________________________________________

ITS Extraction (ITSx)
~~~~~~~~~~~~~~~~~~~~~

The ITS extraction step uses **ITSx** software to identify and extract the ITS rRNA gene region from your sequences. 
This is particularly valuable for fungal identification because:

- **Improves taxonomic accuracy**: The ITS region is the standard barcode for fungi
- **Removes flanking regions**: Eliminates 18S and 5.8S rRNA genes that may bias taxonomy assignment
- **Standardizes sequences**: Produces comparable sequences for database matching

**Pipeline option:**

- **use_itsx** (default: true) - Set to false if you want to skip ITS extraction (useful for non-ITS sequences)

**Important considerations:**

- ITSx works best on full or near-full length amplicons
- Sequences without detectable ITS regions will produce empty output files
- The tool extracts both ITS1 and ITS2 regions when present

+---------------------------------------+---------------------------------------------------+
| Output directory |output_icon|        | ``05_its_extracted``                              |
+=======================================+===================================================+
| \*.its.fasta                          | extracted ITS sequences per sample                |
+---------------------------------------+---------------------------------------------------+

____________________________________________________

Taxonomy Assignment (BLAST)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Taxonomy assignment uses **BLAST** to compare your sequences against a reference database and assign taxonomic classifications.

**Configurable parameters:**

- **database_file** (required) - Reference database file in FASTA format. The pipeline will automatically create a BLAST database from this file.
- **run_id** (default: "funbaront_run") - Unique identifier for this analysis run. Used for naming output directories and files.
- **strands** (default: "both") - Query strand to search against database. "both" includes reverse complement.
- **e_value** (default: 10) - E-value threshold. Lower values indicate more significant matches.
- **word_size** (default: 11) - Initial word size for BLAST alignment.

**Additional pipeline options:**

- **output_all_polished_seqs** (default: false) - Output all polished sequences even those without database hits (useful for non-ITS sequences).
- **rel_abu_threshold** (default: 10) - Output only clusters with barcode-wise relative abundance above this percentage (0-100).

**For fungal ITS barcoding:**

Use a fungal ITS reference database such as:

- **UNITE database** (https://unite.ut.ee/) - The standard reference database for fungal ITS sequences

**Interpreting taxonomy output:**

The BLAST results include taxonomic assignments at multiple ranks when available in the reference database.
Not all sequences may have reliable classifications; sequences without database hits will have empty results.

+---------------------------------------+-------------------------------------------+
| Output directory |output_icon|        | ``06_blast_results``                      |
+=======================================+===========================================+
| \*.blast.tsv                          | BLAST results in tabular format per sample|
+---------------------------------------+-------------------------------------------+

____________________________________________________

Final Results and JSON Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FunBarONT pipeline produces comprehensive final results in multiple formats:

+---------------------------------------+-----------------------------------------------------------+
| Output directory |output_icon|        | ``07_json_results``                                       |
+=======================================+===========================================================+
| \*.results.json                       | JSON formatted results per sample with all analysis data  |
+---------------------------------------+-----------------------------------------------------------+

**Main results file:**

+---------------------------------------+-----------------------------------------------------------+
| Output file |output_icon|             | Root output directory                                     |
+=======================================+===========================================================+
| funbaront_run.results.xlsx            | Excel spreadsheet with all results (taxonomy, quality)    |
+---------------------------------------+-----------------------------------------------------------+
| README.md                             | Summary of the pipeline run with parameters and citations |
+---------------------------------------+-----------------------------------------------------------+

____________________________________________________

Data interpretation and post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the FunBarONT workflow is complete, you have several files for further analysis:

**Primary output files:**

- **funbaront_run.results.xlsx** - Excel spreadsheet with comprehensive results including sequence info, taxonomy, and quality metrics
- **\*.blast.tsv** - BLAST results per sample in tabular format
- **\*.its.fasta** - Extracted ITS sequences per sample
- **\*.medaka.consensus.fasta** - Polished consensus sequences per sample
- **\*.results.json** - JSON formatted results for programmatic access

**Output directory structure:**

.. code-block::
   :caption: FunBarONT output directory structure

    <run_id>_results/
    ├── 01_quality_reports/          # NanoPlot quality reports per sample
    │   └── <sample>_NanoPlot_results/
    ├── 02_filtered_sequences/       # Chopper-filtered sequences
    ├── 03_clusters/                 # VSEARCH clustering centroids
    ├── 04_polished_sequences/       # Racon and Medaka polished sequences
    ├── 05_its_extracted/            # ITSx extracted ITS sequences
    ├── 06_blast_results/            # BLAST taxonomy results
    ├── 07_json_results/             # JSON formatted results
    ├── funbaront_run.results.xlsx   # Main results spreadsheet
    └── README.md                    # Run summary and citations

**Recommended post-processing steps:**

1. **Filter for target kingdom/phylum** - Remove non-fungal OTUs if present
2. **Examine unclassified OTUs** - BLAST search or phylogenetic analysis
3. **Summary statistics** - Calculate diversity metrics, rarefaction curves
4. **Downstream analyses** - Community composition, differential abundance, etc.

**Quality checks before downstream analyses:**

- Verify that the majority of OTUs are classified as fungi
- Check for contamination (e.g., unusual taxa)
- Review the rarefaction curves to assess sampling depth
- Examine OTU abundance distributions

____________________________________________________

Save workflow configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have configured all parameters for your FunBarONT analysis, you can save the configuration file by pressing ``save workflow`` button on the right ribbon.

|save|

The configuration file will be saved as a JSON file (e.g., ``pipecraft2_last_run_configuration.json``) in your working directory. 
This file stores all your parameter choices and can be reloaded into PipeCraft2 to reproduce the exact same analysis in future runs.

.. admonition:: automatic configuration backup

  If you forget to manually save the configuration, PipeCraft2 will automatically generate 
  ``pipecraft2_last_run_configuration.json`` upon starting the workflow. 
  This serves as a backup of your analysis parameters for this working directory.

____________________________________________________

Start the workflow
~~~~~~~~~~~~~~~~~~

Once all parameters have been configured, press ``START`` on the left ribbon to begin the FunBarONT analysis.

|workflow_finished|

The workflow will proceed through each step in sequence, with progress displayed in the PipeCraft2 interface.

.. admonition:: first-time execution notes

  ... when running the FunBarONT pipeline for the first time, Docker will automatically pull the required container image. 
  This may take several minutes depending on your internet connection and the image size.

  |pulling_image|

  Subsequent runs will use the cached image and will start more quickly.

**Monitoring progress:**

- Each completed step will display a checkmark
- Error messages will appear if any step fails
- Processing time depends on the size of your dataset and computational resources available

____________________________________________________

Troubleshooting
~~~~~~~~~~~~~~~

**Common issues and solutions:**

.. admonition:: Docker image pull fails

  **Issue**: Error message about pulling the FunBarONT Docker image
  
  **Solution**: Check your internet connection and ensure sufficient disk space for the image. 
  Ensure Docker daemon is running properly.

.. admonition:: Insufficient reads in output

  **Issue**: Output OTU table contains very few sequences after processing
  
  **Solution**: 
  
  - Check quality filtering thresholds; they may be too stringent
  - Verify primer sequences are correct
  - Ensure input data quality is adequate
  - Review quality control metrics for the sequencing run

.. admonition:: High proportion of unclassified OTUs

  **Issue**: Many OTUs receive "Unclassified" taxonomy assignments
  
  **Solution**:
  
  - Lower the minimum identity threshold in taxonomy assignment
  - Verify the reference database is appropriate for your samples
  - Consider BLAST search of unclassified OTUs against public databases (NCBI-nr)

.. admonition:: Unexpected taxa in results

  **Issue**: Non-fungal taxa appearing in taxonomy assignments
  
  **Solution**:
  
  - Filter output based on Kingdom classification (keep only Fungi)
  - Review primer specificity for potential off-target amplification
  - Consider whether environmental contamination is likely

____________________________________________________

Example output analysis
~~~~~~~~~~~~~~~~~~~~~~~~

After the workflow completes, you can examine the output files to understand your fungal community composition.

**Example taxonomy summary visualization:**

The taxonomy output typically includes a comma-separated or tab-delimited table with columns for:

- OTU ID or sequence identifier
- Full sequence
- Kingdom, Phylum, Class, Order, Family, Genus, Species
- Confidence scores or assignment method detailsSetting

Tooltip

database_file

select a database file in fasta format. Fasta format will be
automatically converted to BLAST database
fasta_file

select a fasta file to be used as a query for BLAST search

task

BLAST search settings according to blastn or megablast

strands

query strand to search against database. Both = search also reverse
complement
e_value


a parameter that describes the number of hits one can expect to see
by chance when searching a database of a particular size. The lower
the e-value the more ‘significant’ the match is
word_size

the size of the initial word that must be matched between the
database and the query sequence
reward

reward for a match

penalty

penalty for a mismatch

gap_open

cost to open a gap

gap_extend

cost to extend a g