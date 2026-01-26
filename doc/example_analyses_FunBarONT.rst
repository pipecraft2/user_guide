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
        PipeCraft manual. FunBarONT workflow tutorial

|

FunBarONT pipeline, ITS |PipeCraft2_logo|
-------------------------------------------

This example data analyses follows **FunBarONT** workflow as implemented in PipeCraft2's pre-compiled pipelines panel. 

FunBarONT is a specialized pipeline for processing **Oxford Nanopore Technologies (ONT) fungal barcoding data**, 
specifically targeting the **ITS2 rRNA gene region**. This pipeline is optimized for long-read sequencing data and incorporates 
quality filtering, demultiplexing, sequence polishing, and taxonomic assignment to generate high-confidence fungal identifications.

| Example data set details: **Oxford Nanopore long-read ITS2 amplicon sequences** (see nano_test_data folder).
| This is a sample dataset for **fungal identification and barcoding** using ITS2 amplicons.

____________________________________________________

Starting point 
~~~~~~~~~~~~~~

This example dataset consists of **ITS2 rRNA gene amplicon sequences**; targeting fungi:

- **single-end** Oxford Nanopore sequencing data;
- **demultiplexed** set (per-sample fastq files, typically demultiplexed using Guppy, MinKNOW, or similar);
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
| 
| **To select input data**, press ``SELECT WORKDIR``
| and specify
| ``sequence files extension`` as **\*.fastq** or **\*.fastq.gz**;  
| ``sequencing read types`` as **single-end**.

____________________________________________________

Workflow overview
~~~~~~~~~~~~~~~~~

The FunBarONT pipeline consists of several key processing steps designed to handle the characteristics of Oxford Nanopore long-read sequencing data:

1. **Input validation and quality assessment** - Evaluates read quality and length distribution
2. **Length filtering** - Removes reads that do not meet minimum length requirements for reliable ITS2 amplicon detection
3. **Adapter/primer removal** - Optional trimming of known sequences from the terminal ends of reads
4. **Quality filtering** - Applies quality thresholds to remove low-quality reads
5. **Sequence polishing** - Corrects sequencing errors common in long-read data
6. **Clustering/OTU generation** - Groups similar sequences into operational taxonomic units (OTUs) or ASVs
7. **Chimera filtering** - Identifies and removes potential chimeric sequences
8. **ITS2 extraction** - Extracts and isolates the ITS2 subregion for improved taxonomic accuracy
9. **Taxonomy assignment** - Assigns taxonomic classification based on reference databases

.. note::

  The exact processing steps and their order may be customized in the PipeCraft2 interface. 
  The sequence described here represents a typical recommended workflow for fungal ITS2 barcoding with ONT data.

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

- **Read length**: Ensure that the expected amplicon length (including ITS2 and flanking regions) matches your read lengths. 
  The default minimum length filtering is typically set to accommodate full-length ITS2 amplicons (~500-700 bp for fungi).

- **Mixed samples**: If your data contains mixed fungal species or environmental samples, the clustering steps will group similar sequences together.

____________________________________________________

Quality control
~~~~~~~~~~~~~~~

Before processing, it is advisable to assess the quality of your sequencing data. 
You can use the ``QualityCheck`` panel in PipeCraft2 to visualize read quality distribution, length distribution, and other metrics.

Key quality metrics to consider:

- **Mean quality score**: Typically > Q10 is acceptable for ONT basecalling data
- **Read length distribution**: Should center around your expected amplicon length
- **Quality score trends**: Quality often decreases towards the end of reads in ONT data

____________________________________________________

Length filtering
~~~~~~~~~~~~~~~~

The FunBarONT pipeline can filter reads based on minimum and/or maximum length thresholds. 

This step is particularly important for ONT data because:

- Very short reads may represent sequencing artifacts or incomplete amplicons
- Very long reads may indicate concatenated or chimeric sequences
- Length filtering helps focus on full-length, high-quality amplicons

**Recommended settings for ITS2 fungal barcoding:**

- **Minimum length**: 400-500 bp (adjust based on your specific primers and target amplicon size)
- **Maximum length**: 1500-2000 bp (adjust based on expected maximum read length)

____________________________________________________

Cut primers
~~~~~~~~~~~

If your input sequences still contain primers or adapter sequences, this step removes them.

This is important because:

- Primer sequences can bias downstream analyses and clustering
- Removing primers focuses analysis on the variable metabarcode region
- Some taxonomy assignment databases expect primer-free sequences

**For ITS2 fungal amplicons:**

Specify the forward and reverse primer sequences used in your PCR amplification. 
For example, common ITS2 primers targeting fungi include:

- **Forward**: GCATCGATGAAGAACGCAGC (fITS7)
- **Reverse**: TCCTCCGCTTATTGATATGC (ITS4)

Adjust the ``min overlap`` and ``mismatches`` parameters based on:
- Primer length
- Expected sequence quality
- Tolerance for primer sequence variations

.. admonition:: when working with your own ITS2 data ...

  ... if you are using **ITSx** for ITS2 subregion extraction, you may optionally skip the primer cutting step, 
  as ITSx will extract the ITS2 region and remove flanking sequences (18S and 5.8S rRNA genes) automatically.

____________________________________________________

Quality filtering
~~~~~~~~~~~~~~~~~

Quality filtering removes low-quality reads that do not meet specified error thresholds. 
This step is critical for Oxford Nanopore data, which typically has higher error rates than short-read sequencing.

Two main approaches are commonly used:

1. **Percentage-based filtering**: Keep only reads with a certain percentage of bases above a quality threshold
2. **Expected error filtering**: Remove reads exceeding a maximum number of expected errors

For ONT data processing:

- Use more lenient quality thresholds compared to Illumina data (due to the inherent characteristics of long-read sequencing)
- Consider the read length when setting quality thresholds; longer reads with uniform quality are acceptable
- The default settings in FunBarONT are pre-optimized for typical long-read ITS2 data

+-----------------------+-------------------------------------------------------+
| Output directory |output_icon|          ``qualFiltered_out``                  |
+=======================+=======================================================+
| \*.fastq              | quality filtered sequences per sample in FASTQ format |
+-----------------------+-------------------------------------------------------+
| seq_count_summary.txt | summary of sequence counts per sample                 |
+-----------------------+-------------------------------------------------------+

____________________________________________________

Sequence polishing
~~~~~~~~~~~~~~~~~~

Oxford Nanopore long reads often contain random errors that can be corrected using sequence polishing algorithms. 
The FunBarONT pipeline may include optional polishing steps to improve sequence quality.

Common polishing approaches include:

- **Consensus-based polishing**: Multiple reads are aligned and a consensus sequence is generated
- **Machine learning-based polishing**: Uses trained models to correct systematic errors
- **Reference-based polishing**: Aligns reads against a reference database and corrects errors

.. note::

  Polishing is optional and may increase processing time. It is particularly beneficial when:
  
  - Working with low-quality sequencing runs
  - Processing mixed environmental samples with variable coverage per species
  - Requiring high-confidence sequences for downstream analyses

____________________________________________________

Length filtering (post-processing)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After quality filtering and polishing, an additional length filtering step may remove sequences that have become 
too short during the quality trimming process.

This ensures that downstream analyses (clustering, taxonomy assignment) only use sequences of appropriate length.

____________________________________________________

Clustering
~~~~~~~~~~

The clustering step groups similar sequences into operational taxonomic units (OTUs). 
The FunBarONT pipeline can use several clustering methods:

- **vsearch clustering**: Fast, deterministic clustering at a specified similarity threshold
- **SWARM clustering**: Abundance-based clustering that does not require pre-defined similarity thresholds
- **ASV-based approaches**: Each unique sequence is treated as an amplicon sequence variant (ASV)

**Recommended settings for ITS2 fungal barcoding:**

- **Similarity threshold**: 97% is commonly used for ITS2 fungal OTU clustering (adjust based on your taxonomic resolution needs)
- **Strand specification**: Set to "both" to cluster sequences regardless of orientation
- **Abundance weighting**: Consider enabling if using SWARM clustering to account for sequence abundance in the sample

**Why clustering for ONT data:**

- Reduces impacts of sequencing errors by grouping similar sequences
- Produces more stable OTU tables compared to treating every sequence as unique (ASV approach)
- Improves computational efficiency for downstream analyses

+-------------------------------------+-------------------------------------------+
| Output directory |output_icon|  ``clustering_out``                              |
+=====================================+===========================================+
| OTU_table.txt                       | OTU-by-sample abundance table             |
+-------------------------------------+-------------------------------------------+
| OTUs.fasta                          | representative sequences per OTU in FASTA |
+-------------------------------------+-------------------------------------------+
| OTUs.uc                             | clustering results mapping file           |
+-------------------------------------+-------------------------------------------+

____________________________________________________

Chimera filtering
~~~~~~~~~~~~~~~~~

Chimera detection and filtering is a crucial quality control step for amplicon sequencing data. 
Chimeras are hybrid sequences formed from two or more biological sequences, typically arising during PCR amplification.

The FunBarONT pipeline offers multiple chimera detection methods:

- **De novo chimera detection** (uchime_denovo): Identifies chimeras based on sequence patterns without requiring a reference database
- **Reference-based chimera detection** (uchime_ref): Compares against a reference database to identify chimeric sequences
- **Combined approach**: De novo filtering followed by reference-based filtering for maximum sensitivity

**For fungal ITS2 data:**

We recommend using the **de novo approach** as the primary method, with optional reference-based filtering using a fungal ITS database such as:

- UNITE database (https://unite.ut.ee/)
- RDP Fungal ITS database

.. important:: 

  Ensure that primers have been removed from your amplicons before chimera filtering; 
  otherwise, legitimate sequences may be incorrectly flagged as chimeras. 

.. admonition:: when working with ITS2 fungal data ...

  ... download the appropriate reference database for your target organisms. 
  The UNITE database is widely used for fungal ITS analysis and is available in multiple formats 
  including UCHIME/USEARCH format suitable for reference-based chimera filtering.

+---------------------------+---------------------------------------------+
| Output directory  |output_icon|  ``chimera_Filtered_out``               |
+===========================+=============================================+
| OTU_table.txt or \*.fasta | chimera-filtered OTU table and/or sequences |
+---------------------------+---------------------------------------------+
| seq_count_summary.txt     | summary of sequence counts per sample       |
+---------------------------+---------------------------------------------+
| ``chimeras``/\*.fasta     | discarded sequences identified as chimeras  |
+---------------------------+---------------------------------------------+

____________________________________________________

Extract ITS2 subregion
~~~~~~~~~~~~~~~~~~~~~~

The ITS2 subregion extraction step uses ITSx software to identify and extract the ITS2 rRNA gene region from your sequences. 
This is particularly valuable for fungal identification because:

- **Improves taxonomic accuracy**: The ITS2 region is more conserved at certain taxonomic levels, improving classification
- **Removes non-functional regions**: Eliminates 18S and 5.8S rRNA genes that may bias clustering and taxonomy assignment
- **Standardizes sequence length**: Produces more uniform sequence lengths for improved clustering

**Setting for FunBarONT ITS2 extraction:**

- **Region for clustering**: Set to ``ITS2`` (the specific region of interest)
- **Organisms**: Select ``fungi`` to limit the search to fungal ITS regions
- **Complement**: Enable to also search the reverse complement strand (important if sequences have mixed orientation)

**Important considerations:**

- ITSx works best on full or near-full length amplicons
- Very short or partial sequences may not be reliably detected
- The tool may produce both "full" and "partial" ITS2 sequences; review the output to decide which to use

.. note::

  If you have already performed primer cutting in earlier steps and your sequences are in consistent orientation,
  you may disable the ``complement`` search to reduce processing time.

+-----------------------------------------+-------------------------------------------------------------+
| Output directory |output_icon| ``ITSx_out``                                                           |
+=========================================+=============================================================+
| ``ITS2``/\*.fasta                       | ITS2 sequences (without flanking regions) per sample        |
+-----------------------------------------+-------------------------------------------------------------+
| ``ITS2``/``full_and_partial``/\*.fasta  | full and partial ITS2 sequences per sample                  |
+-----------------------------------------+-------------------------------------------------------------+
| seq_count_summary.txt                   | summary of sequence counts per sample                       |
+-----------------------------------------+-------------------------------------------------------------+

____________________________________________________

Taxonomy assignment
~~~~~~~~~~~~~~~~~~~

Taxonomy assignment annotates your OTUs with taxonomic classifications based on comparison against reference databases. 
The FunBarONT pipeline supports multiple taxonomy assignment methods:

1. **BLAST** - Basic Local Alignment Search Tool (uses top BLAST hits)
2. **SINTAX** - Fast, simple classifier (requires formatted database)
3. **UTAX** - Naive Bayes classifier (vsearch UTAX)
4. **Custom methods** - Can be integrated based on your specific needs

**For fungal ITS2 barcoding, recommended approaches:**

- **BLAST with UNITE database** - Most widely used combination, allows easy interpretation of hits
- **SINTAX with UNITE database** - Faster than BLAST with reasonable accuracy

**Setting up taxonomy assignment for ITS2 fungi:**

1. Download a suitable reference database (e.g., UNITE fungal ITS database)
2. Specify the database location in the ``taxonomy assignment`` panel
3. Adjust ``minimum identity`` threshold (typically 85-97% for ITS2)
4. If using BLAST, select whether you want ``1st best hit`` or ``top 10 hits``

**Interpreting taxonomy output:**

The taxonomy assignment process produces classifications at multiple taxonomic ranks:

- Kingdom
- Phylum
- Class
- Order
- Family
- Genus
- Species

Not all ranks may have reliable classifications for every OTU; "Unclassified" entries indicate where confidence is insufficient.

.. admonition:: when working with your fungal barcoding data ...

  ... note that ITS2 provides excellent discrimination at family and genus levels for fungi, 
  but species-level identification may require additional information (such as morphology) 
  for accurate determination. Some OTUs may represent:
  
  - **Species complexes** - Closely related species that are difficult to distinguish by ITS2 alone
  - **Rare or undescribed species** - Not represented in the reference database
  - **Heterogeneous OTUs** - Multiple species grouped due to sequence similarity

+-------------------------------------+--------------------------------+
| Output directory     |output_icon|    ``taxonomy_out``               |
+=====================================+================================+
| BLAST_1st_best_hit.txt              | BLAST 1st hit per OTU/ASV      |
+-------------------------------------+--------------------------------+
| BLAST_10_best_hits.txt              | First 10 BLAST hits per OTU    |
+-------------------------------------+--------------------------------+
| taxonomy.txt or SINTAX_taxonomy.txt | Taxonomic classification table |
+-------------------------------------+--------------------------------+

____________________________________________________

Data interpretation and post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the FunBarONT workflow is complete, you have several files for further analysis:

**Primary output files:**

- **OTU_table.txt** - Abundance table (OTUs × samples)
- **OTUs.fasta** - Representative sequences for each OTU
- **ITS2.fasta** - ITS2-extracted sequences (if ITSx step was performed)
- **taxonomy.txt** - Taxonomic classification for each OTU

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
- Confidence scores or assignment method details

**Next steps for your analysis:**

1. Load the OTU table into R or Python for statistical analysis
2. Perform rarefaction analysis to assess sampling completeness
3. Conduct community composition analysis (e.g., bar plots, PCA)
4. Calculate alpha and beta diversity metrics
5. Perform differential abundance testing if comparing groups
6. Generate phylogenetic trees for evolutionary insights

____________________________________________________

References and further reading
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **ITSx**: `Bengtsson-Palme et al. 2013 <https://microbiology.se/software/itsx/>`_ - 
  Identifies and extracts ITS regions from rRNA gene sequences

- **UNITE database**: `Koljalg et al. 2013 <https://unite.ut.ee/>`_ - 
  Fungal ITS sequence reference database for taxonomic classification

- **Oxford Nanopore Technologies**: Comprehensive sequencing platform documentation and basecalling guidelines

- **vsearch**: `Rognes et al. 2016 <https://github.com/torognes/vsearch>`_ - 
  Fast and versatile open-source tool for metagenomics

- **BLAST+**: `Camacho et al. 2009 <https://blast.ncbi.nlm.nih.gov/>`_ - 
  Suite of programs for sequence similarity searches

____________________________________________________

Citation
~~~~~~~~

If you use the FunBarONT pipeline in your research, please cite:

- PipeCraft2: `Anslan et al. 2020+ <https://github.com/pipecraft2/>`_
- Oxford Nanopore long-read sequencing technology
- ITSx software used for ITS2 extraction
- Reference databases (UNITE) used for taxonomy assignment

____________________________________________________

|
