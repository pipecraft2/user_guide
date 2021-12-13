.. image:: _static/PipeCraft2_icon_v2.png
  :width: 80
  :alt: logo

.. |main_interface| image:: _static/main_interface.png
  :width: 2000
  :alt: Alternative text

.. |asv_main| image:: _static/asv_main.png
  :width: 1500
  :alt: Alternative text

.. |otu_main| image:: _static/otu_main.png
  :width: 1500
  :alt: Alternative text

.. |console| image:: _static/console.png
  :width: 1500
  :alt: Alternative text

==========
User guide
==========

The interface
==============

The startup panel:

|main_interface|

____________________________________________________

FULL PIPELINE PANELS
====================

.. _asvpipe:

ASVs workflow panel (with `DADA2 <https://benjjneb.github.io/dada2/index.html>`_)
----------------------------------------------------------------------------------

.. note::
  Current (v2.0.1) ASVs workflow works only with **PAIRED-END** reads!

|asv_main|

This automated workflow is based on DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html 
 | Note that ``demultiplexing``, ``reorient`` and ``remove primers`` steps are optional and do not represent parts from DADA2 tutorial. Nevertheless, it is advisable to :ref:`reorient <reorinet>` your reads (to 5'-3') and :ref:`remove primers <remove_primers>` before proceeding with ASV generation with DADA2.

| `DADA2 manual is here <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_
 
.. _dada2_defaults:

**Default options:**

================================================= =========================
Analyses step                                     Default setting
================================================= =========================
:ref:`DEMULTIPLEX <demux>` (optional)              --
:ref:`REORIENT <reorinet>` (optional)              --
:ref:`REMOVE PRIMERS <remove_primers>` (optional)  --
:ref:`QUALITY FILTERING <dada2_qual_filt>`         | ``read R1`` = _R1
                                                   | ``read R2`` = _R2
                                                   | ``samp ID`` = _
                                                   | ``maxEE`` = 1
                                                   | ``maxN`` = 0
                                                   | ``minLen`` = 32
                                                   | ``truncQ`` = 2
                                                   | ``truncLen`` = 0
                                                   | ``maxLen`` = 600
                                                   | ``minQ`` = 2
:ref:`DENOISE <dada2_denoise>`                     | ``pool`` = FALSE
                                                   | ``selfConsist`` = FASLE
                                                   | ``qualityType`` = Auto
:ref:`MERGE PAIRED-END READS <dada2_merge_pairs>`  | ``minOverlap`` = 12
                                                   | ``maxMismatch`` = 0
                                                   | ``returnRejects`` = FALSE
:ref:`CHIMERA FILTERING <dada2_chimeras>`          | ``method`` = consensus
:ref:`ASSGIN TAXONOMY <dada2_taxonomy>`            | ``minBoot`` = 50
                                                   | ``tryRC`` = FALSE
                                                   | ``refFasta`` = select a database
================================================= =========================

____________________________________________________

.. _dada2_qual_filt:

QUALITY FILTERING [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `filterAndTrim <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to 
perform quality filtering on input fastq files based on user selected criteria. Outputs filtered fastq files into ``qualFiltered_out.dada2`` directory.

==================== ============
Setting              Tooltip
==================== ============
``read R1``          | identifier string for R1 reads. Default = _R1, 
                     | which means that all R1 reads in a directory may be identified via latter string
``read R2``          | identifier string for R2 reads. Default = _R2
``samp ID``          | identifyer string that separates the sample name from redundant charachters 
                     | (e.g. file name = sampl84_S73_L001_R1_001.fastq, then underscore '_' would be 
                     | the 'identifier string' (sample name = sampl84))
``maxEE``            | discard sequences with more than the specified number of expected errors
``maxN``             | discard sequences with more than the specified number of N’s (ambiguous bases)
``minLen``           | remove reads with length less than minLen. minLen is enforced after all other 
                     | trimming and truncation
``truncQ``           | truncate reads at the first instance of a quality score less than or equal to truncQ
``truncLen``         | truncate reads after truncLen bases (applies to R1 reads when working with paired-end data). 
                     | Reads shorter than this are discarded. Explore quality profiles (with QualityCheck module) 
                     | see whether poor quality ends needs to truncated
``truncLen_R2``      | truncate R2 reads after truncLen bases. 
                     | Reads shorter than this are discarded. Explore quality profiles 
                     | (with QualityCheck module) see whether poor quality ends needs to truncated
``maxLen``           | remove reads with length greater than maxLen. maxLen is enforced on the raw reads. 
                     | In dada2, the default = Inf, but here set as 9999
``minQ``             | after truncation, reads contain a quality score below minQ will be discarded
==================== ============

see :ref:`default settings <dada2_defaults>`

____________________________________________________

.. _dada2_denoise:

DENOISING [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `dada <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to remove sequencing errors.
Outputs filtered fasta files into ``denoised_assembled.dada2`` directory.

==================== ============
Setting              Tooltip
==================== ============
``pool``             | if TRUE, the algorithm will pool together all samples prior to sample inference. 
                     | Pooling improves the detection of rare variants, but is computationally more expensive. 
                     | If pool = 'pseudo', the algorithm will perform pseudo-pooling between individually processed samples. 
                     | This argument has no effect if only 1 sample is provided, and pool does not affect error rates, 
                     | which are always estimated from pooled observations across samples.
``selfConsist``      | if TRUE, the algorithm will alternate between sample inference and error rate estimation until convergence
``qualityType``      | means to attempt to auto-detect the fastq quality encoding. 
                     | This may fail for PacBio files with uniformly high quality scores, in which case use 'FastqQuality'
==================== ============

see :ref:`default settings <dada2_defaults>`

____________________________________________________

.. _dada2_merge_pairs:

MERGE PAIRS [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `mergePairs <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to merge paired-end reads. 
Outputs merged fasta files into ``denoised_assembled.dada2`` directory.

==================== ============
Setting               Tooltip
==================== ============
``minOverlap``       | the minimum length of the overlap required for merging the forward and reverse reads.
``maxMismatch``      | the maximum mismatches allowed in the overlap region
``trimOverhang``     | if TRUE, overhangs in the alignment between the forwards and reverse read are trimmed off. 
                     | Overhangs are when the reverse read extends past the start of the forward read, 
                     | and vice-versa, as can happen when reads are longer than the amplicon and read 
                     | into the other-direction primer region
``justConcatenate``  | if TRUE, the forward and reverse-complemented reverse read are concatenated  
                     | rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them
==================== ============

see :ref:`default settings <dada2_defaults>`

.. _dada2_chimeras:

____________________________________________________

CHIMERA FILTERING [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `removeBimeraDenovo <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to remove chimeras. 
Outputs filtered fasta files into ``chimeraFiltered_out.dada2`` and final ASVs to ``ASVs_out.dada2`` directory.

==================== ============
Setting               Tooltip
==================== ============
``method``           | 'consensus' - the samples are independently checked for chimeras, and a consensus 
                     | decision on each sequence variant is made. 
                     | If 'pooled', the samples are all pooled together for chimera identification. 
                     | If 'per-sample', the samples are independently checked for chimeras
==================== ============

see :ref:`default settings <dada2_defaults>`

____________________________________________________

.. _dada2_taxonomy:

ASSIGN TAXONOMY [ASVs workflow] 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 `assignTaxonomy <https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf>`_ function to classify ASVs. 
Outputs classified fasta files into ``taxonomy_out.dada2`` directory.

==================== ============
Setting               Tooltip
==================== ============
``minBoot``          | the minimum bootstrap confidence for assigning a taxonomic level
``tryRC``            | the reverse-complement of each sequences will be used for classification 
                     | if it is a better match to the reference sequences than the forward sequence
``refFasta``         | select a reference database fasta file for taxonomy annotation
                     | `Download DADA2-formatted reference databases here <https://benjjneb.github.io/dada2/training.html>`_
==================== ============

see :ref:`default settings <dada2_defaults>`

____________________________________________________

.. _otupipe:

OTUs workflow panel
--------------------

.. note::
  This OTU workflow works with paired-end (e.g. Illumina, MGI-Tech) as well as single-end reads (e.g. PacBio, assembled Illumina reads)

|otu_main|

This automated workflow is mostly based on `vsearch <https://github.com/torognes/vsearch>`_ (`Rognes et. al 2016 <https://peerj.com/articles/2584/>`_) [`manual <_static/vsearch_2.18.0_manual.pdf>`_.]
 | Note that ``demultiplexing``, ``reorient`` and ``remove primers`` steps are optional. Nevertheless, it is advisable to :ref:`reorient <reorinet>` your reads (to 5'-3') and :ref:`remove primers <remove_primers>` before proceeding.

 
.. _otupipe_defaults:

**Default options:**

================================================= =========================
Analyses step                                     Default setting
================================================= =========================
:ref:`DEMULTIPLEX <demux>` (optional)              --
:ref:`REORIENT <reorinet>` (optional)              --
:ref:`REMOVE PRIMERS <remove_primers>` (optional)  --
:ref:`MERGE READS <merge_pairs>`                   | ``min overlap`` = 12
                                                   | ``min length`` = 32
                                                   | ``allow merge stagger`` = TRUE 
                                                   | ``include only R1`` = FALSE 
                                                   | ``max diffs`` = 20
                                                   | ``max Ns`` = 0
                                                   | ``max len`` = 600
                                                   | ``keep disjoined`` = FALSE 
                                                   | ``fastq qmax`` = 41
:ref:`QUALITY FILTERING <qual_filt>`               | ``maxEE`` = 1
                                                   | ``maxN`` = 0
                                                   | ``minLen`` = 32
                                                   | ``max length`` = undefined
                                                   | ``qmax`` = 41
                                                   | ``qmin`` = 0
                                                   | ``maxee rate`` = undefined
                                                   | ``minsize`` = 1
:ref:`CHIMERA FILTERING <chimeras>`                | ``pre-cluster`` = 0.98
                                                   | ``min unique size`` = 1
                                                   | ``denovo`` = TRUE 
                                                   | ``reference based`` = undefined
                                                   | ``abundance skew`` = 2
                                                   | ``min h`` = 0.28
:ref:`ITS Extractor <itsx>`                        | ``organisms`` = Fungi 
                                                   | ``regions`` = all
                                                   | ``partial`` = 50
                                                   | ``e-value`` = 1e-5
                                                   | ``scores`` = 0
                                                   | ``domains`` = 2
                                                   | ``complement`` = TRUE 
                                                   | ``only full`` = FALSE
                                                   | ``truncate`` = TRUE 
:ref:`CLUSTERING <clustering>`                     | ``OTU type`` = centroid
                                                   | ``similarity threshold`` = 0.97
                                                   | ``strands`` = both
                                                   | ``min OTU size`` = 2
                                                   | ``similarity type`` = 2
                                                   | ``sequence sorting`` = cluster_size
                                                   | ``centroid type`` = similarity
                                                   | ``max hits`` = 1
                                                   | ``relabel`` = sha1
                                                   | ``mask`` = dust
                                                   | ``dbmask`` = dust
                                                   | ``output UC`` = FALSE
:ref:`ASSGIN TAXONOMY <taxonomy>`                  | ``database file`` = select a database
                                                   | ``task`` = blastn
                                                   | ``strands`` = both
================================================= =========================

____________________________________________________

.. _panels:

ANALYSES PANELS
===============

.. _demux:

DEMULTIPLEX
------------
If data is **multiplexed, the first step would be demultiplexing** (using `cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_ (`Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_)).
This is done based on the user specified :ref:`indexes file <indexes>`, which includes molecular identifier sequences (so called indexes/tags/barcodes) per sample. 
Note that reverse complementary matches will also be searched. 

Supported file format for the input data are **fastq** or **fasta**.
**Outputs** are fastq/fasta files per sample in ``demultiplexed_out`` directory. Indexes are **truncated** from the sequences. 

.. _indexes:

Indexes file example (fasta formatted)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note::
  Only **IUPAC codes** are allowed.

1. **Demultiplexing using single indexes**:

 | >sample1
 | AGCTGCACCTAA
 | >sample2
 | AGCTGTCAAGCT
 | >sample3
 | AGCTTCGACAGT
 | >sample4
 | AGGCTCCATGTA
 | >sample5
 | AGGCTTACGTGT
 | >sample6
 | AGGTACGCAATT

2. **Demultiplexing using dual (paired) indexes:**

.. note::
 **IMPORTANT!** reverse indexes will be automatically oriented to 5'-3' (for the search); so you can simply copy-paste the indexes from your lab protocol.


| >sample1
| AGCTGCACCTAA...AGCTGCACCTAA
| >sample2
| AGCTGTCAAGCT...AGCTGTCAAGCT
| >sample3
| AGCTTCGACAGT...AGCTTCGACAGT
| >sample4
| AGGCTCCATGTA...AGGCTCCATGTA
| >sample5
| AGGCTTACGTGT...AGGCTTACGTGT
| >sample6
| AGGTACGCAATT...AGGTACGCAATT

.. note::
 Anchored indexes (https://cutadapt.readthedocs.io/en/stable/guide.html#anchored-5adapters) with ^ symbol are not supported in PipeCraft demultiplex GUI panel. 

 DO NOT USE, e.g. 

 | >sample1
 | ^AGCTGCACCTAA
 | 
 | >sample1
 | ^AGCTGCACCTAA...AGCTGCACCTAA

|

How to compose indexes.fasta 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In Excel (or any alternative program); 
first column represents sample sames,
second (and third) column represent indexes (or index combinations) per sample:

Exaples::

     sample1	AGCTGCACCTAA
     sample2	AGCTGTCAAGCT
     sample3	AGCTTCGACAGT 
     sample4	AGGCTCCATGTA
     sample5	AGGCTTACGTGT
     sample6	AGGTACGCAATT

or ::

     sample1	AGCTGCACCTAA	AGCTGCACCTAA
     sample2	AGCTGTCAAGCT	AGCTGTCAAGCT
     sample3	AGCTTCGACAGT	AGCTTCGACAGT
     sample4	AGGCTCCATGTA	AGGCTCCATGTA
     sample5	AGGCTTACGTGT	AGGCTTACGTGT
     sample6	AGGTACGCAATT	AGGTACGCAATT

Copy those two (or three) columns to text editor that support regular expressions, such as NotePad++ or Sublime Text.
If using **PAIRED** indexes (three columns), proceed to bullet no. 5

* single-end indexes:

  #. Open 'find & replace'
     Find ^   (which denotes the beginning of each line).
     Replace with >  (and DELETE THE LAST > in the beginning of empty row).

  #. Find \\t   (which denotes tab).
     Replace with \\n   (which denotes the new line).

     **FASTA FORMATTED (single-end indexes) indexes.fasta file is ready; SAVE the file.**


* Only for paired-indexes:

  #. Open 'find & replace':
     Find ^   (denotes the beginning of each line);
     replace with >  (and DELETE THE LAST > in the beginning of empty row).

  #. Find .*\\K\\t (which captures the second tab);
     replace with ... (to mark the linked paired-indexes). 

  #. Find \\t (denotes the tab);
     replace with \\n (denotes the new line).

     **FASTA FORMATTED (paired indexes) indexes.fasta file is ready; SAVE the file.**

____________________________________________________

.. _reorinet:

REORIENT
--------

`fqgrep <https://github.com/indraniel/fqgrep>`_   `Indraniel Das 2011 <https://github.com/indraniel/fqgrep>`_

Sequences are often in both, 5’-3’ and 3’-5’, orientations in the raw sequencing data sets. If the data still contains PCR primers that were used to generate amplicons, then by specifying these PCR primers (up to 13 pairs allowed), this panel will perform sequence reorientation of all sequences to 5’-3’. For reorienting, first the forward primer will be searched (fwd specified in 5’-3’ orientation as for PCR) and if detected then the read is considered as forward complementary (5’-3’). Then the reverse primer (specified in 3’-5’ orientation, as for PCR) will be searched from the same input data and if detected, then the read is considered to be in reverse complementary orientation (3’-5’). Latter reads will be transformed to 5’-3’ orientation and merged with other 5’-3’ reads. Note that for paired-end data, R1 files will be reoriented to 5’-3’ but R2 reads will be reoriented to 3’-5’ in order to merge paired-end reads (see below, Merging paired end sequences).
At least one of the PCR primers must be found in the sequence. For example, read will be recorded if forward primer was found even though reverse primer was not found (and vice versa). Sequence is discarded if none of the PCR primers are found. Sequences that contain multiple forward or reverse primers (multi-primer artefacts) are discarded as it is highly likely that these are chimeric sequences. Reorienting sequences will not remove primer strings from the sequences. Primers may be removed in the “Cut primers” panel (see below). Note that for single-end data, sequences will be reoriented also during the ‘cut primers’ process (see below); therefore this step may be skipped when working with single-end data (such as data from PacBio machines OR already assembled paired-end data).

Reorienting reads may be relevant for generating ASVs with DADA2 (this step is performed by default when selecting DADA2 full workflow to process the data in PipeCraft) as reverse complement sequences will represent separate ASVs. In the clustering step of an OTU pipeline, both strands of the sequences can be compared prior forming OTUs; thus this step may be skipped in the OTU pipeline. 

Supported file formats for paired-end input data are only fastq (extensions must be .fastq or .fq), but also fasta (extensions must be .fasta, .fa, .fas) for single-end data.

|

Specifics of the panel workflow: user has to specify the PCR primers that were used to generate amplicons; IUPAC codes for degenerate bases are allowed (example: CCTCCSCTTANTDAT). ‘Any base’ can be marked with N or I. If the PCR primer strings are not found in the input data, then the warning is displayed and no output for that file is generated. Fastq formatted files are supported for paired-end data (Illumina, MGI-Tech); fasta and fastq formatted files for single-end data (PacBio, Nanopore, Ion Torrent, 454). Paired-end data must contain strings ‘R1’ and ‘R2’ in corresponding file names. Gz or zip compressed files are supported, but decompressed using pigz prior analyses.  Fqgrep is used for searching the primers in the input files; fastx-toolkit is used to reverse complement the reads when needed; seqkit is used to synchronize the paired-end reads. Running the process several times in the same directory will overwrite all outputs.
Summary of the sequence counts can be found in the ‘seq_count_summary.txt’ file.

Supported file format for the input data are **fastq** or **fasta**.
**Outputs** are fastq/fasta files in **reoriented_out** directory. Primers are **not truncated** from the sequences. 

____________________________________________________

.. _remove_primers:

CUT PRIMERS
-----------

If the input data contains PCR primers (or e.g. adapters), these can be removed in the ‘CUT PRIMERS’ panel. `Cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_ (`Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_) is used for that. 

Up to 13 forward and reverse primers may be specified. For generating OTUs or ASVs, it is recommended to truncate the primers from the reads. Sequences where PCR primer strings were not detected are discarded by default (but stored in ‘untrimmed’ directory). 
Reverse complementary search of the primers in the sequences is also performed. Thus, primers are clipped from both 5’-3’ and 3’-5’ oriented reads. However, note that **paired-end reads will not be reoriented** to 5’-3’ during this process, 
but **single-end reads will be reoriented** to 5’-3’ (thus no extra reorient step needed for single-end data).

For paired-end data, the seqs_to_keep option should be left as default (‘keep_all’). This will output sequences where at least one primer has been clipped. 
‘keep_only_linked’ option outputs only sequences where both the forward and reverse primers are found (i.e. 5’-forward…reverse-3’).

Supported file format for the input data are **fastq** or **fasta**.
**Outputs** are fastq/fasta files in **primersCut_out** directory. Primers are **truncated** from the sequences. 

____________________________________________________

.. _qual_filt:

QUALITY FILTERING
------------------

Quality filtering of the **fastq** formatted files can be conducted using `vsearch <https://github.com/torognes/vsearch>`_ (`Rognes et. al 2016 <https://peerj.com/articles/2584/>`_) 
or `trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_ (`Bolger et al. 2014 <https://doi.org/10.1093/bioinformatics/btu170>`_). 

Supported file format for the input data is **fastq**.
**Outputs** are fastq files in **qualFiltered_out** directory.

================================ =========================
**vsearch** setting              Tooltip
================================ =========================
``maxEE``                        | maximum number of expected errors per sequence. Sequences with higher error rates will be discarded
``maxN``                         | discard sequences with more than the specified number of Ns
``minLen``                       | minimum length of the filtered output sequence
``max length``                   | discard sequences with more than the specified number of bases
``qmax``                         | specify the maximum quality score accepted when reading FASTQ files. 
                                 | The default is 41, which is usual for recent Sanger/Illumina 1.8+ files. For PacBio data use 93
``qmin``                         | the minimum quality score accepted for FASTQ files. The default is 0, which is 
                                 | usual for recent Sanger/Illumina 1.8+ files. Older formats may use scores between -5 and 2
``maxee rate``                   | discard sequences with more than the specified number of expected errors per base
``minsize``                      | discard sequences with an abundance lower than the specified value
================================ =========================

================================ =========================
**trimmomatic** setting          Tooltip
================================ =========================
``window size``                  | the number of bases to average base qualities. 
                                 | Starts scanning at the 5'-end of a sequence and trimms the read once the 
                                 | average required quality (required_qual) within the window size falls below the threshold
``required quality``             | the average quality required for selected window size
``min length``                   | minimum length of the filtered output sequence
``leading_qual_threshold``       | quality score threshold to remove low quality bases from the beginning of the read. 
                                 | As long as a base has a value below this threshold the base is removed and the next base will be investigated
``trailing_qual_threshold``      | quality score threshold to remove low quality bases from the end of the read. 
                                 | As long as a base has a value below this threshold the base is removed and the next base will be investigated
``phred``                        | phred quality scored encoding. 
                                 | Use phred64 if working with data from older Illumina (Solexa) machines
================================ =========================

____________________________________________________

.. _merge_pairs:

ASSEMBLE PAIRED-END reads 
--------------------------

Paired-end sequences (such as those from Illumina or MGI-Tech platforms) may be merged using `vsearch <https://github.com/torognes/vsearch>`_ (`Rognes et. al 2016 <https://peerj.com/articles/2584/>`_)

As an additional in-built module of PipeCraft, users may include unassembled R1 reads to the set of assembled reads per sample. 
This may be relevant when working with e.g. ITS2 sequences, because the ITS2 region in some taxa is too long for paired-end assembly, 
therefore discarded completely after the assembly process. Thus, including also unassembled R1 reads, partial ITS2 sequences for 
these taxa will be represented in the final output (in that case, when applying ITSx module, include ‘partial’ option).

Supported file format for the input data is **fastq**.
**Outputs** are fastq files in **assembled_out** directory.

================================ =========================
Setting                          Tooltip
================================ =========================
``min overlap``
``min length``
``allow merge stagger``
``include only R1``
``max diffs``
``max Ns``
``max len``
``keep disjoined``
``fastq qmax``
================================ =========================

____________________________________________________

.. _chimeras:

CHIMERA FILTERING
-----------------

De-novo as well as reference database based chimera filtering is supported through `vsearch <https://github.com/torognes/vsearch>`_ (`Rognes et. al 2016 <https://peerj.com/articles/2584/>`_). 
Users may include any custom fasta formatted database for reference based chimera filtering.

Chimera filtering is performed by **sample-wise approach** (each sample is treated separately). 

Supported input file format:
    Supported file format for the input data are **fastq** or **fasta**. [Fastq inputs will be converted to fasta]
    **Outputs** are fasta files in **chimera_Filtered_out** directory.

================================ =========================
Setting                          Tooltip
================================ =========================
``pre-cluster``
``min unique size``
``denovo``
``reference based``
``abundance skew``
``min h``
================================ =========================

____________________________________________________

.. _itsx:

ITS Extractor
-------------

`ITS Extractor <https://microbiology.se/software/itsx/>`_ (`Bengtsson-Palme et al. 2013 <https://doi.org/10.1111/2041-210X.12073>`_)

Supported file format for the input data are **fastq** or **fasta**. [Fastq inputs will be converted to fasta]
**Outputs** are fasta files in **ITSx_out** directory.

================================ =========================
Setting                          Tooltip
================================ =========================
``organisms`` 
``regions``
``partial``
``e-value``
``scores``
``domains``
``complement`` 
``only full``
``truncate`` 
================================ =========================

.. _clustering:

CLUSTERING
----------

`vsearch <https://github.com/torognes/vsearch>`_ (`Rognes et. al 2016 <https://peerj.com/articles/2584/>`_)

Supported file format for the input data is **fasta**.
**Outputs** are OTUs.fasta and OTU_table.txt files in **clustering_out** directory.

================================ =========================
Setting                          Tooltip
================================ =========================
``OTU type``
``similarity threshold``
``strands``
``min OTU size``
``similarity type``
``sequence sorting``
``centroid type``
``max hits``
``relabel``
``mask``
``dbmask``
``output UC``
================================ =========================

____________________________________________________

.. _taxonomy:

ASSIGN TAXONOMY
---------------

Implemented tools for taxonomy annotation:

`BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_ (`Camacho et al. 2009 <https://doi.org/10.1186/1471-2105-10-421>`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BLAST search query sequences againt selected database. 

================================ =========================
Setting                          Tooltip
================================ =========================
 ``database file``               | select a database file in fasta format.
                                 | Fasta format will be automatically converted to BLAST database
``task``                         | BLAST search settings according to blastn or megablast
``strands``                      | query strand to search against database. Both = search also reverse complement
``e-value``                      | a parameter that describes the number of hits one can expect to see 
                                 | by chance when searching a database of a particular size. 
                                 | The lower the e-value the more 'significant' the match is
``word size``                    | the size of the initial word that must be matched between the database and the query sequence
``reward``                       | reward for a match
``penalty``                      | penalty for a mismatch
``gap open``                     | cost to open a gap
``gap extend``                   | cost to extend a gap
================================ =========================

____________________________________________________

.. _expert_mode:

Expert-mode (PipeCraft console)
===============================

**UNDER CONSTRUCTION**

|console|





