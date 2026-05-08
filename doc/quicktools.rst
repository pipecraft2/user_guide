.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 50
  :target: https://github.com/pipecraft2/user_guide 


.. |cut_primers_expand_example| image:: _static/cut_primers_expand_example.png
  :width: 600

.. meta::
    :description lang=en:
        PipeCraft manual. PipeCraft in a Graphical User Interface software for metabarcoding data analyses

.. _quicktools:

================================================
Individual steps (Quick Tools) |PipeCraft2_logo|
================================================

QuickTools provide a list of processes that can be used to perform individual steps of the analysis (i.e., perform a custom pipeline).
They are accessed by pressing the ``Quick Tools`` button on the right-ribbon interface.

.. |quicktools_button| image:: _static/quicktools_button.png
  :width: 600

|quicktools_button|

.. _demux:

__________________________________________________

DEMULTIPLEXING
==============

| `Download example set here for trying demultiplexing <https://zenodo.org/records/18770850/files/demux_example.zip?download=1>`_ and unzip it. 

If the data is **multiplexed, the first step would be demultiplexing** (using `cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_ (`Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_)).
This is done based on the user specified :ref:`indexes file <indexes>`, which includes molecular identifier sequences 
(so-called indexes/tags/barcodes) per sample. 
Note that reverse complementary matches will also be searched. 

| **Fastq/fasta** formatted paired-end and single-end data are supported.
| **Outputs** are fastq/fasta files per sample in ``demultiplexed_out`` directory. Indexes are **truncated** from the sequences. 
| Paired-end samples get ``.R1`` and ``.R2`` read identifiers.
| **unknown.fastq** file contain sequences where specified index combinations were not found. 

.. note:: 

  **When using paired indexes**, then sequences with all possible index combinations will be outputted to 'unnamed_index_combinations' dir.
  That means, if, for example, your sample_1 is indexed with *indexFwd_1-indexRev_1* and 
  sample_2 with *indexFwd_2-indexRev_2*, then files with *indexFwd_1-indexRev_2* and *indexFwd_2-indexRev_1*
  are also written (although latter index combinations were not used in the lab to index any sample [i.e. represent tag-switches]). 
  Simply remove those files if not needed or use to estimate tag-switching error if relevant. 


.. _demux_settings:

+--------------------+----------------------------------------------------------------------+
| Setting            | Tooltip                                                              |
+====================+======================================================================+
|| ``index file``    || select your fasta formatted indexes file for demultiplexing         |
||                   || (:ref:`see guide here <indexes>`), where fasta headers are sample   |
||                   || names, and sequences are sample specific index or index combination |
+--------------------+----------------------------------------------------------------------+
| ``index mismatch`` | allowed mismatches during the index search                           |
+--------------------+----------------------------------------------------------------------+
|| ``overlap``       || number of overlap bases with the index. Recommended overlap is the  |
||                   || maximum length of the index for confident sequence assignments to   |
||                   || samples                                                             |
+--------------------+----------------------------------------------------------------------+
|| ``search window`` || the index search window size. The default 35 means that the forward |
||                   || index is searched among the first 35 bp and the reverse index among |
||                   || the last 35 bp. This search restriction prevents random index       |
||                   || matches in the middle of the sequence                               |
+--------------------+----------------------------------------------------------------------+
|| ``no indels``     || do not allow insertions or deletions is primer search. Mismatches   |
||                   || are the only type of errors accounted in the error rate parameter   |
+--------------------+----------------------------------------------------------------------+
| ``min length``     | minimum length of the output sequence                                |
+--------------------+----------------------------------------------------------------------+


.. note::

 Heterogenity spacers or any redundant base pairs attached to index sequences do not affect demultiplexing. 
 Indexes are trimmed from the best matching position.

.. _indexes:

Indexes file example (fasta formatted)
--------------------------------------

.. note::

  Only **IUPAC codes** are allowed in the sequences. Avoid using '.' in the sample names (e.g. instead of sample.1, use sample_1)

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

2. **Demultiplexing using paired (dual) indexes:**

.. important::

 **IMPORTANT!** The reverse indexes **must be in the 3'-5' orientation** in the indexes file when doing demultiplexing in PipeCraft, because
 reverse indexes are automatically oriented to 5'-3' under the hood. This facilitates the simple copy-paste of the indexes from the lab protocol.
 Therefore, **if you have pre-compliled indexes file**, so, 
 that you have reverse indexes already reverse-comlemented (5'-3' orientation), then the demultiplexing will fail (all will be unknown.fastq).


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

 Anchored indexes (https://cutadapt.readthedocs.io/en/stable/guide.html#anchored-5adapters) with ^ symbol are **not supported** 
 in PipeCraft demultiplex GUI panel. 

 DO NOT USE, e.g. 

 | >sample1
 | ^AGCTGCACCTAA
 | 
 | >sample1
 | ^AGCTGCACCTAA...AGCTGCACCTAA

 Instead, specify ``search window`` = 0.

|

How to compose indexes.fasta 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In Excel (or any alternative program); 
first column represents sample names,
second (and third) column represent indexes (or index combinations) per sample:

Example of **single-end indexes** ::

     sample1	AGCTGCACCTAA
     sample2	AGCTGTCAAGCT
     sample3	AGCTTCGACAGT 
     sample4	AGGCTCCATGTA
     sample5	AGGCTTACGTGT
     sample6	AGGTACGCAATT

Example of **paired indexes** ::

     sample1	AGCTGCACCTAA	AGCTGCACCTAA
     sample2	AGCTGTCAAGCT	AGCTGTCAAGCT
     sample3	AGCTTCGACAGT	AGCTTCGACAGT
     sample4	AGGCTCCATGTA	AGGCTCCATGTA
     sample5	AGGCTTACGTGT	AGGCTTACGTGT
     sample6	AGGTACGCAATT	AGGTACGCAATT

Copy those two (or three) columns to text editor that support regular expressions, such as **NotePad++ or Sublime Text**.

* single-end indexes:

  #. Open 'find & replace'
     Find ^   (which denotes the beginning of each line).
     Replace with >  (and DELETE THE LAST > in the beginning of empty row).

  #. Find \\t   (which denotes tab).
     Replace with \\n   (which denotes the new line).

     **FASTA FORMATTED (single-end indexes) indexes.fasta file is ready; SAVE the file.**

* Paired indexes:

  #. Open 'find & replace':
     Find ^   (denotes the beginning of each line);
     replace with >  (and DELETE THE LAST > in the beginning of empty row).

  #. Find .*\\K\\t (which captures the second tab);
     replace with ... (to mark the linked paired-indexes). 

  #. Find \\t (denotes the tab);
     replace with \\n (denotes the new line).

     **FASTA FORMATTED (paired indexes) indexes.fasta file is ready; SAVE the file.**

____________________________________________________

.. _remove_primers:

CUT PRIMERS
===========

If the input data contains PCR primers (or e.g. adapters), these can be removed in the ``CUT PRIMERS`` panel.
CUT PRIMERS processes relies on `cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_ (`Martin 2011 <https://doi.org/10.14806/ej.17.1.200>`_). 

For generating OTUs or ASVs, it is recommended to truncate the primers from the reads 
(**unless ITS Extractor is used** to remove flanking primer binding regions from ITS1/ITS2/full ITS; 
in that case keep the primers for better detection of the 18S, 5.8S and/or 28S regions). 
Sequences where PCR primer strings were not detected are discarded by default (but stored in 'untrimmed' directory). 
Reverse complementary search of the primers in the sequences is also performed. 
Thus, primers are clipped from both 5'-3' and 3'-5' oriented reads. However, note that 
**paired-end reads will not be reoriented** to 5'-3' during this process, 
but **single-end reads will be reoriented** to 5'-3'.

.. note::

 For paired-end data, the **seqs_to_keep option should be left as default ('keep_all')**. 
 This will output sequences where at least one primer has been clipped. 
 'keep_only_linked' option outputs only sequences where both the forward and reverse primers are found (i.e. 5'-forward…reverse-3'). 
 'keep_only_linked' may be used for single-end data to keep only **full-length amplicons**.


|cut_primers_expand_example|

Example above: Forward primer has 19 bp and reverse 20 bp - to keep a bit of flexibility in the primer search, we are requesting the ``min overlap`` of **18 bp** and are allowing maximum of 2 ``mismatches`` . 
Note that too low ``min overlap`` may lead to random matches.

| **Fastq**/**fasta** formatted paired-end and single-end data are supported.
| **Outputs** are fastq/fasta files in ``primersCut_out`` directory. Primers are **truncated** from the sequences. 


.. admonition:: when working with your own ITS data ... 

  ... and applying the **ITSx** step, then note that cutting primers process may be skipped, since those regions are removed in the ITS subregion extraction process. 
  

+----------------------+-----------------------------------------------------------------------+
| Setting              | Tooltip                                                               |
+======================+=======================================================================+
|| ``forward primers`` || specify forward primer **(5'-3')**; IUPAC codes allowed; add up to   |
||                     || 13 primers                                                           |
+----------------------+-----------------------------------------------------------------------+
|| ``reverse primers`` || specify reverse primer **(3'-5')**; IUPAC codes allowed; add up to   |
||                     || 13 primers                                                           |
+----------------------+-----------------------------------------------------------------------+
| ``mismatches``       | allowed mismatches in the primer search                               |
+----------------------+-----------------------------------------------------------------------+
|| ``min overlap``     || number of overlap bases with the primer sequence. Partial matches    |
||                     || are allowed, but short matches may occur by chance, leading to       |
||                     || erroneously clipped bases. Specifying higher overlap than the length |
||                     || of primer sequnce will still clip the primer (e.g. primer length is  |
||                     || 22 bp, but overlap is specified as 25 - this does not affect the     |
||                     || identification and clipping of the primer as long as the match is    |
||                     || in the specified mismatch error range)                               |
+----------------------+-----------------------------------------------------------------------+
|| ``seqs to keep``    || keep sequences where at least one primer was found (fwd or rev);     |
||                     || recommended when cutting primers from paired-end data (unassembled), |
||                     || when individual R1 or R2 read lengths are shorther than the expected |
||                     || amplicon length. 'keep_only_linked' = keep sequences if primers are  |
||                     || found in both ends (fwd…rev); discards the read if both primers were |
||                     || not found in this read                                               |
+----------------------+-----------------------------------------------------------------------+
|| ``pair filter``     || **applies only for paired-end data.** 'both', means that a read is   |
||                     || discarded only if both, corresponding R1 and R2, reads do not        |
||                     || contain primer strings (i.e. a read is kept if R1 contains primer    |
||                     || string, but no primer string found in R2 read). Option 'any'         |
||                     || discards the read if primers are not found in both, R1 and R2 reads  |
+----------------------+-----------------------------------------------------------------------+
|| ``no indels``       || do not allow insertions or deletions is primer search. Mismatches    |
||                     || are the only type of errors accounted in the error rate parameter    |
+----------------------+-----------------------------------------------------------------------+

____________________________________________________

|

.. _qual_filt:

QUALITY FILTERING
=================

Quality filtering removes low-quality sequences before downstream analysis.
Keeping only high-quality sequences prevents noisy data from
creating erroneous OTUs/ASVs. Different tools
implement quality filtering in slightly different ways, but the goal is the same: retain
sequences that meet the specified threshold(s). 

Before running quality filtering, it is best to inspect the sequence quality profiles
to see where quality starts to decline and which trimming settings could be appropriate.
Checkout the :ref:`Inspect quality profiles <qualitycheck>` section for a walkthrough
of this step.

| **Fastq** formatted paired-end and single-end data are supported.
| **Outputs** are fastq files in ``qualFiltered_out`` directory.

Below lists the different quality filtering tools implemented in PipeCraft.

.. _qfilt_vsearch:

`vsearch <https://github.com/torognes/vsearch>`_
--------------------------------------------------

Vsearch (*fastq_filter* function) filters reads by calculating the **expected errors** per read 
(``maxee``; sum of per-base error probabilities derived from Phred scores) 
and discards reads exceeding the threshold. 
In addition, reads can be removed based on **ambiguous bases** (``maxNs``) 
and **length constraints** (``min length`` / ``max length``), 
and optionally **truncated** to a fixed length (``trunc length``). 
Applying ``trunc length`` may be helpful to remove low-quality ends of reads before filtering 
(``maxee`` filtering is applied to the truncated reads).

+---------------------+-------------------------------------------------------------------------+
| **vsearch** setting | Tooltip                                                                 |
+=====================+=========================================================================+
|| ``maxee``          || maximum number of expected errors per sequence                         |
||                    || (`see here <https://drive5.com/usearch/manual/exp_errs.html>`_).       |
||                    || Sequences with higher error rates will be discarded                    |
+---------------------+-------------------------------------------------------------------------+
| ``maxNs``           | discard sequences with more than the specified number of Ns             |
+---------------------+-------------------------------------------------------------------------+
| ``min length``      | minimum length of the filtered output sequence                          |
+---------------------+-------------------------------------------------------------------------+
|| ``trunc length``   || truncate sequences to the specified length. Shorter sequences are      |
||                    || discarded; thus if specified, check that 'min length' setting is       |
||                    || lower than 'trunc length' ('min length' therefore has basically no     |
||                    || effect) [empty field = no action taken]                                |
+---------------------+-------------------------------------------------------------------------+
|| ``qmax``           || specify the maximum quality score accepted when reading FASTQ files.   |
||                    || The default is 41, which is usual for recent Sanger/Illumina 1.8+      |
||                    || files. **For PacBio data use 93**                                      |
+---------------------+-------------------------------------------------------------------------+
|| ``max length``     || discard sequences with more than the specified number of bases. Note   |
||                    || NOT be lower than 'trunc length' (otherwise all reads are discared)    |
||                    || [empty field = no action taken] Note that if 'trunc length' setting    |
||                    || is specified, then 'min length' SHOULD BE lower than 'trunc length'    |
||                    || (otherwise all reads are discared)                                     |
||                    ||                                                                        |
+---------------------+-------------------------------------------------------------------------+
|| ``qmin``           || the minimum quality score accepted for FASTQ files. The default is 0,  |
||                    || which is usual for recent Sanger/Illumina 1.8+ files. Older formats    |
||                    || may use scores between -5 and 2                                        |
+---------------------+-------------------------------------------------------------------------+
|| ``maxee rate``     || discard sequences with more than the specified number of expected      |
||                    || errors per base                                                        |
+---------------------+-------------------------------------------------------------------------+
|| ``truncqual``      || tuncate sequences starting from the first base with the specified      |
||                    || base quality score value or lower (0 or empty field = no action taken) |
+---------------------+-------------------------------------------------------------------------+
|| ``truncee``        || truncate sequences so that their total expected error is not higher    |
||                    || than the specified value (0 or empty field = no action taken)          |
||                    ||                                                                        |
+---------------------+-------------------------------------------------------------------------+

| 

.. _qfilt_trimmomatic:

`trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
---------------------------------------------------------------

Trimmomatic 

Trimmomatic trims and filters reads based on **base-quality scores**. 
The main trimming is performed with a **sliding-window** approach: 
Trimmomatic scans from the 5' end and cuts the read once the mean quality in a 
window (``window_size``) drops below the threshold (``required_quality``). 
Optional additional steps remove low-quality bases at the **start** (``leading_qual_threshold``) 
and **end** (``trailing_qual_threshold``) of reads. 
Reads shorter than ``min_length`` after trimming are discarded. 

+------------------------------+-----------------------------------------------------------------------+
| **trimmomatic** setting      | Tooltip                                                               |
+==============================+=======================================================================+
|| ``window_size``             || the number of bases to average base qualities. Starts scanning at    |
||                             || the 5'-end of a sequence and trimms the read once the average        |
||                             || required quality (required_qual) within the window size falls below  |
||                             || the threshold                                                        |
+------------------------------+-----------------------------------------------------------------------+
| ``required_quality``         | the average quality required for selected window size                 |
+------------------------------+-----------------------------------------------------------------------+
| ``min_length``               | minimum length of the filtered output sequence                        |
+------------------------------+-----------------------------------------------------------------------+
|| ``leading_qual_threshold``  || quality score threshold to remove low quality bases from the         |
||                             || beginning of the read. As long as a base has a value below this      |
||                             || threshold the base is removed and the next base will be investigated |
+------------------------------+-----------------------------------------------------------------------+
|| ``trailing_qual_threshold`` || quality score threshold to remove low quality bases from the end of  |
||                             || the read. As long as a base has a value below this threshold the     |
||                             || base is removed and the next base will be investigated               |
+------------------------------+-----------------------------------------------------------------------+
|| ``phred``                   || phred quality scored encoding. Default = 33. Use 64 if working       |
||                             || with data from older Illumina (Solexa) machines                      |
+------------------------------+-----------------------------------------------------------------------+


| 

.. _qfilt_fastp:

`fastp <https://github.com/OpenGene/fastp>`_
--------------------------------------------

fastp uses a **sliding-window** trimming approach, similar to Trimmomatic, (``window_size`` + ``required_qual``) 
to trim reads when local mean quality drops below the threshold. 
It scans reads from the 5' end toward the 3' end and trims the read from the first low-quality window onward (i.e. removes the low-quality 3' tail).
It additionally filters reads based on the fraction of low-quality bases 
(``min_qual`` + ``min_qual_thresh``), **ambiguous bases** (``maxNs``), **minimum/maximum length** 
(``min_length`` / ``max_length``).

fastp can also trim/remove reads affected by two common artifacts:
**polyG trimming** removes artificial poly-G tails that can appear in some 
Illumina NextSeq/NovaSeq two-colour chemistry runs when signal drops and bases are miscalled as long runs of ``G``.
However, note that when clipping primers, those poly-G artifacts are removed, as the primer sequences are recorded before the signal drops.
**low-complexity filter** removes reads dominated by repetitive/low-information sequence 
(e.g. homopolymers like ``AAAAAA`` or simple repeats). 

+----------------------------+-----------------------------------------------------------------------+
| **fastp** setting          | Tooltip                                                               |
+============================+=======================================================================+
| ``window_size``            | the window size for calculating mean quality                          |
+----------------------------+-----------------------------------------------------------------------+
| ``required_qual``          | the mean quality requirement per sliding window (window_size)         |
+----------------------------+-----------------------------------------------------------------------+
|| ``min_qual``              || the quality value that a base is qualified. Default 15 means phred   |
||                           || quality >=Q15 is qualified                                           |
+----------------------------+-----------------------------------------------------------------------+
| ``min_qual_thresh``        | how many percents of bases are allowed to be unqualified (0-100)      |
+----------------------------+-----------------------------------------------------------------------+
| ``maxNs``                  | discard sequences with more than the specified number of Ns           |
+----------------------------+-----------------------------------------------------------------------+
|| ``min_length``            || minimum length of the filtered output sequence. Shorter sequences    |
||                           || are discarded                                                        |
+----------------------------+-----------------------------------------------------------------------+
|| ``max_length``            || reads longer than 'max length' will be discarded, default 0 means no |
||                           || limitation                                                           |
+----------------------------+-----------------------------------------------------------------------+
|| ``trunc_length``          || truncate sequences to specified length. Shorter sequences are        |
||                           || discarded; thus check that 'min length' setting is lower than 'trunc |
||                           || length'                                                              |
+----------------------------+-----------------------------------------------------------------------+
|| ``aver_qual``             || if one read's average quality score <'aver_qual', then this          |
||                           || read/pair is discarded. Default 0 means no requirement               |
+----------------------------+-----------------------------------------------------------------------+
|| ``low_complexity_filter`` || enables low complexity filter and specify the threshold for low      |
||                           || complexity filter. The complexity is defined as the percentage of    |
||                           || base that is different from its next base (base[i] != base[i+1]).    |
||                           || E.g. vaule 30 means then 30% complexity is required. Not specified = |
||                           || filter not applied                                                   |
+----------------------------+-----------------------------------------------------------------------+


| 

.. _qfilt_dada2:

`DADA2 <https://github.com/benjjneb/dada2>`_ ('filterAndTrim' function)
-----------------------------------------------------------------------

DADA2 (*filterAndTrim* function) filters reads based on `expected errors <https://drive5.com/usearch/manual/exp_errs.html>`_ 
(``maxEE``; same as vsearch) and **ambiguous bases** (``maxN``).
It also allows for **truncation** of reads at the first instance of a quality 
score less than or equal to ``truncQ`` (applied to both R1 and R2 reads).
Reads shorter than ``minLen`` after truncation are discarded.
``truncLen`` / ``truncLen_R2`` options truncate reads to a **fixed number of bases** before ``maxEE`` filtering. 
This may be helpful to remove low-quality ends of reads before filtering.

+-------------------+-------------------------------------------------------------------------+
| **DADA2** setting | Tooltip                                                                 |
+===================+=========================================================================+
|| ``maxEE``        || discard sequences with more than the specified number of expected      |
||                  || errors                                                                 |
+-------------------+-------------------------------------------------------------------------+
|| ``maxN``         || discard sequences with more than the specified number of N's           |
||                  || (ambiguous bases)                                                      |
+-------------------+-------------------------------------------------------------------------+
|| ``minLen``       || remove reads with length less than minLen. minLen is enforced after    |
||                  || all other trimming and truncation                                      |
+-------------------+-------------------------------------------------------------------------+
|| ``truncQ``       || truncate reads at the first instance of a quality score less than or   |
||                  || equal to truncQ                                                        |
+-------------------+-------------------------------------------------------------------------+
|| ``truncLen``     || truncate reads after truncLen bases (applies to **R1 reads** when      |
||                  || working with **paired-end** data). Reads shorter than this are         |
||                  || discarded. Explore quality profiles (with QualityCheck module) and     |
||                  || see whether poor quality ends needs to be truncated                    |
+-------------------+-------------------------------------------------------------------------+
|| ``truncLen_R2``  || applies only for **paired-end** data. Truncate **R2 reads** after      |
||                  || truncLen bases. Reads shorter than this are discarded. Explore         |
||                  || quality profiles (with QualityCheck module) and see whether poor       |
||                  || quality ends needs to truncated                                        |
+-------------------+-------------------------------------------------------------------------+
|| ``maxLen``       || remove reads with length greater than maxLen. maxLen is enforced on    |
||                  || the raw reads. In dada2, the default = Inf, but here set as 9999       |
+-------------------+-------------------------------------------------------------------------+
|| ``minQ``         || after truncation, reads contain a quality score below minQ will be     |
||                  || discarded                                                              |
+-------------------+-------------------------------------------------------------------------+
|| ``matchIDs``     || applies only for **paired-end** data. If TRUE, then double-checking    |
||                  || (with seqkit pair) that only paired reads that share ids are outputted |
+-------------------+-------------------------------------------------------------------------+

____________________________________________________

| 

.. _merge_pairs:

ASSEMBLE PAIRED-END reads 
=========================

Assemble/merge paired-end sequences (such as those from Illumina or MGI-Tech platforms). 

**Fastq** formatted paired-end data is supported.
**Outputs** are fastq files in ``assembled_out`` directory.


.. _merge_vsearch:

`vsearch <https://github.com/torognes/vsearch>`_
--------------------------------------------------

vsearch (*--fastq_mergepairs* function) merges paired-end reads based on **overlap** between the reads.
The best-supported overlap is accepted only if it meets the spacified  
constraints (e.g., ``min_overlap`` and ``max_diffs``). 
In the overlap region, base conflicts are resolved using 
the read quality scores (higher-quality base is preferred). 

``include_only_R1`` represents additional in-built option in PipeCraft. If TRUE, 
unassembled R1 reads will be included to the set of assembled reads per sample. 
This may be relevant when working with e.g. ITS2 sequences, because the ITS2 region in some 
taxa is too long for paired-end assembly using short-read sequencing technology. 
Therefore longer ITS2 sequences are discarded completely after the assembly process. 
Thus, including also unassembled R1 reads (``include_only_R1`` = TRUE), partial ITS2 sequences for 
these taxa will be represented in the final output. But when using :ref:`ITSx <itsextractor>`, 
keep ``only_full`` = FALSE and include ``partial`` = 50.

+--------------------------+-------------------------------------------------------------------------+
| Setting                  | Tooltip                                                                 |
+==========================+=========================================================================+
| ``min overlap``          | minimum overlap between the merged reads                                |
+--------------------------+-------------------------------------------------------------------------+
| ``min length``           | minimum length of the merged sequence                                   |
+--------------------------+-------------------------------------------------------------------------+
|| ``allow merge stagger`` || when TRUE, vsearch will also attempt to merge **staggered** read pairs |
||                         || (pairs with an overhang rather than a clean overlap). This can occur   |
||                         || when the insert/fragment is very short relative to the sequencing run. |
||                         || In that situation,                                                     |
||                         || after reverse-complementing R2 the reads can “pass” each other, so one |
||                         || read extends beyond the start of the other (an overhang). Enabling     |
||                         || this option allows such short-insert pairs to be merged. Short inserts |
||                         || often come with adapter read-through; adapter/primer trimming helps,   |
||                         || but very short inserts can still occur, so this option can be useful.  |
+--------------------------+-------------------------------------------------------------------------+
|| ``include only R1``     || Include unassembled R1 reads to the set of assembled reads per sample. |
||                         || This may be relevant when working with e.g. ITS2 sequences,            |
||                         || because the ITS2 region in some taxa is too long for assembly,         |
||                         || therefore discarded completely after assembly process. Thus, including |
||                         || also unassembled R1 reads, partial ITS2 sequences for these            |
||                         || taxa will be represented in the final output                           |
+--------------------------+-------------------------------------------------------------------------+
|| ``max diffs``           || the maximum number of non-matching nucleotides allowed in the overlap  |
||                         || region                                                                 |
+--------------------------+-------------------------------------------------------------------------+
| ``max Ns``               | discard sequences with more than the specified number of Ns             |
+--------------------------+-------------------------------------------------------------------------+
| ``max length``           | maximum length of the merged sequence                                   |
+--------------------------+-------------------------------------------------------------------------+
| ``keep disjoined``       | output reads that were not merged into separate FASTQ files             |
+--------------------------+-------------------------------------------------------------------------+
|| ``fastq qmax``          || maximum quality score accepted when reading FASTQ files. The default   |
||                         || is 41, which is usual for recent Sanger/Illumina 1.8+ files            |
+--------------------------+-------------------------------------------------------------------------+

|
     

.. _merge_dada2:

`DADA2 <https://github.com/benjjneb/dada2>`_
--------------------------------------------

DADA2 (*mergePairs* function) merges paired-end reads based on **overlap** between the reads.
It allows for **trimming** of overhangs in the alignment between the forwards and reverse read, 
and **concatenation** of the forward and reverse-complemented reverse read with a spacer of 10 Ns.

.. important::

  Here, dada2 will perform also denoising (function 'dada') before assembling paired-end data. 
  Because of that, input sequences (in **fastq** format) must consist of 
  only A/T/C/Gs (**no ambiguous bases (Ns)**); theerefore apply DADA2 merge on quality-filtered reads. 

+----------------------+-----------------------------------------------------------------------+
| Setting              | Tooltip                                                               |
+======================+=======================================================================+
|| ``minOverlap``      || the minimum length of the overlap required for merging the forward   |
||                     || and reverse reads                                                    |
+----------------------+-----------------------------------------------------------------------+
| ``maxMismatch``      | the maximum mismatches allowed in the overlap region                  |
+----------------------+-----------------------------------------------------------------------+
|| ``trimOverhang``    || if TRUE, overhangs in the alignment between the forwards and reverse |
||                     || read are trimmed off. Overhangs are when the reverse read extends    |
||                     || past the start of the forward read, and vice-versa, as can happen    |
||                     || when reads are longer than the amplicon and read into the            |
||                     || other-direction primer region                                        |
+----------------------+-----------------------------------------------------------------------+
|| ``justConcatenate`` || if TRUE, the forward and reverse-complemented reverse read are       |
||                     || concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer    |
||                     || inserted between them                                                |
+----------------------+-----------------------------------------------------------------------+
|| ``pool``            || denoising setting. If TRUE, the algorithm will pool together all     |
||                     || samples prior to sample inference. Pooling improves the detection of |
||                     || rare variants, but is computationally more expensive. If pool =      |
||                     || 'pseudo', the algorithm will perform pseudo-pooling between          |
||                     || individually processed samples.                                      |
+----------------------+-----------------------------------------------------------------------+
|| ``selfConsist``     || denoising setting. If TRUE, the algorithm will alternate between     |
||                     || sample inference and error rate estimation until convergence         |
+----------------------+-----------------------------------------------------------------------+
|| ``qualityType``     || 'Auto' means to attempt to auto-detect the fastq quality encoding.   |
||                     || This may fail for PacBio files with uniformly high quality scores,   |
||                     || in which case use 'FastqQuality'                                     |
+----------------------+-----------------------------------------------------------------------+


.. _chimFilt:

____________________________________________________

|

CHIMERA FILTERING
=================

Chimeras are PCR artifacts that are a combination of two (or more) biological sequences. 
In PipeCraft2 (via vsearch UCHIME/UCHIME3), sequences are first **dereplicated** (identical sequences collapsed),
optionally **pre-clustered** (``pre_cluster``) so that very similar reads are merged and their **size annotations**
reflect the combined abundance (helping to account for residual sequencing errors and providing more robust abundance
information for chimera detection), and can be filtered by a minimum abundance (``min_unique_size``). 

Chimera filtering is performed by **sample-wise approach** (i.e. each sample (input file) is treated separately). 

For **de-novo** detection (``uchime_denovo``), 
candidate sequences are evaluated against more abundant sequences in the same sample; a sequence 
is flagged as chimeric if it can be explained as a mosaic of two "parent" sequences. 

For **reference-based** detection (``uchime_ref``), sequences are compared against a reference database (user-provided). 
Sequences are flagged as chimeras when they are better explained as a mosaic of two reference sequences than by any single reference match.

| **Fastq/fasta** formatted single-end data is supported [fastq inputs will be converted to fasta].
| **Outputs** are fasta files in ``chimera_Filtered_out`` directory.

.. _chimFilt_vsearch:

uchime_denovo
-------------

| Perform chimera filtering with **uchime_denovo** and **uchime_ref** algorithms in `vsearch <https://github.com/torognes/vsearch>`_ 

+----------------------+----------------------------------------------------------------------+
| Setting              | Tooltip                                                              |
+======================+======================================================================+
|| ``pre_cluster``     || identity percentage when performing 'pre-clustering' with           |
||                     || --cluster_size for denovo chimera filtering with --uchime_denovo    |
+----------------------+----------------------------------------------------------------------+
|| ``min_unique_size`` || minimum amount of a unique sequences in a fasta file. If value = 1, |
||                     || then no sequences are discarded after dereplication; if value = 2,  |
||                     || then sequences, which are represented only once in a given file are |
||                     || discarded; and so on                                                |
+----------------------+----------------------------------------------------------------------+
| ``denovo``           | if TRUE, then perform denovo chimera filtering with --uchime_denovo  |
+----------------------+----------------------------------------------------------------------+
|| ``reference_based`` || perform reference database based chimera filtering with             |
||                     || --uchime_ref. Select fasta formatted reference database (e.g. UNITE |
||                     || for ITS reads).                                                     |
||                     || If denovo = TRUE, then reference based chimera filtering will       |
||                     || be performed after denovo.                                          |
||                     ||                                                                     |
+----------------------+----------------------------------------------------------------------+
|| ``abundance_skew``  || the abundance skew is used to distinguish in a threeway alignment   |
||                     || which sequence is the chimera and which are the parents. The        |
||                     || assumption is that chimeras appear later in the PCR amplification   |
||                     || process and are therefore less abundant than their parents. The     |
||                     || default value is 2.0, which means that the parents should be at     |
||                     || least 2 times more abundant than their chimera. Any positive value  |
||                     || equal or greater than 1.0 can be used                               |
+----------------------+----------------------------------------------------------------------+
|| ``min_h``           || minimum score (h). Increasing this value tends to reduce the number |
||                     || of false positives and to decrease sensitivity. Values ranging from |
||                     || 0.0 to 1.0 included are accepted                                    |
+----------------------+----------------------------------------------------------------------+


.. _chimFilt_vsearch_uchime3:

uchime3_denovo
--------------

| Perform chimera filtering with **uchime3_denovo** algorithm in `vsearch <https://github.com/torognes/vsearch>`_ 
| Designed for denoised amplicons. 
| uchime3_denovo can be applied also in :ref:`UNOISE3 clustering <clustering_unoise3>`

+----------------------+----------------------------------------------------------------------+
| Setting              | Tooltip                                                              |
+======================+======================================================================+
|| ``pre_cluster``     || identity percentage when performing 'pre-clustering' with           |
||                     || --cluster_size for denovo chimera filtering with --uchime_denovo    |
+----------------------+----------------------------------------------------------------------+
|| ``min_unique_size`` || minimum amount of a unique sequences in a fasta file. If value = 1, |
||                     || then no sequences are discarded after dereplication; if value = 2,  |
||                     || then sequences, which are represented only once in a given file are |
||                     || discarded; and so on                                                |
+----------------------+----------------------------------------------------------------------+
| ``denovo``           | if TRUE, then perform denovo chimera filtering with --uchime_denovo  |
+----------------------+----------------------------------------------------------------------+
|| ``reference_based`` || perform reference database based chimera filtering with             |
||                     || --uchime_ref. Select fasta formatted reference database (e.g. UNITE |
||                     || for ITS reads.                                                      |
||                     || If denovo = TRUE, then reference based chimera filtering will       |
||                     || be performed after denovo.                                          |
||                     ||                                                                     |
+----------------------+----------------------------------------------------------------------+
|| ``abundance_skew``  || the abundance skew is used to distinguish in a threeway alignment   |
||                     || which sequence is the chimera and which are the parents. The        |
||                     || assumption is that chimeras appear later in the PCR amplification   |
||                     || process and are therefore less abundant than their parents. The     |
||                     || default value is 2.0, which means that the parents should be at     |
||                     || least 2 times more abundant than their chimera. Any positive value  |
||                     || equal or greater than 1.0 can be used                               |
+----------------------+----------------------------------------------------------------------+
|| ``min_h``           || minimum score (h). Increasing this value tends to reduce the number |
||                     || of false positives and to decrease sensitivity. Values ranging from |
||                     || 0.0 to 1.0 included are accepted                                    |
+----------------------+----------------------------------------------------------------------+

.. _itsextractor:

____________________________________________________

|

`ITS Extractor <https://microbiology.se/software/itsx/>`_
==========================================================

ITSx (`Bengtsson-Palme et al. 2013 <https://doi.org/10.1111/2041-210X.12073>`_) detects **ITS regions** by searching for conserved rRNA gene fragments 
(18S, 5.8S, 28S) using profile HMMs (HMMER). 
When these boundaries are found, ITSx extracts the requested 
region(s) (e.g. **ITS1**, **ITS2**, or the full **ITS1-5.8S-ITS2**) and outputs new FASTA files. 
Parameters such as ``e-value``, ``scores``, and the required number of 
matched ``domains`` control how strict the rRNA-gene detection is (stricter settings reduce false positives but may remove divergent sequences).

If ``truncate`` is FALSE, then ITSx will identify the ITS sequences but 
does not trim the flanking regions (default ``truncate`` is TRUE).

ITSx is may be useful as it **standardizes what part of the rDNA amplicon you cluster and compare**:

- It removes conserved flanking rRNA gene segments (18S/5.8S/28S) so clustering is driven by the ITS barcode region.
- It improves comparability across taxa and studies by extracting the same region (ITS1/ITS2/full ITS) even when reads contain different amounts of flanking sequence.

You *can* cluster ITS reads with flanking regions, 
but it is often suboptimal because conserved rRNA segments can inflate sequence 
similarity and bias clustering (e.g., over-merging distinct ITS variants or producing inconsistent distances when flanking lengths differ).

.. note::

  Note that if the primer binding sites are close to the ITS region, then for better detection of the 18S, 5.8S and/or 28S regions 
  it may be beneficial to keep the primers (i.e. do not use 'CUT PRIMERS') .

| **Fastq/fasta** formatted single-end data is supported [fastq inputs will be converted to fasta].
| **Outputs** are fasta files in ``ITSx_out`` directory.

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR`` and the ``sequence files extension`` (fasta or fastq), 
  but the ``read types`` (paired-end or single-end) does not matter here (just click 'Confirm').

+----------------+-----------------------------------------------------------------------+
| Setting        | Tooltip                                                               |
+================+=======================================================================+
|| ``organisms`` || set of profiles to use for the search. Can be used to restrict the   |
||               || search to only a few organism groups types to save time, if one or   |
||               || more of the origins are not relevant to the dataset under study      |
+----------------+-----------------------------------------------------------------------+
|| ``regions``   || ITS regions to output (note that 'all' will output also full ITS     |
||               || region [ITS1-5.8S-ITS2])                                             |
+----------------+-----------------------------------------------------------------------+
|| ``partial``   || if larger than 0, ITSx will save additional FASTA-files for full and |
||               || partial ITS sequences longer than the specified value. This can be   |
||               || beneficial when some taxa have ITS regions that are too long to be   |
||               || fully covered/assembled (or reads are truncated by quality), so      |
||               || keeping partial ITS sequences helps retain those taxa in downstream  |
||               || clustering and diversity analyses. If this setting is left to 0      |
||               || (zero), it means OFF                                                 |
+----------------+-----------------------------------------------------------------------+
|| ``e-value``   || domain e-value cutoff a sequence must obtain in the HMMER-based step |
||               || to be included in the output                                         |
+----------------+-----------------------------------------------------------------------+
|| ``scores``    || domain score cutoff that a sequence must obtain in the HMMER-based   |
||               || step to be included in the output                                    |
+----------------+-----------------------------------------------------------------------+
|| ``domains``   || the minimum number of domains (different HMM gene profiles) that     |
||               || must match a sequence for it to be included in the output (detected  |
||               || as an ITS sequence). Setting the value lower than two will increase  |
||               || the number of false positives, while increasing it above two will    |
||               || decrease ITSx detection abilities on fragmentary data                |
+----------------+-----------------------------------------------------------------------+
| ``complement`` | if TRUE, ITSx checks both DNA strands for matches to HMM-profiles     |
+----------------+-----------------------------------------------------------------------+
|| ``only full`` || If TRUE, the output is limited to full-length ITS1 and ITS2 regions  |
||               || only                                                                 |
+----------------+-----------------------------------------------------------------------+
|| ``truncate``  || removes ends of ITS sequences if they are outside of the ITS region. |
||               || If FALSE, the whole input sequence is saved                          |
+----------------+-----------------------------------------------------------------------+

____________________________________________________

|

.. _clustering:

CLUSTERING
==========

Cluster sequences, generate OTUs (with :ref:`vsearch <clustering_vsearch>`), 
swarm-clusters (with :ref:`SWARM <clustering_swarm>`) 
or zOTUs (with :ref:`UNOISE3 <clustering_unoise3>`).

Clustering groups similar sequences into units that are treated as the same biological feature. 
Reads that are sufficiently similar (method-dependent) are assigned to the same cluster, 
producing a representative sequence for each cluster (fasta file) and a corresponding abundance table per sample ("OTU table").

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR`` and the ``sequence files extension`` (must be fasta/fa), 
  but the ``read types`` (paired-end or single-end) does not matter here (just click 'Confirm').

.. _clustering_vsearch:

`vsearch <https://github.com/torognes/vsearch>`_ 
------------------------------------------------

vsearch performs a similarity-threshold clustering (e.g., 97% identity threshold). 
Sequences are grouped if they meet a chosen global identity cutoff. 
Output units are "traditional" OTUs and depend strongly on the selected identity threshold and input data.

+---------------------------+----------------------------------------------------------------------+
| Tooltip                   |                                                                      |
+===========================+======================================================================+
|| ``OTU_type``             || centroid" = output centroid sequences; "consensus" = output         |
||                          || consensus sequences                                                 |
+---------------------------+----------------------------------------------------------------------+
|| ``similarity_threshold`` || define OTUs based on the sequence similarity threshold; 0.97 = 97%  |
||                          || similarity threshold                                                |
+---------------------------+----------------------------------------------------------------------+
|| ``strands``              || when comparing sequences with the cluster seed, check both strands  |
||                          || (forward and reverse complementary) or the plus strand only         |
+---------------------------+----------------------------------------------------------------------+
|| ``remove_singletons``    || if TRUE, then singleton OTUs will be discarded (OTUs with only one  |
||                          || sequence)                                                           |
+---------------------------+----------------------------------------------------------------------+
|| ``similarity_type``      || pairwise sequence identity definition                               |
||                          || `--iddef <_static/vsearch_manual_2.22.1.pdf>`_                      |
+---------------------------+----------------------------------------------------------------------+
|| ``sequence_sorting``     || size = sort the sequences by decreasing abundance; "length" = sort  |
||                          || the sequences by decreasing length (--cluster_fast); "no" = do not  |
||                          || sort sequences (--cluster_smallmem --usersort)                      |
+---------------------------+----------------------------------------------------------------------+
|| ``centroid_type``        || "similarity" = assign representative sequence to the closest (most  |
||                          || similar) centroid (distance-based greedy clustering); "abundance" = |
||                          || assign representative sequence to the most abundant centroid        |
||                          || (abundance-based greedy clustering; --sizeorder), ``max_hits``      |
||                          || should be > 1                                                       |
+---------------------------+----------------------------------------------------------------------+
|| ``max_hits``             || maximum number of hits to accept before stopping the search (should |
||                          || be > 1 for abundance-based selection of centroids [centroid type])  |
+---------------------------+----------------------------------------------------------------------+
|| ``mask``                 || mask regions in sequences using the "dust" method, or do not mask   |
||                          || ("none")                                                            |
+---------------------------+----------------------------------------------------------------------+

.. _clustering_swarm:

`SWARM <https://github.com/torognes/swarm>`_ 
---------------------------------------------

Cluster sequences using SWARM (`Mahé et al. 2021 <https://doi.org/10.1093/bioinformatics/btab493>`_), a 
robust and scalable clustering method that does not rely on an global clustering threshold.
Sequences are clustered using a distance parameter (``d``: maximum number of differences between reads (local linking threshold)), 
where clusters are resilient to input-order changes, therefore, forming stable features.

+----------------------+--------------------------------------------------------------------+
| Setting              | Tooltip                                                            |
+======================+====================================================================+
|| ``d``               || the maximum number of differences allowed between two amplicons.  |
||                     || Resolution of 1 is recommended for denoising. Higher values group |
||                     || sequences more loosely into swarm-clusters. Default = 1           |
+----------------------+--------------------------------------------------------------------+
|| ``no_otu_breaking`` || if TRUE, prevents the so-called 'OTU-breaking', which refines     |
||                     || clusters by eliminating links between amplicons with specific     |
||                     || abundance patterns. Recommended for denoising                     |
+----------------------+--------------------------------------------------------------------+
|| ``fastidious``      || if TRUE (and resolution = 1), SWARM uses fastidious mode which    |
||                     || grafts low-abundance swarms onto larger ones, reducing the number |
||                     || of small clusters. Highly recommended for resolution=1            |
+----------------------+--------------------------------------------------------------------+
|| ``boundary``        || minimum mass of a large swarm for fastidious mode. Only applies   |
||                     || when fastidious mode is enabled and resolution = 1. Default = 3   |
+----------------------+--------------------------------------------------------------------+

**Advanced options (for resolution > 1):**

+------------------------------+----------------------------------------------------------------------+
| Setting                      | Tooltip                                                              |
+==============================+======================================================================+
|| ``match_reward``            || reward for a nucleotide match in pairwise alignment. Only applies   |
||                             || when resolution > 1. Default = 5                                    |
+------------------------------+----------------------------------------------------------------------+
|| ``mismatch_penalty``        || penalty for a nucleotide mismatch in pairwise alignment. Only       |
||                             || applies when resolution > 1. Default = 4                            |
+------------------------------+----------------------------------------------------------------------+
|| ``gap_opening_penalty``     || penalty for opening a gap in pairwise alignment. Only applies when  |
||                             || resolution > 1. Default = 12                                        |
+------------------------------+----------------------------------------------------------------------+
|| ``gap_extension_penalty``   || penalty for extending a gap in pairwise alignment. Only applies     |
||                             || when resolution > 1. Default = 4                                    |
+------------------------------+----------------------------------------------------------------------+

.. _clustering_unoise3:

`UNOISE3, with vsearch <https://github.com/torognes/vsearch>`_ 
---------------------------------------------------------------

UNOISE3 is a denoising-based approach that outputs ASVs (zOTUs) that that represent sequence variants, not similarity-binned OTUs.

+---------------------+-----------------------------------------------------------------------+
| Tooltip             |                                                                       |
+=====================+=======================================================================+
|| ``strands``        || when comparing sequences with the cluster seed, check both strands   |
||                    || (forward and reverse complementary) or the plus strand only          |
+---------------------+-----------------------------------------------------------------------+
| ``minsize``         | minimum abundance of sequences for denoising                          |
+---------------------+-----------------------------------------------------------------------+
| ``remove_chimeras`` | perform chimera removal with **uchime3_denovo** algoritm              |
+---------------------+-----------------------------------------------------------------------+
|| ``unoise_alpha``   || alpha parameter to the vsearch --cluster_unoise command. default =   |
||                    || 2.0.                                                                 |
+---------------------+-----------------------------------------------------------------------+
|| ``denoise_level``  || at which level to perform denoising; global = by pooling samples,    |
||                    || individual = independently for each sample (if samples are denoised  |
||                    || individually, reducing minsize to 4 may be more reasonable for       |
||                    || higher sensitivity)                                                  |
+---------------------+-----------------------------------------------------------------------+
|| ``abskew``         || the abundance skew of chimeric sequences in comparsion with parental |
||                    || sequences (by default, parents should be at least 16 times more      |
||                    || abundant than their chimera)                                         |
+---------------------+-----------------------------------------------------------------------+
| ``maxaccepts``      | maximum number of hits to accept before stopping the search           |
+---------------------+-----------------------------------------------------------------------+
|| ``maxrejects``     || maximum number of non-matching target sequences to consider before   |
||                    || stopping the search                                                  |
+---------------------+-----------------------------------------------------------------------+
|| ``mask``           || mask regions in sequences using the "dust" method, or do not mask    |
||                    || ("none")                                                             |
+---------------------+-----------------------------------------------------------------------+

.. _assign_taxonomy:

____________________________________________________

ASSIGN TAXONOMY
===============

Implemented tools for taxonomy annotation:

.. _assign_taxonomy_blast:

`BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_ 
---------------------------------------------------

| BLAST search (`Camacho et al. 2009 <https://doi.org/10.1186/1471-2105-10-421>`_) sequences againt selected :ref:`database <databases>`. 

.. important::

 **BLAST database needs to be an unzipped fasta file in a separate folder** (fasta will be automatically converted to BLAST database files). 
 If converted BLAST database files (.ndb, .nhr, .nin, .not, .nsq, .ntf, .nto) already exist, then just SELECT **one** of those files as BLAST database in 
 'ASSIGN TAXONOMY' panel.

| Supported file format for the input data is **fasta**.
| 
| **Output** files in``taxonomy_out`` directory:
| # BLAST_1st_best_hit.txt = BLAST results for the 1st best hit in the used database.
| # BLAST_10_best_hits.txt = BLAST results for the 10 best hits in the used database.

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR`` (will be the output directory),
  but the ``sequence files extension`` and ``read type`` (single-end or paired-end) does not matter here (just click 'Next').

.. important::

  Make sure you do not have any other BLAST database files is the same directory as the database you are using.
  That is, use dedicated directory for the BLAST database.

.. note::

 BLAST values filed separator is '+'. When pasting the taxonomy results to e.g. Excel, then first denote '+' as 
 as filed separator to align the columns.


 **Check** :ref:`this section <parse_BLAST_results>` for additional parsing of the BLAST results.

+--------------------+----------------------------------------------------------------------+
| Setting            | Tooltip                                                              |
+====================+======================================================================+
|| ``database_file`` || select a database file in fasta format. Fasta format will be        |
||                   || automatically converted to BLAST database                           |
+--------------------+----------------------------------------------------------------------+
| ``fasta_file``     | select a fasta file to be used as a query for BLAST search           |
+--------------------+----------------------------------------------------------------------+
| ``task``           | BLAST search settings according to blastn or megablast               |
+--------------------+----------------------------------------------------------------------+
|| ``strands``       || query strand to search against database. Both = search also reverse |
||                   || complement                                                          |
+--------------------+----------------------------------------------------------------------+
|| ``e_value``       || a parameter that describes the number of hits one can expect to see |
||                   || by chance when searching a database of a particular size. The lower |
||                   || the e-value the more 'significant' the match is                     |
+--------------------+----------------------------------------------------------------------+
|| ``word_size``     || the size of the initial word that must be matched between the       |
||                   || database and the query sequence                                     |
+--------------------+----------------------------------------------------------------------+
| ``reward``         | reward for a match                                                   |
+--------------------+----------------------------------------------------------------------+
| ``penalty``        | penalty for a mismatch                                               |
+--------------------+----------------------------------------------------------------------+
| ``gap_open``       | cost to open a gap                                                   |
+--------------------+----------------------------------------------------------------------+
| ``gap_extend``     | cost to extend a gap                                                 |
+--------------------+----------------------------------------------------------------------+

____________________________________________________

|

.. _assign_taxonomy_rdp:

RDP classifier
---------------

Classify sequences with RDP classifier (`Wang et al. 2007 <https://doi.org/10.1128/aem.00062-07>`_) againt trained RDP database.

.. important::

 **RDP classifier database needs to be an a trained database** 
 Check section "Trained classifiers that work with MetaWorks and the RDP Classifier" from `MetaWorks <https://terrimporter.github.io/MetaWorksSite/>`_ 
 for the list of trained databases.

| 
| **Output** files in ``taxonomy_out.rdp`` directory:
| # taxonomy.txt = classifier results with bootstrap values.

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR`` (will be the output directory),
  but the ``sequence files extension`` and ``read type`` (single-end or paired-end) does not matter here (just click 'Next').

+----------------+--------------------------------------------------------------+
| Setting        | Tooltip                                                      |
+================+==============================================================+
| ``database``   | select a trained RDP classifier database                     |
+----------------+--------------------------------------------------------------+
| ``fasta_file`` | select a fasta file to be used as a query for RDP classifier |
+----------------+--------------------------------------------------------------+
| ``confidence`` | confidence threshold for assigning a taxonomic level         |
+----------------+--------------------------------------------------------------+
| ``mem``        | the amount of memory to allocate for the RDP classifier      |
+----------------+--------------------------------------------------------------+

__________________________________________________

|

.. _assign_taxonomy_sintax:

SINTAX
------

| Classify sequences with SINTAX (`Edgar 2016 <https://www.biorxiv.org/content/10.1101/074161v1>`_) againt selected :ref:`database <databases>` in fasta format.

.. important::

  Note that the database sequence headers need to be in the following format: 
  >CP002711;tax=d:Fungi,p:Ascomycota,c:Saccharomycetes,o:Saccharomycetales,
  f:Saccharomycetaceae,g:Eremothecium,s:gossypii;

  | In this format:
  | - d denotes the domain
  | - p denotes the phylum
  | - c denotes the class
  | - o denotes the order
  | - f denotes the family
  | - g denotes the genus
  | - s denotes the species

  This structured header allows SINTAX to accurately interpret the taxonomic hierarchy of each reference sequence.

| **Output** files in ``taxonomy_out.sintax`` directory:
| # taxonomy.sintax.txt = classifier results with bootstrap values.


.. note::

  To **START**, specify working directory under ``SELECT WORKDIR`` (will be the output directory),
  but the ``sequence files extension`` and ``read type`` (single-end or paired-end) does not matter here (just click 'Confirm').

+----------------+---------------------------------------------------------------------+
| Setting        | Tooltip                                                             |
+================+=====================================================================+
| ``database``   | select database file (following the format above)                   |
+----------------+---------------------------------------------------------------------+
| ``fasta_file`` | select a fasta file to be used as a query for SINTAX                |
+----------------+---------------------------------------------------------------------+
| ``cutoff``     | confidence threshold for assigning a taxonomic level                |
+----------------+---------------------------------------------------------------------+
|| ``strand``    || check both strands (forward and reverse complementary) or the plus |
||               || strand (fwd) only                                                  |
+----------------+---------------------------------------------------------------------+
| ``wordlength`` | length of k-mers for database indexing (default is 8)               |
+----------------+---------------------------------------------------------------------+


____________________________________________________

|

.. _assign_taxonomy_dada2:

`DADA2 classifier <https://github.com/benjjneb/dada2>`_ 
-------------------------------------------------------

| Classify sequences with DADA2 RDP naive Bayesian classifier (function *assignTaxonomy*) againt selected :ref:`database <databases>`.

| Supported file format for the input data is **fasta**.
| 
| **Output** files in``taxonomy_out.dada2`` directory:
| # taxonomy.txt = classifier results with bootstrap values.

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR`` (will be the output directory),
  but the ``sequence files extension`` and ``read type`` (single-end or paired-end) does not matter here (just click 'Confirm').

+---------------------+--------------------------------------------------------------------------------------------------------+
| Setting             | Tooltip                                                                                                |
+=====================+========================================================================================================+
|| ``dada2_database`` || select a reference database fasta file for taxonomy annotation.                                       |
||                    || Download DADA2-formatted reference databases `here <https://benjjneb.github.io/dada2/training.html>`_ |
+---------------------+--------------------------------------------------------------------------------------------------------+
| ``fasta_file``      | select a fasta file to be used as a query for DADA2 classifier                                         |
+---------------------+--------------------------------------------------------------------------------------------------------+
| ``minBoot``         | the minimum bootstrap confidence for assigning a taxonomic level                                       |
+---------------------+--------------------------------------------------------------------------------------------------------+
|| ``tryRC``          || the reverse-complement of each sequences will be used for classification                              |
||                    || if it is a better match to the reference sequences than the forward sequence                          |
+---------------------+--------------------------------------------------------------------------------------------------------+


____________________________________________________

|

.. _assign_taxonomy_boldigger3:

`BOLDigger3 <https://github.com/DominikBuchner/BOLDigger3>`_
------------------------------------------------------------

| Identify sequences with BOLDigger3 (`Buchner & Leese 2020 <https://doi.org/10.3897/mbmg.4.53535>`_) 
| against the `BOLD Systems v5 <https://www.boldsystems.org/>`_ online identification engine.

BOLDigger3 is an automated tool designed for DNA sequence identification through BOLD Systems v5. 
It provides high-performance processing with up to 10,000 identifications per hour (depending on settings), 
and features an intelligent top-hit selection algorithm that considers similarity thresholds at different 
taxonomic levels.

.. important::

  **No local database download required.** BOLDigger3 queries the BOLD Systems v5 online database directly.
  Because data are retrieved online for each run, total download time can remain substantial regardless of query FASTA file size.
  The tool automatically manages the identification process, including queuing requests, downloading results, 
  and selecting the best-fitting taxonomic assignment.

| Supported file format for the input data is **fasta**.
| 
| **Output** files in ``taxonomy_out.boldigger3`` directory:
| # \*_identification_result.xlsx = taxonomy assignments with detailed metadata (Excel format)
| # \*_identification_result.parquet.snappy = taxonomy assignments (Parquet format for large datasets)
| 
| Additional BOLDigger3-generated files (metadata, cache, intermediate files) are also retained in the output directory.

.. note::

  To **START**, specify working directory under ``SELECT WORKDIR`` (will be the output directory),
  but the ``sequence files extension`` and ``read type`` (single-end or paired-end) does not matter here (just click 'Confirm').

+----------------+--------------------------------------------------------------------------------------+
| Setting        | Tooltip                                                                              |
+================+======================================================================================+
| ``fasta_file`` | select a fasta file to be used as a query for BOLDigger3                             |
+----------------+--------------------------------------------------------------------------------------+
|| ``database``  || BOLD v5 database number (1-8). See database list below                              |
+----------------+--------------------------------------------------------------------------------------+
|| ``mode``      || operating mode (1-3) that determines identification speed and thoroughness.         |
||               || See operating modes below                                                           |
+----------------+--------------------------------------------------------------------------------------+
|| ``thresholds``|| similarity thresholds (space-separated) for taxonomic levels: Species Genus Family  |
||               || Order [Class]. Up to 5 values can be specified. Default: '97 95 90 85'.             |
||               || Example: '99 97' sets Species=99%, Genus=97%, remaining levels use defaults         |
+----------------+--------------------------------------------------------------------------------------+

**BOLD v5 Databases:**

| 1 = ANIMAL LIBRARY (PUBLIC)
| 2 = ANIMAL SPECIES-LEVEL LIBRARY (PUBLIC + PRIVATE)
| 3 = ANIMAL LIBRARY (PUBLIC + PRIVATE)
| 4 = VALIDATED CANADIAN ARTHROPOD LIBRARY
| 5 = PLANT LIBRARY (PUBLIC)
| 6 = FUNGI LIBRARY (PUBLIC)
| 7 = ANIMAL SECONDARY MARKERS (PUBLIC)
| 8 = VALIDATED ANIMAL RED LIST LIBRARY

**Operating Modes:**

| 1 = **Rapid Species Search** (fastest, up to 1000 sequences/batch, ~10,000 sequences/hour)
| 2 = **Genus and Species Search** (200 sequences/batch, moderate speed)
| 3 = **Exhaustive Search** (most thorough, 100 sequences/batch, slowest)

.. note::

  **Top-Hit Selection Algorithm:** BOLDigger3 uses an intelligent algorithm to select the best taxonomic 
  assignment from up to 100 hits. It applies similarity thresholds at different taxonomic levels 
  (default: 97% for species, 95% for genus, 90% for family, 85% for order), and selects the most 
  common hit that has complete taxonomic information. The algorithm also implements a flagging system 
  to highlight uncertain assignments (e.g., reverse BIN taxonomy, private data, unique hits, multiple BINs).


____________________________________________________

|

.. _databases:

Sequence databases
------------------

A *(noncomprehensive)* list of **public databases available for taxonomy annotation**:

+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+
| Database                                                                                                                                      | Description                                                |
+===============================================================================================================================================+============================================================+
| `EUKARYOME <https://eukaryome.org/>`_                                                                                                         | 18S rRNA (SSU), ITS, and 28S rRNA (LSU) for all eukaryotes |
+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+
| `UNITE <https://unite.ut.ee/repository.php>`_                                                                                                 | ITS rRNA, Fungi and all Eukaryotes                         |
+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+
| `SILVA <https://www.arb-silva.de/>`_                                                                                                          | 16S/18S (SSU), Bacteria, Archaea and Eukarya               |
+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+
| `MIDORI2 <www.reference-midori.info/index.html>`_                                                                                             | Eukaryota mitochondrial genes (including COI)              |
+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+
| `CO1 Classifier <https://github.com/terrimporter/CO1Classifier>`_                                                                             | Metazoa COI (includes outgroups)                           |
+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+
| `BOLD (distilled) database <https://boldsystems.org/data/boldistilled/>`_                                                                     | Metazoa COI (includes outgroups)                           |
+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+
| `DADA2-formatted reference databases <https://benjjneb.github.io/dada2/training.html>`_                                                       | Multiple third-party databases                             |
+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+
| `DIAT.BARCODE database <https://carrtel-collection.hub.inrae.fr/barcoding-databases/diat.barcode/pipelines-and-aligned-and-trimed-database>`_ | Diatoms rbcL/18S                                           |
+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+
| `PR2 database <https://github.com/pr2database/pr2database/releases>`_                                                                          | 18S rRNA, all Eukaryotes                                  |
+-----------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------+


____________________________________________________

.. _postprocessing:

POSTPROCESSING
==============

Post-processing tools. :ref:`See this page <postprocessingtools>`

____________________________________________________

.. _utilities:

UTILITIES
=========

Utility tools for sequence processing and manipulation.

.. _utilities_reorient:

reorient
--------

Sequences are often (if not always) in both, 5'-3' and 3'-5', orientations in the raw sequencing data sets. 
If the data still contains PCR primers that were used to generate amplicons, 
then by specifying these PCR primers, this panel will perform sequence reorientation 
of all sequences. 

**Generally, this step is not needed** when following **vsearch OTUs** or **UNOISE ASVs** pipeline, 
because both strands of the sequences can be compared prior forming OTUs (``strand=both``). 
This is automatically handled also in **NextITS** pipeline.
In the **DADA2 ASVs** pipeline, if working with mixed orientation data (seqs in 5'-3' and 3'-5' orientations), 
then select ``PAIRED-END MIXED`` mode to account for mixed orientation data. 

**Process description:** for reorienting, 
first the forward primer will be searched (using `fqgrep <https://github.com/indraniel/fqgrep>`_)  
and if detected then the read is considered as forward complementary (5'-3'). 
Then the reverse primer will be searched (using `fqgrep <https://github.com/indraniel/fqgrep>`_) 
from the same input data and if detected, then the read is considered to be in 
reverse complementary orientation (3'-5'). Latter reads will be transformed to 5'-3' 
orientation and merged with other 5'-3' reads. 
Note that for paired-end data, R1 files will be reoriented to 5'-3' 
but R2 reads will be reoriented to 3'-5' in order to merge paired-end reads.

At least one of the PCR primers must be found in the sequence. 
For example, read will be recorded if forward primer was found even 
though reverse primer was not found (and vice versa). 
**Sequence is discarded if none of the PCR primers are found.** 

Sequences that contain **multiple forward or reverse primers (multi-primer artefacts) 
are discarded** as it is highly likely that these are chimeric sequences. 
Reorienting sequences **will not remove** primer strings from the sequences. 

.. note::

 For single-end data, sequences will be reoriented also during 
 the 'cut primers' process (see below); therefore this step may be skipped
 when working with single-end data (such as data from PacBio machines OR already assembled paired-end data).

Supported file formats for paired-end input data are only **fastq**,
but also **fasta** for single-end data.
**Outputs** are fastq/fasta files in ``reoriented_out`` directory. 
Primers are **not truncated** from the sequences; this can be done using :ref:`CUT PRIMER panel <remove_primers>`

+----------------------------------+----------------------------------------------------------------------+
| Setting                          | Tooltip                                                              |
+==================================+======================================================================+
| ``mismatches``                   | allowed mismatches in the primer search                              |
+----------------------------------+----------------------------------------------------------------------+
| ``forward_primers``              | specify forward primer **(5'-3')**; IUPAC codes allowed; add up to   |
|                                  | 13 primers                                                           |
+----------------------------------+----------------------------------------------------------------------+
| ``reverse_primers``              | specify reverse primer **(3'-5')**; IUPAC codes allowed; add up to   |
|                                  | 13 primers                                                           |
+----------------------------------+----------------------------------------------------------------------+

____________________________________________________


|

.. _utilities_seqkit_stats:

seqkit stats
------------

Get sequence statistics with `seqkit stats <https://bioinf.shenwei.me/seqkit/>`_. 
Works with **fasta(.gz)/fastq(.gz)** files in the WORKING DIRECTORY. 

**Output** is the tab-delimited text file *seqkit_stats.$fileFormat.txt* with the following content:

+-----------+---------------------------+
| Statistic | Description               |
+===========+===========================+
| file      | Input file name           |
+-----------+---------------------------+
| format    | File format (FASTA/FASTQ) |
+-----------+---------------------------+
| type      | Sequence type (DNA/RNA)   |
+-----------+---------------------------+
| num_seqs  | Number of sequences       |
+-----------+---------------------------+
| sum_len   | Total sequence length     |
+-----------+---------------------------+
| min_len   | Minimum sequence length   |
+-----------+---------------------------+
| avg_len   | Average sequence length   |
+-----------+---------------------------+
| max_len   | Maximum sequence length   |
+-----------+---------------------------+

____________________________________________________

.. _utilities_self-comparison:

Self-comparison
---------------

You can run self-comparison of sequences in a fasta file to find identical or similar sequences within the same file. 
There are two methods implemented: BLAST and vsearch. This tool is useful for identifying duplicate, near-duplicate, 
or highly similar sequences within your dataset.

| **Supported file format** for input data is **fasta**.
| **Outputs** are tab-delimited text files in ``self_comparison_out`` directory.


+------------------------------------+-------------------------------------+
|              Setting               |            Description              |
+====================================+=====================================+
|          method                    |     Choose between 'vsearch' or     |
|                                    |     'blast' for sequence comparison |
+------------------------------------+-------------------------------------+
|          fasta_file                |     Select input fasta file for     |
|                                    |     self-comparison analysis        |
+------------------------------------+-------------------------------------+
|      identity_threshold            |     Minimum sequence identity       |
|                                    |     percentage to report matches    |
|                                    |     (default: 60%)                  |
+------------------------------------+-------------------------------------+
|      coverage_threshold            |     Minimum sequence coverage       |
|                                    |     percentage to report matches    |
|                                    |     (default: 60%)                  |
+------------------------------------+-------------------------------------+
|           strand                   |              both or plus           |
+------------------------------------+-------------------------------------+



**vsearch output:**


+---------+---------------------------------+
| Column  | Description                     |
+=========+=================================+
| query   | Query sequence identifier       |
+---------+---------------------------------+
| target  | Target sequence identifier      |
+---------+---------------------------------+
| id      | Sequence identity percentage    |
+---------+---------------------------------+
| alnlen  | Alignment length                |
+---------+---------------------------------+
| qcov    | Query coverage percentage       |
+---------+---------------------------------+
| tcov    | Target coverage percentage      |
+---------+---------------------------------+
| ql      | Query sequence length           |
+---------+---------------------------------+
| tl      | Target sequence length          |
+---------+---------------------------------+
| ids     | Number of identical positions   |
+---------+---------------------------------+
| mism    | Number of mismatches            |
+---------+---------------------------------+
| gaps    | Number of gap openings          |
+---------+---------------------------------+
| qilo    | Query alignment start position  |
+---------+---------------------------------+
| qihi    | Query alignment end position    |
+---------+---------------------------------+
| qstrand | Query strand orientation (+/-)  |
+---------+---------------------------------+
| tstrand | Target strand orientation (+/-) |
+---------+---------------------------------+


**BLAST output:**

+----------+----------------------------------+
| Column   | Description                      |
+==========+==================================+
| qseqid   | Query sequence identifier        |
+----------+----------------------------------+
| sseqid   | Subject sequence identifier      |
+----------+----------------------------------+
| pident   | Percentage of identical matches  |
+----------+----------------------------------+
| length   | Alignment length                 |
+----------+----------------------------------+
| mismatch | Number of mismatches             |
+----------+----------------------------------+
| gapopen  | Number of gap openings           |
+----------+----------------------------------+
| qstart   | Query alignment start position   |
+----------+----------------------------------+
| qend     | Query alignment end position     |
+----------+----------------------------------+
| sstart   | Subject alignment start position |
+----------+----------------------------------+
| send     | Subject alignment end position   |
+----------+----------------------------------+
| evalue   | Expect value                     |
+----------+----------------------------------+
| bitscore | Bit score                        |
+----------+----------------------------------+
| qlen     | Query sequence length            |
+----------+----------------------------------+
| slen     | Subject sequence length          |
+----------+----------------------------------+
| qcovs    | Query coverage per subject       |
+----------+----------------------------------+
|| qcovhsp || Query coverage per high-scoring |
||         || pair                            |
+----------+----------------------------------+
| sstrand  | Subject strand orientation       |
+----------+----------------------------------+

____________________________________________________

.. _expert_mode:

Expert-mode (PipeCraft2 console)
================================

Bioinformatic tools used by PipeCraft2 are stored on `Dockerhub <https://hub.docker.com/u/pipecraft>`_ as Docker images. 
These images can be used to launch any tool with the Docker CLI to utilize the compiled tools.
Especially useful in Windows OS, where majority of implemented modules are not compatible. 

:ref:`See list of docker images with implemented software here <dockerimages>`

Show a list of all images in your system (using e.g. **Expert-mode**):

.. code-block::

  docker images 

Download an image if required (from `Dockerhub <https://hub.docker.com/u/pipecraft>`_):

.. code-block::
  :caption: docker pull pipecraft/IMAGE:TAG
  
  docker pull pipecraft/vsearch:2.18

Delete an image

.. code-block::
  :caption: docker rmi IMAGE 

  docker rmi pipecraft/vsearch:2.18

Run docker container in your working directory to access the files. Outputs will be generated into the specified working directory.
Specify the working directory under the -v flag:

.. code-block::

  docker run -i --tty -v users/Tom/myFiles/:/Files pipecraft/vsearch:2.18

Once inside the container, move to /Files directory, which represents your working directory in the container; and run analyses

.. code-block::

  cd Files
  vsearch --help
  vsearch *--whateversettings*
      

Exit from the container:

.. code-block:: 

  exit
