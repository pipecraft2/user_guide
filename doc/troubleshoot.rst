.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text

.. |learnErrors| image:: _static/troubleshoot/learnErrors.png
  :width: 250
  :alt: Alternative text

.. |dimnames| image:: _static/troubleshoot/dimnames.png
  :width: 250
  :alt: Alternative text

.. |learnErrors_fewReads| image:: _static/troubleshoot/learnErrors_fewReads.png
  :width: 250
  :alt: Alternative text

.. |workflow_stopped| image:: _static/troubleshoot/workflow_stopped.png
  :width: 250
  :alt: Alternative text

.. |DADA2_read_identifiers| image:: _static/troubleshoot/DADA2_read_identifiers.png
  :width: 250
  :alt: Alternative text

  
 
|PipeCraft2_logo|
  `github <https://github.com/SuvalineVana/pipecraft>`_

================
Troubleshooting
================

This page is developing based on the user feedback.

____________________________________________________


General
=======

.. error::

 Conflict. The container name XXX is already in use by container "XXX".
 You have to remove (or rename) that container to be able to reuse that name.

**Reason**: Process stopped unexpectedly and docker container was not closed.

**Fix**: Remove the docker container (not image!) that is causing the conflict


.. error::

 No files in the output folder, but PipeCraft said "Workflow finished".

**Reason**: ?

**Fix**: Check if there was a README.txt output and read that. Please :ref:`report <contact>` unexpexted errors. 



ASVs workflow
==============

.. error::
  "Workflow stopped"

 |workflow_stopped|

**Possible reason**: Computer's memory is full, cannot finish the analyses.

**Fix**: Analyse fewer number of samples or increase RAM size.

____________________________________________________

.. error::

 "Error in derepFastq(fls[[i]], qualityType = qualityType) : Not all provided files exist. Calls: learnErrors -> derepFastq. Execution halted"

 |learnErrors| 

**Possible reason**: Some samples have completely discarded by quality filtering process. 

**Fix**: Examine **seq_count_summary.txt** file in ``qualFiltered_out`` folder and discard samples, which had 0 quality filtered sequences (poor quality samples). Or edit the quality filtering settings.

____________________________________________________

.. error::

 Error in filterAndTrim. Every input file must have a corresponding output file.

  |DADA2_read_identifiers|

**Possible reason**: wrong read identifiers for ``read R1`` and ``read R2`` in QUALITY FILTERING panel. 

**Fix**: Check the input fastq file names and edit the identifiers. 
Specify identifyer string that is common for all R1 reads (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. 
Note that backslash is only needed to escape dot regex; e.g. when all R1 files have '_R1' string, then enter '_R1'.). When demultiplexing data in during ASV (DADA2) workflow, then specify as '\\.R1'
____________________________________________________

.. error::

  "Error rates could not be estimated (this is usually because of very few reads). Error in getErrors(err, enforce = TRUE) : Error matrix is null."

  |learnErrors_fewReads|

**Possible reason**: Too small data set; samples contain too few reads for DADA2 denoising.

**Fix**: use OTU workflow.



