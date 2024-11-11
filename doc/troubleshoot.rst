.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 50
  :alt: Alternative text
  :target: https://github.com/pipecraft2/user_guide

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

.. |debug| image:: _static/debug.png
  :width: 100
  :alt: Alternative text

.. |DADA2_read_identifiers| image:: _static/troubleshoot/DADA2_read_identifiers.png
  :width: 250
  :alt: Alternative text

=================================
Troubleshooting |PipeCraft2_logo| 
=================================

This page is developing based on the user feedback.

____________________________________________________

Debugging mode
==============

Turn on '**debugging mode**' (bottom-right button) to keep temporary (log) files for identifying the cause of the error

|debug|

____________________________________________________


General
=======

.. error::

 Conflict. The container name XXX is already in use by container "XXX".
 You have to remove (or rename) that container to be able to reuse that name.

**Reason**: Process stopped unexpectedly and docker container was not closed.

**Fix**: Remove the docker container (not image!) that is causing the conflict

____________________________________________________

.. error::

 No files in the output folder, but PipeCraft said "Workflow finished".

**Possible reason**: Computer's memory (RAM) is full, and process was killed. Cannot finish the analyses with those local resources. 

**Possible fix**: In Windows, try to increase the RAM size accessible to Docker (see :ref:`here <increase_RAM>`).
  Check if there was a README.txt output and read that. Please :ref:`report <contact>` unexpexted errors. 

____________________________________________________

.. error::

 No OTU_table.txt with version v0.1.4

**Reason**: known bug.

**Fix**: Fixed the bug. Reinstall PipeCraft v0.1.4 (or higher)

____________________________________________________

.. error::
  
  "ERROR]: cannot find files with specified extension"

**Reason**: wrongly specified working directory or extension; OR issues with external hard drives in Windows.

**Fix**: Double-check the specified directory and extention; OR restart Windows.


____________________________________________________

.. error::
  "Workflow stopped"

 |workflow_stopped|

**Possible reason**: Computer's memory (RAM) is full, and process was killed. Cannot finish the analyses with those local resources. 

**Possible fix**: In Windows, try to increase the RAM size accessible to Docker (see :ref:`here <increase_RAM>`).

____________________________________________________

.. error::

 Error in filterAndTrim. Every input file must have a corresponding output file.

  |DADA2_read_identifiers|

**Possible reason**: wrong read identifiers for ``read R1`` and ``read R2`` in QUALITY FILTERING panel. 

**Fix**: Check the input fastq file names and edit the identifiers. 
Specify identifyer string that is common for all R1 reads (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. 
Note that backslash is only needed to escape dot regex; e.g. when all R1 files have '_R1' string, then enter '_R1'.). 

____________________________________________________

.. error::

  "Error rates could not be estimated (this is usually because of very few reads). Error in getErrors(err, enforce = TRUE) : Error matrix is null."

  |learnErrors_fewReads|

**Possible reason**: Too small data set; samples contain too few reads for DADA2 denoising.

**Fix**: use OTU workflow.
