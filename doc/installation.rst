.. image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: logo

.. |resources| image:: _static/troubleshoot/resources1.png
  :width: 600
  :alt: Alternative text

.. |openanyway| image:: _static/troubleshoot/openanyway.png
  :width: 400
  :alt: Alternative text
  

.. meta::
    :description lang=en:
        PipeCraft manual. How to install PipeCraft

==============
Installation
==============

| This is **pre-release 0.1.0**.
| Current version **does not work on High Performance Computing (HPC) clusters yet**.

Prerequisites
-------------
The only prerequisite is `Docker <https://www.docker.com/>`_.

See OS-specific (Windows, Mac, Linux) docker installation guidelines below.

.. note:: 

 Modules of PipeCraft are distributed through Docker containers, which will liberate the users from the
 struggle to install/compile various software for metabarcoding data analyses.
 **Thus, all processes are run in Docker containers**.
 Relevant Docker container will be downloaded prior the analysis.

.. warning::

 Your OS might warn that PipeCraft is dangerous software! Please ignore the warning in this case. 

____________________________________________________

Windows
-------

PipeCraft was tested on **Windows 10**. Older Windows versions do not support PipeCraft GUI workflow through Docker.

1. Download `Docker for windows <https://www.docker.com/get-started>`_ 

2. Download PipeCraft for `Windows: pre-release 0.1.0 <https://github.com/SuvalineVana/pipecraft/releases/download/0.1.0/pipecraft.Setup.0.1.0.exe>`_

3. Install PipeCraft via the setup executable

.. note::

 Resource limits for Docker are managed by Windows; 
 but you can configure limits in a **.wslconfig** file (see **Settings** -> **Resources** on your Docker desktop app)

____________________________________________________

MacOS
------

PipeCraft is supported on macOS 10.15+. Older OS versions might not support PipeCraft GUI workflow through Docker.

1. Check your Mac chip (Apple or Intel) and download `Docker for Mac <https://www.docker.com/get-started>`_

2. Download PipeCraft for `Mac: pre-release 0.1.0 <https://github.com/SuvalineVana/pipecraft/releases/download/0.1.0/pipecraft-0.1.0.pkg>`_

3. Install PipeCraft via **pkg** file

4. Currently, this app might be identified as app from an unidentified developer. Grant an exception for a blocked app by clicking the "**Open Anyway**" button in the General panel of **Security & Privacy** preferences

|openanyway|

.. note::

 Manage Docker resource limits in the Docker dashboard:
 |resources|
 
____________________________________________________

Linux
------

PipeCraft was tested with **Ubuntu 20.04** and **Mint 20.1**. Older OS versions might not support PipeCraft GUI workflow through Docker.

1. Install Docker; `follow the guidelines under appropriate Linux distribution <https://docs.docker.com/engine/install/>`_

2. If you are a non-root user complete these `post-install steps <https://docs.docker.com/engine/install/linux-postinstall/>`_

3. Download PipeCraft for `Linux: pre-release 0.1.0 <https://github.com/SuvalineVana/pipecraft/releases/download/0.1.0/pipecraft-0.1.0.AppImage>`_

4. Make AppImage executable as program: 1) right click the AppImage; 2) navigate under **Permissions** tab; 3) tick the **Allow executing file as program** 

5. Run PipeCraft via the AppImage file

.. note::

 On Linux, Docker can use all available host resources.

____________________________________________________

Update PipeCraft
----------------

To avaoid any potential software conflicts, all Docker images of previous PipeCraft version should be removed. 
See :ref:`removing docker images <removedockerimages>` section.

.. note::

 | Currently available versions:
 | :ref:`pre-release 0.1.0 <0.1.0>`

____________________________________________________

.. _uninstalling:

Uninstalling
------------

| **Windows**: uninstall PipeCraft via control panel
| **MacOS**: eject and delete the DMG file
| **Linux**: simply delete AppImage file

____________________________________________________

.. _removedockerimages:

Removing Docker images
----------------------

| On **MacOS** and **Windows**: Docker images and container can be easily managed from the Docker dashboard. For more info visit https://docs.docker.com/desktop/dashboard/

| On **Linux** machines: containers and images are managed via the Docker cli commands:
| ``rmi https://docs.docker.com/engine/reference/commandline/rmi/``
| ``rm  https://docs.docker.com/engine/reference/commandline/rm/``
| ``sudo docker images`` --> to see which docker images exist

