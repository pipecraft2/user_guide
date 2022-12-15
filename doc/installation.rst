.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text

.. |resources| image:: _static/resources1.png
  :width: 600
  :alt: Alternative text

.. |openanyway| image:: _static/openanyway.png
  :width: 400
  :alt: Alternative text

.. |mac_docker_share| image:: _static/Mac_docker_share.png
  :width: 400
  :alt: Alternative text
  

.. meta::
    :description lang=en:
        PipeCraft manual. How to install PipeCraft


|PipeCraft2_logo|
  `github <https://github.com/pipecraft2/pipecraft>`_

==============
Installation
==============

| Current :ref:`versions <releases>` do not work on High Performance Computing (**HPC**) clusters **yet**.

____________________________________________________

.. contents:: Contents
   :depth: 3

____________________________________________________

Prerequisites
-------------
The only prerequisite is `Docker <https://www.docker.com/>`_.

See OS-specific (Windows, Mac, Linux) docker installation guidelines below.

.. note:: 

 Modules of PipeCraft2 are distributed through Docker containers, which will liberate the users from the
 struggle to install/compile various software for metabarcoding data analyses.
 **Thus, all processes are run in Docker containers**.
 Relevant Docker container will be automatically downloaded prior the analysis.

.. warning::

 Your OS might warn that PipeCraft2 is dangerous software! Please ignore the warning in this case. 

____________________________________________________

Windows
-------

PipeCraft2 was tested on **Windows 10** and **Windows 11**. Older Windows versions do not support PipeCraft GUI workflow through Docker.

1. Download `Docker for windows <https://www.docker.com/get-started>`_ 

2. Download PipeCraft for `Windows: v0.1.4 <https://github.com/pipecraft2/pipecraft/releases/download/v0.1.4/pipecraft-Setup-0.1.4.exe>`_

3. Install PipeCraft via the setup executable

.. warning::

  In Windows, please keep you working directory path as short as possible. Maximum path length in Windows is 260 characters. 
  PipeCraft may not be able to work with files, that are buried "deep inside" (i.e. the path is too long).


.. note::

 Resource limits for Docker are managed by Windows; 
 but you can configure limits in a **.wslconfig** file (see **Settings** -> **Resources** on your Docker desktop app)

____________________________________________________

MacOS
------

PipeCraft is supported on macOS 10.15+. Older OS versions might not support PipeCraft GUI workflow through Docker.

1. Check your Mac chip (Apple or Intel) and download `Docker for Mac <https://www.docker.com/get-started>`_

2. Download PipeCraft for `Mac: v0.1.4 <https://github.com/pipecraft2/pipecraft/releases/download/v0.1.4/pipecraft-0.1.4.pkg>`_

3. Install PipeCraft via **pkg** file

4. Currently, this app might be identified as app from an unidentified developer. Grant an exception for a blocked app by clicking the "**Open Anyway**" button in the General panel of **Security & Privacy** preferences

 |openanyway|

5. Open **Docker dashboard**: Settings -> Resources -> File Sharing; and add the directory where **pipecraft.app** was installed (it is usually /Appications)

 |mac_docker_share|

.. note::

 Manage Docker resource limits in the Docker dashboard:
 |resources|
 
____________________________________________________

Linux
------

PipeCraft was tested with **Ubuntu 20.04** and **Mint 20.1**. Older OS versions might not support PipeCraft GUI workflow through Docker.

1. Install Docker; `follow the guidelines under appropriate Linux distribution <https://docs.docker.com/engine/install/>`_

2. If you are a non-root user complete these `post-install steps <https://docs.docker.com/engine/install/linux-postinstall/>`_

3. Download PipeCraft for `Linux: v0.1.4 <https://github.com/pipecraft2/pipecraft/releases/download/v0.1.4/pipecraft_v0.1.4.deb>`_

4. Right click on the pipecraft_*.deb file and "Open With GDebi Package Installer" (Install Package) or ``sudo dpkg -i path_to_deb_file``
   
.. note::

   When you encounter ERROR during installation, then uninstall the previous version of PipeCraft ``sudo dpkg --remove pipecraft-v0.1.3``

5. Run PipeCraft. If PipeCraft shortcut does not appear on the Desktop, then search the app and generate shortcut manually (installed in */opt/pipecraft* directory)

.. note::

 On Linux, Docker can use all available host resources.

____________________________________________________

Updating PipeCraft2
--------------------

To avaoid any potential software conflicts from PipeCraft2 **v0.1.1 to v0.1.4**, all Docker images of older PipeCraft2 version should be removed. 
*Starting from v0.1.5 --> if docker container is updated, it will get a new tag for new PipeCraft2 version*

See :ref:`removing docker images <removedockerimages>` section.

.. note::

 | Currently available versions :ref:`HERE <releases>`

____________________________________________________

.. _uninstalling:

Uninstalling
------------

| **Windows**: uninstall PipeCraft via control panel
| **MacOS**: Move pipecraft.app to Bin
| **Linux**: remove pipecraft via Software Manager/Software Centre or via terminal ``sudo dpkg --remove pipecraft``

____________________________________________________

.. _removedockerimages:

Removing Docker images
----------------------

| On **MacOS** and **Windows**: Docker images and container can be easily managed from the Docker dashboard. For more info visit https://docs.docker.com/desktop/dashboard/
| See **command-line** based way below.

.. |purge_docker_Win| image:: _static/purge_docker_Win.png
  :width: 500
  :alt: Alternative text

|purge_docker_Win|

| 
| On **Linux** machines: containers and images are managed via the Docker cli commands (https://docs.docker.com/engine/reference/commandline/rmi/):
| ``sudo docker images`` --> to see which docker images exist
| ``sudo docker rmi IMAGE_ID_here`` --> to delete selected image

or

| ``sudo docker system prune -a`` --> to delete all unused containers, networks, images 

