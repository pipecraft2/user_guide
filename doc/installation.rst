.. image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: logo

.. meta::
    :description lang=en:
        PipeCraft manual. How to install PipeCraft

==============
Installation
==============

| This is **pre-release 0.1.0**.
| Current version **does not support work in computer clusters**.

Prerequisites
-------------
The only prerequisity is `docker <https://www.docker.com/>`_.

See OS-specific (Windows, Mac, Linux) docker installation guidelines below.

.. note:: 

 Modules of PipeCraft are distributed through docker containers, which will liberate the users from the
 struggle to install/compile various software for metabarcoding data analyses.
 **Thus, all processes are run in docker containers**.
 Relevant docker container will be downloaded prior the analysis.

.. warning::

 Your OS might warn that PipeCraft is dangerous software! Please ignore the warning in this case. 

____________________________________________________

Windows
-------

Tested with **Windows 10**. Older Windows versions do not support PipeCraft GUI workflow through docker.

1. Download `docker for windows <https://www.docker.com/get-started>`_ 

2. Download PipeCraft for `Windows: pre-release 0.1.0 <https://github.com/SuvalineVana/pipecraft/releases/download/0.1.0/pipecraft.Setup.0.1.0.exe>`_

3. Install PipeCraft via the setup executable

.. note::

 Resource limits for docker are managed by Windows; 
 but You can configure limits in a .wslconfig file (see **Settings** -> **Resources** on your docker desktop app)

____________________________________________________

MacOS
------

There are currently unknown issues with the permissions on macOS to run the underlying tasks for the processes. Fix is on the way.

____________________________________________________

Linux
------

Tested with **Ubuntu 20.04** and **Mint 20.1**. Older OS versions might not support PipeCraft GUI workflow through docker.

1. Install docker; `follow the guidelines under appropriate Linux distribution <https://docs.docker.com/engine/install/>`_

2. If you are a non-root user complete these `post-install steps <https://docs.docker.com/engine/install/linux-postinstall/>`_

3. Download PipeCraft for `Linux: pre-release 0.1.0 <https://github.com/SuvalineVana/pipecraft/releases/download/0.1.0/pipecraft-0.1.0.AppImage>`_

4. Make AppImage executable as program: 1) right click the AppImage; 2) navigate under **Permissions** tab; 3) tick the **Allow executing file as program** 

5. Run PipeCraft via the AppImage file

.. note::

 On Linux, docker can use all available host resources.

____________________________________________________

Update PipeCraft
----------------

To avaoid any potential software conflicts, all docker images of previous PipeCraft version should be removed. 
See :ref:`removing docker images <removedockerimages>` section.

.. note::

 | Currently available versions:
 | :ref:`pre-release 0.1.0 <0.1.0>`

____________________________________________________

.. _uninstalling:

Uninstalling
------------

| **Windows**: uninstall pipecraft via control panel
| **MacOS**: eject and delete the DMG file
| **Linux**: simply delete AppImage file

____________________________________________________

.. _removedockerimages:

Removing docker images
----------------------

| On **MacOS** and **Windows**: docker images and container can be easily managed from the docker dashboard. For more info visit https://docs.docker.com/desktop/dashboard/

| On **Linux** machines: containers and images are managed via the docker cli commands:
| rmi https://docs.docker.com/engine/reference/commandline/rmi/ 
| rm  https://docs.docker.com/engine/reference/commandline/rm/ 
| ``sudo docker images`` --> to see which docker images exist

