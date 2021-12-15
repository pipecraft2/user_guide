.. image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: logo

.. meta::
    :description lang=en:
        PipeCraft manual. How to install PipeCraft

==============
Installation
==============

Prerequisites
-------------
The only prerequisity is `docker <https://www.docker.com/>`_.

See platform-specific (Windows, Mac, Linux) docker installation guidelines below.

*Modules of PipeCraft are distributed through docker containers, which will liberate the users from the
struggle to install/compile various software for metabarcoding data analyses.*

.. note::

 Your OS might warn that this is dangerous software! Please ignore the warning in this case. 

____________________________________________________

Windows
-------
1. `Download docker for windows <https://www.docker.com/get-started>`_ 

2. Download PipeCraft for `Windows: version 0.1.0 <https://github.com/SuvalineVana/pipecraft/releases/download/0.1.0/pipecraft.Setup.0.1.0.exe>`_

3. Install PipeCraft via the setup executable

____________________________________________________

MacOS
------
1. `Download docker for Mac <https://www.docker.com/get-started>`_

2. Download PipeCraft for `Mac: version 0.1.0 <https://github.com/SuvalineVana/pipecraft/releases/download/0.1.0/pipecraft-0.1.0.dmg>`_

3. Install PipeCraft via DMG file

____________________________________________________

Linux
------

1. `Install docker; follow the guidelines under appropriate Linux distribution <https://docs.docker.com/engine/install/>`_

2. If you are a non-root user complete these `post-install steps <https://docs.docker.com/engine/install/linux-postinstall/>`_

3. Download PipeCraft for `Linux: version 0.1.0 <https://github.com/SuvalineVana/pipecraft/releases/download/0.1.0/pipecraft-0.1.0.AppImage>`_

4. Run PipeCraft via the AppImage file 

____________________________________________________

Uninstalling
------------

| **Windows**: uninstall pipecraft via control panel
| **MacOS**: eject and delete the DMG file
| **Linux**: simply delete AppImage file

____________________________________________________

Removing docker images
----------------------

| On **MacOS** and **Windows**: docker images and container can be easily managed from the docker dashboard. For more info visit https://docs.docker.com/desktop/dashboard/

| On **Linux** machines: containers and images are managed via the docker cli commands:
| rmi https://docs.docker.com/engine/reference/commandline/rmi/ 
| rm  https://docs.docker.com/engine/reference/commandline/rm/ 

