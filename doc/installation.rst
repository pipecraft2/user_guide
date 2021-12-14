.. image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: logo

=============
Installation
=============


Prerequisites
-------------
The only prerequisity is `docker <https://www.docker.com/>`_.

See platform-specific (Windows, Mac, Linux) docker installation guidelines below.

*Modules of PipeCraft are distributed through docker containers, which will liberate the users from the
struggle to install/compile various software for metabarcoding data analyses.*

Windows
-------
1. Download docker for windows `here <https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe?utm_source=docker&amp;utm_medium=webreferral&amp;utm_campaign=dd-smartbutton&amp;utm_location=module>`_ (direct download link form https://www.docker.com/get-started) 

2. Download PipeCraft for Windows (`link <https://github.com/SuvalineVana/pipecraft-vue/releases/tag/0.2.0-beta>`_)

3. Install PipeCraft via the setup executable

MacOS
------
1. Download docker for Mac `here <https://www.docker.com/get-started>`_

2. Download PipeCraft for Mac (link)

3. Install PipeCraft via DMG file

Linux
-----
1. Install docker; follow the guidelines under appropriate Linux distribution `here <https://docs.docker.com/engine/install/>`_

2. If you are a non-root user complete these post-install steps https://docs.docker.com/engine/install/linux-postinstall/

3. Download PipeCraft for Linux (link)

4. Run PipeCraft via the AppImage file 


Uninstalling
------------
Windows, uninstall pipecraft via control panel
MacOS, eject and delete the DMG file
Linux, simply delete AppImage file



Removing docker images
----------------------
On MacOS and Windows docker images and container can be easily managed from the docker dashboard. For more info visit https://docs.docker.com/desktop/dashboard/
On Linux machines containers and images are managed via the docker cli commands:
rmi https://docs.docker.com/engine/reference/commandline/rmi/ 
rm  https://docs.docker.com/engine/reference/commandline/rm/ 

