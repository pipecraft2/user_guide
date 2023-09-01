.. |PipeCraft2_logo| image:: _static/PipeCraft2_icon_v2.png
  :width: 100
  :alt: Alternative text

|PipeCraft2_logo|
  `github <https://github.com/SuvalineVana/pipecraft>`_

==============
For Developers
==============

Prerequisites
-------------

| Docker Desktop
| https://www.docker.com/products/docker-desktop/


| NodeJS 14 
| https://nodejs.org/download/release/latest-v14.x/

| Yarn

.. code-block::

  npm install --global yarn

| Git
| (https://git-scm.com/downloads)

____________________________________________________

Setting up the environment
--------------------------

| Clone the repository:

.. code-block::

  git clone https://github.com/pipecraft2/pipecraft

| Install and run:

.. code-block::

  cd pipecraft
  yarn run install_pipe
  yarn electron:serve 

____________________________________________________

Developer tools
---------------

| The apps front-end desing, data storage and navigation are built using the Vue framework and its plugins. To build new or to modify exiting components will reiqure a some proficiency in JavaScript and the Vue framework.
| Other important tools include electron which is mainly used for interacting with the file system and dockerode which is used for controlling docker.

| Electron 13:
| https://www.electronjs.org/
| Vue 2: 
| https://v2.vuejs.org/v2/guide/
| Vuetify 2:
| https://v2.vuetifyjs.com/en/introduction/why-vuetify/#what-is-vuetify3f
| Vue Router 3:
| https://v3.router.vuejs.org/guide/
| Vuex 3:
| https://v3.vuex.vuejs.org/#what-is-vuex
| dockerode:
| https://github.com/apocas/dockerode


.. note::

  All of these tools are automatically installed during the install_pipe command. 


____________________________________________________

.. hide

    TROUBLESHOOTING build
    ---------------------

    - pages dissapear: too long header underline. And/or no blank line at the end of the page. 

    Project structure
    -----------------

  
    | pipecraft/image_development: 
    | This folder contains dockerfiles which are instructions for building docker images for pipecraft. Check out the docker_commads file to view instructions on how to build, edit, run and publish docker images.

  
    | pipecraft/src/components/: 
    | This folder contains vue componenets for every input field available in pipecraft, some additional navigational components and a Run components for controlling docker.

  
    | pipecraft/src/pipecraft-core/service_scripts: 
    | This folder hosts the core scipts that are executed during workflows (in docker containers).

  
    | pipecraft/src/router/index.js: 
    | This file cotains instructions for routing and navigation (route names and according components).

  
    | pipecraft/src/store/index.js: 
    | This is an extensive storage file that is accessible by all components. The store holds static workflow data which is used for rendering front-end components, data for application state tracking and key functions for setting up workflow execution.

  
    | pipecraft/src/views:  
    | This folder contains components used by the router, these components are displayed in the center viewport of the app and often themselves use many components from the src/components folder.

  
    | pipecraft/src/App.vue: 
    | This file sets the main layout for pipecraft (navigation panels on both sides and the router-view in the middle).

  
    | pipecraft/src/background.js: 
    | This file cotains settings for app start-up, update, shutdown and window size. (These are mostly electron parameters).


.. hide

    Running a workflow 
    ------------------

    To execute a workflow Pipecraft will run multiple docker containers one-by-one (one container for every step in the workflow). The functions and controls for running a workflow are located in the src/componenets/Run.vue file.
    Running a workflow start with pre-run checks such as making sure docker-desktop is running, if proper inputs files were selected and if all mandatory inputs were filled in (check out store/index.js getters: customWorkflowReady and selectedStepsReady).
    Once these checks are complete we move on to the runWorkFlow(in Run.vue) function. The core of this function is a simple for loop which runs a docker container for each step in the workflow, this core function in acompanied by many others which account for important
    setup for the workflow such as: attaching files and folder to the container, setting user inputs as environmental variables in the container, tracking execution time, logging container outputs and setting up for the next step.
