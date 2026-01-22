# PipeCraft2

This is repo for PipeCraft2 user guide
https://pipecraft2-manual.rtfd.io/

PipeCraft2 development repo -> https://github.com/pipecraft2/pipecraft

______________________________________________________


# For developers
This [documentation](https://bioscanflow.readthedocs.io/en/latest/index.html#) is built throuh [Read the Docs](https://about.readthedocs.com), an open-sourced software documentation hosting platform. 

**Requirements**: git, python, sphinx (documentation generator).

**1. Check if git is installed** (through command line)
```bash
git --version
```

If git is not installed, then install it
For Windows: https://git-scm.com/install/windows 
For Linux: sudo apt install git-all
For MacOS: git --version (it will prompt the install)

**2. Download github repository**
```bash
cd C:/Users/user/Desktop  # go to directory where you want to place the github repo
git clone https://github.com/anslan/BGE_biomonitoring_wf.git

cd C:/Users/user/Desktop/BGE_biomonitoring_wf/docs  # go to downloaded github repo
git checkout develop # Switch to the develop branch!
```

**3. Install sphinx and other requirements:**
   (make sure python is installed)

   ```bash
   cd C:/Users/user/Desktop/BGE_biomonitoring_wf/docs  # go to downloaded github repo
   python -m pip install -U sphinx
   python -m pip install -r requirements.txt


   ### If you have problems installing as an admin then try:
   cd C:/Users/user/Desktop/BGE_biomonitoring_wf/docs  # go to github repo
   python -m pip install --user sphinx
   python -m pip install --user -r requirements.txt
   ```

**4. Build local page for testing**
*(in BGE_biomonitoring_wf/docs)*

For Windows:
```bash
# in docs dir
.\make.bat html
```

If you installed sphinx with `--user` and `make.bat` doesn't work, use:

```bash
# in docs dir
python -m sphinx -M html . _build
```

For Linux and macOS:

```bash
# in docs dir
make html
```

Open any html file in the "BGE_biomonitoring_wf/docs/_build/html" directory to check the page build.

**6.** Edit the *.rst files as needed. Then build locally (point 4 above) to check the page before pushing. 

If changes do not appear, **then remove the "_build" folder and build again:**

```bash
# in docs dir
rm -r _build
.\make.bat html
```

**7.** When edis are done, then push changes to github
```bash
git add .    # adds all changes 
git commit -m "describe my edits" # add a brief message what was changed
git push     # push changes to github
```

After the "git push" the webpage will automatically update itself.
