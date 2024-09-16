# SDFixer

WORK IN PROGRESS

The aim is to make a tool that can easily browse thorugh an SDfile for identifying molecules with issues, e.g. unsanitizable, and provide the user with options to fix them (either automated fixes or manual editing).


## Roadmap
Create a module and package with installer and cli command

Possibility to use embedded rdEditor to fix errors manually.

Make upstream changes to rdEditor:
- Better support for embedded usage
- Visualization and reporting of chemistry problems on the molecule


Load/save of SMILES and SDFiles

Buttons to skip to next problematic molecule (both with problems, and fully unsanitizable).

Functions with standard fixes:
- charge on quarternary and uncharged nitrogen
- explicit H on aromatic 5 membered rings containing nitrogen (But how to choose which one?)


## Install

Currently only possible to run directly from directory
```
cd sdfixer
python sdfixer.py yoursdfile.sdf
```
