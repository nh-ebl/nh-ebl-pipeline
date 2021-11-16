This is a program that was originally written in IDL that has now been ported over to python. The purpose of this program is to rotate and allign fits images to create a mosaic of a particular field interpolated to the astrometry of a given header.

To get started on using this script, see example.py as an example for how to run the program. You will need to update the config.py file to match your directories and you may want to update the fields.txt file to contain the fields that you are interested in.

Additionally this script has some dependencies. One such dependency is the IRAS catalogue which may be found here: https://www.cita.utoronto.ca/~mamd/IRIS/IrisDownload.html

You will also need the dependencies listed in the req.txt file. To install these requirements run the command :

pip3 install -r req.txt
