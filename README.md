# xrd
Generate XRD patterns from a vasp POSCAR file

Two versions: *_ser_fast (serial code) and *_par (multithreaded).

Ironically, the serial version is faster than the multithreaded version at the moment.  I'll fix this later - it could prob. be done quickly by editing the main triple loop's index starting points to match those of the serial version, but I have not tested this yet.  

To use, link the code to the database (AFF.csv, through the -a command line option) and to a POSCAR file (through the -p command line option).   

The code outputs a list of (hkl), their bragg angle, and the normalized intensity.  To compare two+ different XRD datas, you should 'un-normalize' the intensies (i.e. change the normalization factor from 'norm' to '1.0f' at the bottom of the main .cpp file).  
It will also output a csv file made for easy graphing, which only includes the bragg angle and intensities, which for this file have had a gaussian smoothing profile applied to them.  

To match with experimental work, find the type of radiation used in the experimental paper and use it's wavelength in this code (wavelength input = angstroms).  
There is also the option to adjust the temperature factor, but it is only applied as a constant for now.  

There are a lot of command line options - I recommend just reading the beginning of the main .cpp file for instructions.  

Resources used:
Many equations are from 'SPREADSHEET SIMULATION OF X-RAY POWDER DIFFRACTION' - H. W. G. SPRAGET 
and
https://en.wikipedia.org/wiki/Atomic_form_factor#cite_note-2


Atomic form factors taken directly from 
http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
