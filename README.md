# MATLAB-LOS-algorithm-to-C
This is a generated C-code version of several complex MATLAB functions

Similar to ![my previous MATLAB project](https://github.com/Antiochian/infrared-sounding), this repo will most likely be extremely boring to anyone who isn't a physicist with a highly niche interest in thermal modelling.

Its basically just a C-code version of a bundle of MATLAB functions based around raytracing and Line of Sight (LOS) algorithms for use in a sophisticated 3D thermal model I am working on as part of a research team. Compiling this C-code and linking it to MATLAB as it runs has resulted in order-of-magnitude computational speedup, which is very nice to see.

Important Note:
------------

There are lots of very clever algorithms in here that I wish I could take credit for. However, in all cases except for one I merely translated or tweaked existing algorithms to get them MEX-compatible, and the credit for these algorithms goes to their original authors (commented in the code).

The one exception to this is that I did in fact implement my own bilinear interpolation algorithm inside the mexable_LOS.c file, if you are curious.

