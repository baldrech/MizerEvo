# MizerEvo

This is the R code for the MizerEvo model to accompany the paper:

Forestier R., Blanchard J.L., Nash K.L., Fulton E.A, Johnson C., Audzijonyte A. 2020. Interacting forces of predation and fishing affect species’ maturation size. Ecology and Evolution. In press. 

All data and code to reproduce the figures in the paper are also included.

To run the MizerEvo:

open "main.r"
run "setting things up" section (line ~20)
install packages that do not load

everything is run from the main.r script
after "setting things up" there is a "multiple runs" comment (line 42)
le loop below is an example of a run. Just change "path_to_save" to where you want to save.

to check all the parameters go to myModel.r
obs: the canibalism parameters now determine the value in the whole interaction matrix, not just the diagonal.
the different plots are in plotFunction.r 

Please contact me if you have any questions.

Romain Forestier

P.S. I am developing a formal R extension of this model to work with mizer2.0 and updates. 
Work in progress here:
https://github.com/baldrech/mizerEvolution

