Saucerman model in Matlab and CellML
====================================

In order to compare the CellML version of the Saucerman model with the original Matlab code, I have generated data from them both, and wrote a script to generate plots of the data for visual comparison. This will hopefully allow future revalidation to be less painfull than it was for me.

The plot generation depeds on [Tioga](http://www.kitp.ucsb.edu/~paxton/tioga.html), and only works under Ruby 1.8. Only some of the variables can actually be compared, as the Matlab version does not record all relevant variables.
