Latent Evolutionary Signatures: A General Framework for Analyzing Music and Cultural Evolution

© 2024 Jonathan Warrell, Leonidas Salichos, Michael Ganczc, Mark B. Gerstein

Program in Computational Biology and Bioinformatics, Yale University, New Haven, CT 06520, USA. 

Department of Molecular Biophysics and Biochemistry, Yale University, New Haven, CT 06520, USA. 

Department of Music, Yale University, New Haven, CT 06520, USA. 

Department of Computer Science, Yale University, New Haven, CT 06520, USA.

Department of Biological and Chemical Sciences, New York Institute of Technology, New York, NY 10023, USA.


This repo contains optimal harmony and form+harmony based models for analyzing popular music as an evolutionary structure. Our model architecture is based on a traditional VAE, but with an energy-based prior that penalizes a measure of ‘evolutionary distance’, in this case informed by temporal distance across song release dates (see schematic below), in the latent space. The key output of each model is a set of latent ‘evolutionary signatures’, or characteristic distributions of chord/form k-mers that can be used to predict the date and genre of each song. We use the McGill Billboard corpus of popular song annotations as our database.

![musevo schematic](https://github.com/gersteinlab/Musevo/assets/98661752/a50e38f1-1ed3-4436-81fc-153ac2c63330)

The subdirectory ‘harmonic_formal’ contains the code for training the model based both on chord progressions and formal features. The subdirectory ‘km4’ contains the code for training the model based on chord progressions of length 4.

In addition to the optimal models, we include some other configurations that we tested, including variants on the coarse-graining of formal units (models suffixed with ‘_A’, ‘_B’, and ‘_None’, referring to projection matrices A and B, where A retains formal categories that comprise the first 99% of the data, and collapses all others into a category of ‘other’; and B bins together semantically similar formal units), and different means of normalizing formal feature vector ‘X_struct,’ which in its raw form contains the counts of each formal category (models suffixed with ‘_binarized’ or ‘_zscore’ utilize these normalizations. These are contained within the ‘models’ folder of each subdirectory. Additional items of interest, including code for processing the McGill Billboard dataset, are included in the ‘supplemental’ folder. The ‘McGill-Billboard’ folder contains all raw data.

All code is currently written for MatLab.

