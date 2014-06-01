electrodeRegister
=================

This code registers electrode positions to head surfaces.

For example of how to use see: exampleScript.m


Explanation:
Simply relying on 3 fiducial measurements to align electrodes to a scalp results in errors because any mismatch between defining fiducials in the electrode space and the mri space results in an a registration error. 
 
But we always have multiple electrode locations that should follow the scalp shape.  These can be used to improve the registration, but because electrodes are fairly thick they are usually digitzed at their top which is some millimeters away from the scalp surface. For this reason simply minimizing the distance to the scalp results in a misfit having all the electrodes move towards the fit by the electrode height.

This function solves this problem in two ways:

1) Electrodes are fit to the scalp with a cost function that allows electrodes to all have a consistent distance from the scalp without penalty

2) Optionally, this function utilizes additional registration points from anywhere on the scalp. These points are fit using a cost function that puts them on the scalp surface. These points are best when they come from the face because the face is the most non-spherical part of the head. Especially the nose, the eyebrows and the cheeks.

Combined together this enables accurate electrode registration with the MRI.





(c) Justin Ales, jma23@st-andrews.ac.uk
