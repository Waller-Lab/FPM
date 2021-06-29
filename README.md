# FPM
Matlab code implements Fourier ptychography reconstruction algorithm from a set of images captured under different illumination angles (e.g. in an LED array microscope), using either sequential or multiplexed coded illumination. The algorithm simultaneously estimates the complex object function (amplitude + phase) withhigh resolution (defined by objective NA+illumination NA) and the pupil function (aberrations). It implements a sequential quasi-Newton’s method with Tikhonov (L2) regularization.

### Dependency:
1. [Nature-order filename sort](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort?s_tid=mwa_osa_a)

### Paper:
Please cite as:

L. Tian, X. Li, K. Ramchandran, and L. Waller, “Multiplexed coded illumination for Fourier ptychography with an LED array microscope,” Biomedical Optics Express 5, 2376-2389 (2014).

[Paper Link](https://www.osapublishing.org/boe/fulltext.cfm?uri=boe-5-7-2376&id=294149)

### Data:
Sample data can be found [here](http://www.laurawaller.com/opensource/).

### Updates:
06/28/2021:
1. Added first version of code to repo.
