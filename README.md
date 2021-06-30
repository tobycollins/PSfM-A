# PSfM-A

Matlab implementation of the 2017 PAMI paper "Planar Structure-from-Motion with Affine Camera Models: Closed-Form Solutions, Ambiguities and Degeneracy Analysis" by Toby Collins and Adrien Bartoli. If you use this code please cite the following two papers:


## The main paper

@ARTICLE{7487003,

  author={Collins, Toby and Bartoli, Adrien},

  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence}, 

  title={Planar Structure-from-Motion with Affine Camera Models: Closed-Form Solutions, Ambiguities and Degeneracy Analysis}, 

  year={2017},

  volume={39},

  number={6},

  pages={1237-1255},

  doi={10.1109/TPAMI.2016.2578333}}

@ARTICLE{7487003, author={Collins, Toby and Bartoli, Adrien}, journal={IEEE Transactions on Pattern Analysis and Machine Intelligence}, title={Planar Structure-from-Motion with Affine Camera Models: Closed-Form Solutions, Ambiguities and Degeneracy Analysis}, year={2017}, volume={39}, number={6}, pages={1237-1255}, doi={10.1109/TPAMI.2016.2578333}} 

## The follow-up paper with an improved closed-form solution to camera resection


@article{10.1007/s10851-018-0794-0,
author = {Bartoli, Adrien and Collins, Toby},
title = {Plane-Based Resection for Metric Affine Cameras},
year = {2018},
issue_date = {September 2018},
publisher = {Kluwer Academic Publishers},
address = {USA},
volume = {60},
number = {7},
issn = {0924-9907},
url = {https://doi.org/10.1007/s10851-018-0794-0},
doi = {10.1007/s10851-018-0794-0},
abstract = {We study the problem of resecting the metric affine camera models from at least three non-colinear point correspondences. A direct application is plane pose estimation. We consider the three most popular metric affine cameras, namely the paraperspective, weak-perspective and orthographic cameras. For each model, we give an algebraic procedure which finds the optimal solution, where by optimal we mean the global minimizer of the reprojection error under the Euclidean norm. Our algebraic procedures cover both the minimal case of three points and the redundant cases of more than three points. They always return two solutions, as the problem has a two-way ambiguity on the rotation and translation for the three cameras in the general case. The scale of the paraperspective and weak-perspective cameras is, however, recovered uniquely. The orthographic case is the most involved and has not been solved analytically in the literature. We characterize its intrinsic complexity by showing that it reduces to finding the roots of an irreducible and non-solvable by radicals sextic polynomial. The previous algorithms for the paraperspective and weak-perspective cases have singularities, while, in contrast, our algebraic procedures do not.},
journal = {J. Math. Imaging Vis.},
month = sep,
pages = {1037–1064},
numpages = {28},
keywords = {Resection, Polynomial, Affine camera, Pose, Plane, Optimal}
}

  

