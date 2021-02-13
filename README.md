# NOISE-ASSISTED MULTIVARIATE VARIATIONAL MODE DECOMPOSITION (Accepted Paper on ICASSP 2021)

**Authors:** Charilaos A. Zisou, Georgios K. Apostolidis, Leontios J. Hadjileontiadis

**Abstract:** The variational mode decomposition (VMD) is a widely applied optimization-based method, which analyzes nonstationary signals concurrently. Correspondingly, its recently proposed multivariate extension, i.e., MVMD, has shown great potentials in analyzing multichannel signals. However, the requirement of presetting the number of extracted components K diminishes the analytic property of both VMD and MVMD methods. This work combines MVMD with the noise injection paradigm to propose an efficient alternative for both VMD and MVMD, i.e., the noise-assisted MVMD (NA-MVMD), that aims at relaxing the requirement of presetting K, as well as improving the quality of the resulting decomposition. The noise is injected by adding noise variables/channels to the initial signal to excite the filter bank property of VMD/MVMD on white Gaussian noise. Moreover, an alternative approach of updating center frequencies is proposed, which uses the centroid of the generalized cross–spectrum instead of a simple average of the individual spectral centroids, showing faster convergence. The NA–MVMD is applied to both univariate and multivariate synthetic signals, showing improved analytical ability, noise intolerance, and less sensitivity in selecting the K parameter. 


**Acknowledgments:** The VMD and MVMD Matlab codes were taken from the publicly available MATLAB Central File Exchange pages below:

VMD [1]:
Dominique Zosso (2021). Variational Mode Decomposition (https://www.mathworks.com/matlabcentral/fileexchange/44765-variational-mode-decomposition), MATLAB Central File Exchange.

MVMD [2]:
Naveed ur Rehman (2021). Multivariate Variational Mode Decomposition (MVMD) (https://www.mathworks.com/matlabcentral/fileexchange/72814-multivariate-variational-mode-decomposition-mvmd), MATLAB Central File Exchange.

The corresponding papers:

[1] K. Dragomiretskiy and D. Zosso, "Variational Mode Decomposition," in IEEE Transactions on Signal Processing, vol. 62, no. 3, pp. 531-544, Feb.1, 2014, doi: 10.1109/TSP.2013.2288675.

[2] N. u. Rehman and H. Aftab, "Multivariate Variational Mode Decomposition," in IEEE Transactions on Signal Processing, vol. 67, no. 23, pp. 6039-6052, 1 Dec.1, 2019, doi: 10.1109/TSP.2019.2951223.
