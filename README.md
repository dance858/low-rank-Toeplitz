# Maximum likelihood low-rank Toeplitz covariance estimation
This repository hosts an implementation of a primal logarithmic barrier method for solving the problem

$$
\begin{array}{r}
\text{minimize} & \text{log det } (T + \sigma^2 I) + \textbf{Tr}((T + \sigma^2 I)^{-1} S) + \lambda |T|_* \\
\text{subject to} & T \text{ being Toeplitz}, \hspace{0.4cm} T \succeq 0, \hspace{0.4cm} \sigma^2_1 \geq \sigma^2 \geq \sigma^2_0 , \hspace{0.8cm} 
\end{array}
$$

with decision variables
$T \in \mathbf{H}^{n+1}$ and $\sigma^2 \in \mathbf{R}$,
and problem data
$S \in \mathbf{H}^{n+1}$ representing the sample covariance matrix and noise bounds $\sigma^2_{\text{1}} > \sigma^2_{\text{0}}$.
(Here $\mathbf{H}^{n+1}$ is the space of Hermitian matrices of dimension $n + 1$ and $|T|_*$ is the nuclear norm of $T$.)
This problem is used for estimating a covariance matrix that is known to be the sum of a scaled identity matrix and a low-rank
positive semidefinite Toeplitz matrix. A description of the method along with possible applications can soon be found in our forthcoming paper.

To recreate Figure 1 of the paper, simply run the script MSE_vs_samples.m.
