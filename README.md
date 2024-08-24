# Maximum likelihood low-rank Toeplitz covariance estimation
This repository hosts an implementation of a primal logarithmic barrier method for solving the problem

$$
\begin{array}{r}
\text{minimize} & \text{log det } (T + \sigma^2 I) + \textbf{Tr}((T + \sigma^2 I)^{-1} S) + \lambda \| T \|_{*} \\
\text{subject to} & T \text{ being Toeplitz}, \h T \succeq 0,  \\
&  \sigma^2_{\text{ub}} \geq \sigma^2 \geq \sigma^2_{\text{lb}}
\end{array}
$$


$$
\begin{array}{r}
\text{minimize} & \text{log det } (T + \sigma^2 I) + \textbf{Tr}((T + \sigma^2 I)^{-1} S) + \lambda \| T \|_{*} \\
\text{subject to} & T \text{ being Toeplitz}, \hspace{0.4cm} T \succeq 0, \hspace{0.4cm} \sigma^2_{\text{ub}} \geq \sigma^2 \geq \sigma^2_{\text{lb}}, \hspace{2cm}
\end{array}
$$

with decision variables
$T \in \mathbf{H}^{n+1}$ and $\sigma^2 \in \reals$,
and problem data
$S \in \mathbf{H}^{n+1}$ representing the sample covariance matrix and noise bounds $\sigma^2_{\text{ub}} > \sigma^2_{\text{lb}}$.
(Here $\mathbf{H}^{n+1}$ is the space of Hermitian matrices of dimension $n + 1$.)
This problem is used for estimating a covariance matrix that is known to be the sum of a scaled identity matrix and a low-rank
positive semidefintie Toeplitz matrix. A
description of the method along with possible
applications can soon be found in our forthcoming paper.
