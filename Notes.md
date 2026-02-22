#### Abstract in words
The paper claims that if I have a weight-4 $X$ check, there will be a case where two data qubits contribute to the same boundary spacelike DEM edge for this weight 4 check, and there I'll see an addition of angles of the coherent errors. 
These two data qubits are also going to be connected by a weight 2 $Z$ stabilizer check. If during the $\ket{+}_L$  state, I measured $Z_1Z_2$ for these qubits to be +1, then I can prove that the angles add. If $Z_1Z_2=-1$ during the initial state prep syndrome cycle, and I dont correct this before moving on to the next syndrome cycle where I inject noise, the coherent errors cancel out completely.
If I look at the case where I just don't take the $Z_1Z_2$ stabilizer checks into account at all, and just lump $Z_1Z_2=\pm1$ cases together, then essentially I'm working with a 50/50 incoherent mixture of these two cases. In this case, I am taking an average of cases where either the coherent errors add completely or cancel completely, so I just get the Pauli twirl stochastic answer. 
#### Analytics
The claim is that if $q_1$ and $q_2$ contribute to one spacelike DEM edge corresponding to an $X$ check, the angles add.
Notationally, I write $\rho_L=\prod_{s} \pi_{s}\ket{+}^{\oplus n}\bra{+}^{\oplus n}\pi_s$ for $n=d^2=9$ . 
$\pi_s$ is essentially the _projector_ for a given $Z$ stabilizer measurement $s$.  $\prod_s$ iterates over all $Z$ checks *except* for one: the stabilizer check for $Z_1Z_2$ . We can see visually from the $d=3$ diagram that $Z_1Z_2$ exists as a weight 2 stabilizer check (e.g. if $q_1,q_2$  are the qubits contributing to the DEM edge $p_{22}$ in my diagram above).
The stabilizer check for $Z_1Z_2$ is measured in the first round of syndrome extraction when preparing the quiescent state, and thus defines the Pauli frame. We denote the projector for this mreasurement as $\Pi_\pm$ in all subsequent notation. $\Pi_\pm$ projects out the case where we measured the weight 2 $Z_1Z_2=\pm1$ when preparing the quiescent state. Our initial state then is $\rho(t=0)=\Pi_\pm\rho_L\Pi_\pm$. The other stabilizer checks and their projectors are included in the definition of $\rho_L$.
We now look at $\rho =U_2U_1\rho(t=0)U_1^\dagger U_2^\dagger = U_2U_1\Pi_\pm \rho_L\Pi_\pm U_1^\dagger U_2^\dagger$ where $U_i=e^{-iZ_i\theta_i}$ performs a coherent error.
$$\begin{align}
\rho = &\left[\cos\theta_1\cos\theta_2 -\sin\theta_1\sin\theta_2 Z_1Z_2-i(Z_1\sin\theta_1\cos\theta_2+Z_2\sin\theta_2\cos\theta_1)\right]\Pi_\pm\rho_L\Pi_\pm\\ &\left[\cos\theta_1\cos\theta_2 -\sin\theta_1\sin\theta_2 Z_1Z_2+i(Z_1\sin\theta_1\cos\theta_2+Z_2\sin\theta_2\cos\theta_1)\right]
\end{align}
$$

Now we use the fact that $\Pi_\pm$ tells us whether $Z_1Z_2=\pm 1$
$$\begin{align}
\rho = &\left[\cos(\theta_1\pm\theta_2)-i(Z_1\sin\theta_1\cos\theta_2+Z_2\sin\theta_2\cos\theta_1)\right]\Pi_\pm\rho_L\Pi_\pm\\ &\left[\cos(\theta_1\pm\theta_2)+i(Z_1\sin\theta_1\cos\theta_2+Z_2\sin\theta_2\cos\theta_1)\right]
\end{align}
$$
$$\begin{align}
\rho = &\left[\cos(\theta_1\pm\theta_2)-iZ_1(\sin\theta_1\cos\theta_2+Z_1Z_2\sin\theta_2\cos\theta_1)\right]\Pi\rho_L\Pi\\ &\left[\cos(\theta_1\pm\theta_2)+iZ_1(\sin\theta_1\cos\theta_2+Z_1Z_2\sin\theta_2\cos\theta_1)\right]
\end{align}
$$
$$\begin{align}
\rho = &\left[\cos(\theta_1\pm\theta_2)-iZ_1(\sin\theta_1\cos\theta_2\pm\sin\theta_2\cos\theta_1)\right]\Pi_\pm\rho_L\Pi_\pm\\ &\left[\cos(\theta_1\pm\theta_2)+iZ_1(\sin\theta_1\cos\theta_2\pm\sin\theta_2\cos\theta_1)\right]
\end{align}
$$
$$\begin{align}
\rho = &\left[\cos(\theta_1\pm\theta_2)-iZ_1(\sin(\theta_1\pm\theta_2))\right]\Pi_\pm\rho_L\Pi_\pm\\ &\left[\cos(\theta_1\pm\theta_2)+iZ_1\sin(\theta_1\pm\theta_2)\right]
\end{align}
$$
$$
\rho = \cos^2(\theta_1\pm\theta_2)\Pi_\pm\rho_L\Pi_\pm+\sin^2(\theta_1\pm\theta_2)Z_1\Pi_\pm\rho_L\Pi_\pm Z_1
$$
So probability of $X$ detector firing $P(a_0=-)=\sin^2(\theta_1\pm\theta_2)$.
If we have $Z_1Z_2=1$, we get that $p=\sin^2(2\theta)$. If we have $Z_1Z_2=-1$, we should get $P(a_0=-)=0$ , thus washing out the coherent error.

if we had _not_ taken the stabilizer measurement into account, we would have a 50 / 50 mixture of the two cases:
$$\rho_L = \frac{1}{2}\left(\Pi_+\rho_L\Pi_++\Pi_-\rho_L\Pi_-\right)$$
Which would give us $P=0.5(\sin^2(2\theta)+0)=\sin^2 (2\theta)/2$ which is the _Pauli twirl stochastic answer_. 