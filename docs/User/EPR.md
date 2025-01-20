# Energy Participation Ratio

The energy participation ratio (EPR) is a method used to quantize electromagnetic circuits. This method extracts the parameteres to fully describe the Hamiltonian of a quantum circuit and it also extracts the dissipative parameters of the circuit corresponding to lossy interfaces or regions.

The EPR method reduces the probelm of quantizing an electromagnetic circuit to the following question:

What fraction of the energy in mode, m, is stored in the non-linear element, j?

In the following sections we will describe exactly what we mean by this question and also show how we can answer it.


## Introduction

The Hamiltonian of a linear distributed system which describes an electromagnetic environment is given by:

$$
\hat{H} = \sum_{m=1}^{M} \hbar \omega_m \hat{a_m}^{\dag} \hat{a_m}
$$

where $m$ represents the electromagnetic mode, $\omega_m$ are the angular frequencies of the modes and $\hat{a}_m^{\dag}$ and $\hat{a}_m$ are the bosonic creation and annihilation operators for the given mode. If we are using linear elements in our circuit, such as coplanar waverguides (CPWs), transmission lines or CPW resonators, we can always describe the total energy in our system with the Hamiltonian stated above.

The energy participaition ratio (EPR) method works by describing a superconducting circuit as a linear electromagnetic environment with $\textbf{non-linear\ elements}$ embedded within it. Usually these non-linear elements will be Josephson junctions (JJs), however other non-linear elements like nanowires or supeconducting weaklinks exist. For the rest of this descritpion of the EPR method the only non-linear element we will refer to will be Jospehson junctions.

When using the EPR method we will first need to linearize our non-linear components which will be, in general, qubits (X-mon, transmon, fluxomium, etc.) and then perform an electromagnetic (eigenmode) simulation which will allow us to compute the effects of the non-linearities on the system. We describe the theory of the EPR method in detail below and show how calculating the EPR lets us quantize our circuit and serves as a bridge between the classical and quantum domain. By asking the question how much of the energy in mode, m, is stored in the josephson junction, j, we find a way to determine the parameters of our Hamiltonian. This idea will become more clear as we discuss the implementation of the EPR method.

## Josephson Junction

The energy contained in a Josephson junction can be calculated as follows:

$$
 \mathcal{E} = \int dt\ V(t)\ I(t)\
$$

where $V(t)$ is the voltage across the junction and $I(t)$ is the current through the junction. Using the Josephson relations for current and voltage:

$$
    I(t) = I_c \sin{\left( \frac{\Phi(t)}{\phi_0} \right)},\ \ \ \ \ \frac{d\Phi(t)}{dt} = V(t),
$$

where $\Phi(t)$ is the magnetic flux, $I_c$ is the critical current, and $\phi_0 = \hbar/2e$ is the reduced flux quantum. We can calculate the energy of the Josephson junction as follows, note that we now drop the time dependence for $\Phi(t)$.

$$
\begin{aligned}
\mathcal{E}(\Phi) &= \int \frac{d\Phi}{dt}\ I_c \sin{\left( \frac{\Phi}{\phi_0} \right)}\ dt\\
&= \int I_c \sin{\left( \frac{\Phi}{\phi_0} \right)}\ d\Phi\\\
&= - I_c \phi_0 \cos{\left( \frac{\Phi}{\phi_0} \right)}
\end{aligned}
$$

The prefactor is known as the Josephson energy, $E_J = I_c \Phi_0$, and this quantity is proportional to the rate of tunneling of cooper pairs back and forth across the junction. Note that the units of $E_J$ are in Joules. For units in frequency (Hertz) the expression is $E_{J-\text{freq}} = E_J / h$, where $h$ is Planck's constant. Therefore, the common expression for the energy stored in a Josephson junction as seen in the literature is given as:

$$
\begin{aligned}
    \mathcal{E}(\Phi)&= - E_J \cos{\left( \frac{\Phi}{\phi_0} \right)} 
\end{aligned}
$$

We can Taylor expand the cosine term as follows: 

$$
\begin{aligned}
    \cos{\left( \frac{\Phi}{\phi_0} \right)} &= 1 - \frac{1}{2} \left( \frac{\Phi}{\phi_0} \right)^2 + \frac{1}{24} \left( \frac{\Phi}{\phi_0} \right)^4 - \frac{1}{120} \left( \frac{\Phi}{\phi_0} \right)^6 + ...
\end{aligned}
$$

Then, plugging the Taylor expansion into our expression we have:

$$
\begin{aligned}
\mathcal{E}(\Phi)&= - E_J + \frac{E_J}{2} \left( \frac{\Phi}{\phi_0} \right)^2 - \frac{E_J}{24} \left( \frac{\Phi}{\phi_0} \right)^4 + \mathcal{O}  \Big( \Phi^6 \Big) 
\end{aligned}
$$

We drop the constant energy term $-E_J$ as it only constitutes a shift of our potential energy, and we are only interested in differences in potential energy. Furthermore, we take only the first non-linear correction (transmon limit) we have:

$$
\begin{aligned}
\mathcal{E}(\Phi)&\approx \frac{E_J}{2} \left( \frac{\Phi}{\phi_0} \right)^2 - \frac{E_J}{24} \left( \frac{\Phi}{\phi_0} \right)^4 
\end{aligned}
$$

We can separate the energy $\mathcal{E}(\Phi)$ into a linear and non-linear piece where the second-order term (first term) is analagous to the linear inductive energy in the LC oscillator and the rest of the higher-order terms correspond to the non-linear inductance of the JJ:

$$
\begin{aligned}
\mathcal{E}(\Phi) &= \mathcal{E}(\Phi)^{\text{lin}} + \mathcal{E}(\Phi)^{\text{nlin}}\\
        &\approx \underbrace{\frac{E_J}{2} \left( \frac{\Phi}{\phi_0} \right)^2}_{\text{linear}} \underbrace{- \frac{E_J}{24} \left( \frac{\Phi}{\phi_0} \right)^4}_{\text{non-linear}}\\
        &\approx \underbrace{\frac{\Phi^2}{2L_J}}_{\text{linear}} \underbrace{- \frac{E_J}{24} \left( \frac{\Phi}{\phi_0} \right)^4}_{\text{non-linear}}\\
\end{aligned}
$$

where the linear josephson inductance is given by:

$$
\begin{align*}
    L_J = \frac{\phi_0^2}{ E_J} = \frac{\phi_0}{I_c},  
\end{align*}
$$

and will be very important in implementing the EPR method.

## Transmon Hamiltonian

The Hamiltonian of the transmon consists not only of the inductive energy of the JJ but also includes the capacitive energy stored in the capacitor pads and the small capacitance of the junction itself. The Hamiltonian, which we label $\hat{H}_4$ because it is an approximation to fourth order, is expressed as:

$$
\begin{aligned}
    \hat{H}_4 =\underbrace{\frac{\hat{Q}^2}{2C_\Sigma}}_{\text{cap energy}}+ \underbrace{\frac{\hat{\Phi}^2}{2L_J} - \frac{E_J}{24} \left( \frac{\hat{\Phi}}{\phi_0} \right)^4 }_{\text{ind energy}}
\end{aligned}
$$

where $\hat{Q}$ is the charge operator, $C_\Sigma = C_{\text{pad}} + C_{J}$, is the total capacitance including the capacitance of the pads, $C_{\text{pad}}$, and the junction itself, $C_{J}$. Note that the magnetic flux, $\hat{\Phi}$, and charge, $\hat{Q}$, are now operators and they obey the commutation relation $[\hat{\Phi}, \hat{Q}] = i \hbar$. We can now separate the Hamiltonian into a linear and non-linear piece:

$$
\begin{aligned}
    \hat{H}_4 = \underbrace{\frac{\hat{Q}^2 }{2C_\Sigma}+ \frac{\hat{\Phi}^2}{2L_J}}_{\text{linear}} \underbrace{- \frac{E_J}{24} \left( \frac{\hat{\Phi}}{\phi_0} \right)^4 }_{\text{non-linear}}.
\end{aligned}
$$

The linear portion is the Hamiltonian for the quantum LC oscillator and the non-linearity is kept only up to the fourth order term as higher order terms have a diminishgly small contribution.

We can express the flux operator, $\hat{\Phi}_j$, for a junction, j, in terms of the creation and annihlation operators as follows:

$$
\hat{\Phi}_j = \Phi_{\text{zpf}}\left(\hat{a}^\dagger + \hat{a} \right),
$$

where $ \Phi_{\text{zpf}}$ is the zero point fluctuations of a given mode in our superconducting circuit. Then we can express our Hamiltonian as follows:

$$
\begin{aligned}
\hat{H}_4 = \underbrace{ \hbar \omega \hat{a}^{\dag} \hat{a}}_{\text{linear}}\  \ \underbrace{- \frac{E_J\ \Phi_{\text{zpf}}^4}{24\ \phi_0^4} \left( \hat{a}^\dagger + \hat{a} \right)^4 }_{\text{non-linear}}.
\end{aligned}
$$

The only unknown quanity in $\hat{H}_4$ is $\Phi_{\text{zpf}}$, which will be obtained from the energy participation ratio. However, we first must describe the Hamiltonian for the whole electromagnetic environment not just the transmon.

## Full Hamiltonian of the System

The Hamiltonian of the full system includes not just a single transmon but also potentially other transmon qubits and readout resonators. Recall that resonators are modelled as LC oscillators and are considered linear elements. We therefore, write the Hamiltonina of the whole supeconducting circuit to include all linear modes.

Note that because we have broken the transmon Hamiltonian up into a linear and non-linear piece, we can now seperate the two, and lump the linear piece of the transmon in with the other linear elements including resonators. This process is called linearizing the qubit. The non-linear pieces arising from the Jospehson junctions in the circuit are kept separate and we can express the full Hamiltonian as follows:

$$
\begin{aligned}
    \hat{H}_{\text{full}} = \sum_{m=1}^{M} \hbar \omega_m \hat{a}^{\dag}_m \hat{a}_m - \frac{E_J\ }{24\ \phi_0^4}  \sum_{m=1}^{M} \sum_{j=1}^{J}  (\Phi_{j,m}^{\text{zpf}})^4 \left( \hat{a}_m^\dagger + \hat{a}_m \right)^4
\end{aligned}
$$

The above expression is broken up into the linear piece: 
$$
\hat{H}_{\text{lin}} =  \sum_{m=1}^{M} \hbar \omega_m \hat{a}^{\dag}_m \hat{a}_m 
$$

and non-linear piece:

$$
\hat{H}_{\text{full}} = - \frac{E_J\ }{24\ \phi_0^4}  \sum_{m=1}^{M} \sum_{j=1}^{J}  (\Phi_{j,m}^{\text{zpf}})^4 \left( \hat{a}_m^\dagger + \hat{a}_m \right)^4
$$


Note that the Hamiltonian we describe here is just for transmon qubits a more general Hamiltonian we include the higher order terms not just the fourth order term. We will now use the energy participation ratio to calculate $\Phi_{\text{zpf}}$.

## Calculating the Energy Participation Ratio

The energy participation ratio of a junction in mode $m$ is given by:

$$
\begin{aligned}
    \text{p}_m = \frac{\text{Inductive energy stored in the junction in mode m}}{\text{Total inductive energy in the circuit in mode m}}
\end{aligned}
$$

The calculation is performed as follows:

$$
\begin{aligned}
    \text{p}_m = \frac{\langle \psi_m | \frac{1}{2 L_J} \hat{\Phi}^2| \psi_m \rangle}{\langle \psi_m | \frac{1}{2}\hat{H}_{\text{lin}}| \psi_m \rangle}
\end{aligned}
$$

where we consider consider $|\psi_m \rangle $ as a number state (fock state) of the mode m. Computing these expectation values we get:

$$
\begin{aligned}
    \text{p}_m = \frac{ \frac{1}{ L_J} (\Phi_{m}^{\text{zpf}})^2}{ \frac{1}{2} \hbar \omega_m},
\end{aligned}
$$

and we can rearrange to solve for $(\Phi_{m}^{\text{zpf}})^2$ :

$$
\begin{aligned}
    (\Phi_{m}^{\text{zpf}})^2  = \frac{ \text{p}_m \hbar \omega_m L_j }{ 2}.
\end{aligned}
$$

We now have an expression for $\Phi_{\text{zpf}}$ in terms of the participation ratio. In the next section we will describe how to obtain the participation ratio from a finite element electromagnetic simulation using PALACE.
