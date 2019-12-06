---
title: "Finding NHIM using LD: N DoF quadratic normal form Hamiltonian"
author: "Shibabrat Naik"
geometry: margin = 1in
output:
  pdf_document:
    fig_caption: yes
    fig_height: 1
  html_document:
    fig_caption: yes
    fig_height: 1
bibliography: ./Ham_LD_reaction.bib
csl: elsevier-without-titles.csl
---

__Abstract__

Phase space structures such as dividing surfaces, normally hyperbolic invariant manifolds, their stable and unstable manifolds have been an integral part of computing quantitative results such as transition fraction, stability erosion in multi-stable mechanical systems, and reaction rates in chemical reaction dynamics. Thus, methods that can reveal their geometry in high dimensional phase space (4 or more dimensions) need to be benchmarked by comparing with known results. In this study we assess the capability of one such method called _Lagrangian descriptor_ for revealing the types of high dimensional phase space structures associated with index-1 saddle in Hamiltonian systems. The Lagrangian descriptor based approach is applied to two and three degree-of-freedom quadratic Hamiltonian systems where the high dimensional phase space structures are known, that is as closed-form analytical expressions. This leads to a direct comparison of features in the Lagrangian descriptor plots and the phase space structures' intersection with an isoenergetic two-dimensional surface, and hence provides a validation of the approach.


## Introduction

Identifying invariant manifolds in high dimensional phase space and
their anchor, normally hyperbolic invariant manifold (NHIM), is the
stepping stone in applying phase space transport methods to a broad
array of problems in physics, chemistry, and
engineering [@wiggins90; @Gottwald1995; @Jaffe1999; @DeOliveira2002; @Dellnitz2005; @Gabern2005; @Gabern2006; @waalkens08; @ww10; @wiggins2013normally; @wig2016].
These problems are typically formulated in phase space of more than 4
dimensions (that is, $N = 2$ or more degrees of freedom) and the
geometric approach of computing invariant manifolds that are
codimension-1 separatrices requires Poincaré-Birkoff normal form theory.
The Poincaré-Birkoff normal form approach for computing NHIMs and their
stable and unstable manifolds has only been developed for Hamiltonian
systems in the neighborhood of index-k saddle points. This approach in
itself is successful and rigorous, but its implementation requires
experience in writing or using algebraic manipulation programs.
Furthermore, these computations become expensive as the dimensions
increase since Poincaré-Birkoff normal form theory is based on Taylor
expansion, a high dimensional vector valued polynomial, of the
Hamiltonian in a neighborhood of the saddle point to sufficiently high
order. Besides, the resulting visualization can be non-obvious leading
to cumbersome interpretation and can not work for high dimensions. The
succcess and shortcoming of the vizualization approach based on
topological methods for four-dimensional space becomes apparent when
trying to extend the results in
Refs. [@Kocak1986; @Banks1992; @Hoffmann1991]. This raises the question
of detecting the signatures of the high dimensional structures, if they
exist, when probed with low dimensional surfaces. If one does manage to
detect these structures, how do they manifest on the low dimensional
surfaces.

One such method is the Lagrangian descriptors (LDs) developed in
Refs. [@madrid2009; @mendoza2010; @mancho2013; @lopesino2017].
Lagrangian descriptors are a class of trajectory diagnostic methods that
can reveal phase space structures by encoding geometric property of
trajectories (such as, phase space arc length, configuration space
distance or displacement, cumulative action or kinetic energy)
initialised on a two dimensional surface. The method was originally
developed in the context of Lagrangian transport studies in fluid
dynamics that requires identifying transport barriers which are the
invariant manifolds in time dependent 2D fluid flow. Besides, the method
applies to both Hamiltonian and non-Hamiltonian systems
([@lopesino2017]) as well as to systems with arbitrary, stochastic and
dissipative, time-dependence
([@balibrea2016lagrangian; @craven2017lagrangian; @junginger2016lagrangian]).
Lagrangian descriptor based detection of transport barriers has also
been applied directly to data sets, such as those obtained from
satellite observations or numerical
simulations ([@amism11; @mendoza2014; @ggmwm15; @ramos2018]).
Furthermore, the method is straightforward to implement computationally
and it provides a “high resolution” method for exploring the influence
of high dimensional phase space structure on trajectory behaviour. The
method of Lagrangian descriptors takes an [*opposite*]{} approach to
that of Lyapunov exponent type calculations by emphasizing the initial
conditions of trajectories, rather than their advected locations that is
involved in calculating normalized rate of divergence. This is achieved
by considering a two dimensional section of the full phase space and
discretizing with a dense grid of initial conditions. Even though the
trajectories wander off in the phase space, as the initial conditions
evolve in time, there is no loss in resolution of the two dimensional
section. Our objective is to clarify the use of Lagrangian descriptors
as a diagnostic on two dimensional sections of high dimensional phase
space structures. This diagnostic is also meant to be used as the
preliminary step in computing the NHIM, their stable and unstable
manifolds using other computational
means [@junginger2016transition; @bardakcioglu2018; @ezra_2018]. In this
article, we assess the capability of Lagrangian descriptors to detect
the type of high dimensional phase space structures in Hamiltonian
systems such as the NHIM, their stable, and unstable manifolds that are
used in computing rates of chemical reactions.

Computing chemical reaction rates is based on the fundamental framework
of transition state theory as formulated in phase space by Polanyi,
Evans, Wigner,
Eyring [@evans1935applications; @wigner1938transition; @eyring1938theory; @garrett_2000]).
Further research has shown that constructing a *locally* recrosssing
free, minimal flux orientable surface called a [*dividing surface*]{}
(DS) is the phase space structure that provides the correct estimate for
reaction rates. This dividing surface of geometry $\mathbb{S}^{2N - 2}$
( a 2N-2 dimensional sphere) is constructed from the NHIM of geometry
$\mathbb{S}^{2N - 3}$ which forms the equator of the dividing surface,
on a constant energy
surface [@waalkens2004direct; @waalkens08; @ww10; @wig2016].
Furthermore, the global dynamics of reactive and non-reactive
trajectories is governed by the invariant (stable and unstable)
manifolds, of the NHIM, are $\mathbb{S}^{2N - 3} \times \mathbb{R}$, and
thus act as codimension-1 separatrices on the energy surface. Thus, the
NHIM acts as the anchor for the local dynamics via the dividing surface
from reactants to products or vice-versa and for the global dynamics via
the stable/unstable invariant manifolds. Thus, detecting and
constructing the NHIM forms a fundamental step in obtaining quantitative
results in reaction dynamics [@ezra_2009; @ezra_2018]. In the present
study, we will focus on detecting and verifying these high dimensional
phase space structures.

Recently the applicability of the Lagrangian descriptor based approach
to time dependent problems with random and dissipative forcing in
chemical reaction dynamics has been shown whereby the transition state
trajectory is calculated using the extremal values in Lagrangian
descriptor values; see
Refs. [@craven2015lagrangian; @craven2016deconstructing; @craven2017lagrangian; @junginger2016uncovering; @junginger2016transition; @junginger2017chemical; @junginger2017variational; @feldmaier2017obtaining; @revuelta2017transition; @patra2018detecting].
The initial conditions for the transition state trajectory is identified
by computing the extrema of the Lagrangian descriptor on a two
dimensional domain. Comparing extremal and singular features in the
Lagrangian descriptor plot with invariant manifolds in three dimensional
vector fields, Ref. [@garcia-garrido_2018] has also provided numerical
evidence for detecting NHIM using the LD on two dimensional surfaces.
The authors presented the comparison using a stacked version of a two
dimensional linear saddle, the duffing oscillator, and a time-perturbed
3D geophysical model. This has the advantage that the NHIM and the
invariant manifolds are known exactly or can be computed with other
established methods. For the problems considered therein, the NHIM (a
curve embedded in 3D) and its associated invariant manifolds (a 2D
surface in 3D) can be visualized in the 3D space since the
dimensionality of phase space structures is less than 3. We take this
approach to the next logical step by applying this to benchmark problems
of two and three degrees-of-freedom Hamiltonian system where the phase
space structures are known exactly, that is they can be written as
closed-form analytical expressions This article is outlined as follows.
In Section \[sec:theory\], we describe the method of Lagrangian
descriptor used in this study and present an analytical result on
identifying invariant manifolds using features in the Lagrangian
descriptor. In Section \[sec:normal\_form\_2dof\]
and \[sec:normal\_form\_3dof\], we discuss the benchmark systems and
present numerical evidence of our claim using isoenergetic
two-dimensional surfaces. In Section \[sec:summ\], we summarize our
results on detecting invariant manifolds and discuss related future
directions.

## Method of Lagrangian Descriptor {#sec:LD}

The Lagrangian descriptor (LD) as presented in Ref.[@madrid2009] is the
arc length of a trajectory calculated on a chosen initial time $t_0$ and
measured for fixed forward and backward integration time, $\tau$. For
continuous time dynamical systems, Ref.[@lopesino2017] gives an
alternative definition of the LD which is useful for proving rigorous
results and can be computed along with the trajectory. It provides a
characterization of the notion of singular features of the LD that
facilitates a proof for detecting invariant manifolds in certain model
situations. In addition, the “additive nature” of this new definition of
LD provides an approach for assessing the influence of each
degree-of-freedom separately on the Lagrangian descriptor. This property
was used in Ref.[@demian2017] which showed that a Lagrangian descriptor
can be used to detect Lyapunov periodic orbits in the two
degrees-of-freedom H[é]{}non-Heiles Hamiltonian system. We will describe
this procedure for two and three degrees-of-freedom linear autonomous
Hamiltonian systems. We begin by establishing notation in the general
setting of a time-dependent vector field where

$$\frac{d\mathbf{x}}{dt} = \mathbf{v}(\mathbf{x},t), \quad \mathbf{x} \in \mathbb{R}^n \;,\; t \in \mathbb{R}$$

where $\mathbf{v}(\mathbf{x},t) \in C^r$ ($r \geq 1$) in $\mathbf{x}$
and continuous in time. The definition of LDs depends on the initial
condition $\mathbf{x}_{0} = \mathbf{x}(t_0)$, on the initial time $t_0$
(trivial for autonomous systems) and the integration time $\tau$, and
the type of norm of the trajectory’s components, and takes the form,

\begin{equation}
M_p(\mathbf{x}_{0},t_0,\tau) = \displaystyle{\int^{t_0+\tau}_{t_0-\tau} \sum_{i=1}^{n} |\dot{x}_{i}(t;\mathbf{x}_{0})|^p \; dt}
\label{eqn:M_function}
\end{equation}

where $p \in (0,1]$ and $\tau \in \mathbb{R}^{+}$ are freely chosen
parameters, and the overdot symbol represents the derivative with
respect to time. It is to be noted here that there are three
formulations of the function $M_p$ in the literature: the arc length of
a trajectory in phase space [@madrid2009], the arc length of a
trajectory projected on the configuration
space  [@junginger2016transition; @junginger2016uncovering; @junginger2017chemical; @junginger2017variational],
and the sum of the $p$-norm of the vector field
components [@lopesino_2015; @lopesino2017]. Although the latter
formulation of the Lagrangian descriptor  developed in
Ref. [@lopesino_2015; @lopesino2017] does not resemble the arc length,
the numerical results using either of these forms have been shown to be
in agreement and promise of predictive capability in geophysical
flows [@amism11; @mendoza2014; @ggmwm15; @ramos2018]. The formulation we
adopt here is motivated by the fact that this allows for proving
rigorous result, which we will discuss in the next section, connecting
the singular features and minimum in the LD plots with NHIM and its
stable and unstable manifolds. It follows from the result that

$$\begin{aligned}
\mathcal{W}^s(\mathbf{x}_0, t_0) & = \text{\rm argmin} \; \mathcal{L}^{(f)}(\mathbf{x}_0, t_0, \tau) \\
\mathcal{W}^u(\mathbf{x}_0, t_0) & = \text{\rm argmin} \; \mathcal{L}^{(b)}(\mathbf{x}_0, t_0, \tau)\end{aligned}$$

where the stable and unstable manifolds
($\mathcal{W}^s(\mathbf{x}_0, t_0)$ and
$\mathcal{W}^u(\mathbf{x}_0, t_0)$) denote the invariant manifolds at
intial time $t_0$ and $\text{\rm argmin} \; (\cdot)$ denotes the
argument that minimizes the function
$\mathcal{L}^{(\cdot)}(\mathbf{x}_0, t_0, \tau)$ in forward and backward
time, respectively. In addition, the coordinates of the NHIM at time
$t_0$ is given by the intersection $\mathcal{W}^s(\mathbf{x}_0, t_0)$
and $\mathcal{W}^u(\mathbf{x}_0, t_0)$ of the stable and unstable
manifolds, and thus given by

$$\begin{aligned}
\mathcal{M}(\mathbf{x}_0, t_0) & = \text{\rm argmin} \; \left( \mathcal{L}^{(f)}(\mathbf{x}_0, t_0, \tau) + \mathcal{L}^{(b)}(\mathbf{x}_0, t_0, \tau) \right) = \text{\rm argmin} \; \mathcal{L}(\mathbf{x}_0, t_0, \tau) \qquad \text{NHIM}\end{aligned}$$



## Two Degrees of Freedom {#sec:normal_form_2dof}


In this section we use the Lagrangian descriptor method to identify the
NHIM and its stable and unstable manifolds for a 2 DoF separable
quadratic Hamiltonian with an index-1 saddle at the origin. The
advantage of using this model Hamiltonian is that we can compare the
analytical expressions for the phase space structures with features in
LD plots.

### Decoupled quadratic Hamiltonian

We consider a quadratic Hamiltonian for a two degree-of-freedom system
given by
$$H(q_1, p_1, q_2, p_2) = \underbrace{\frac{\lambda}{2} (p_1^2 - q_1^2)}_\text{$H_r$} + 
\underbrace{\frac{\omega_2}{2}(q_2^2 + p^2_2)}_\text{$H_b$},
\quad 
\lambda,\omega_2 > 0 
\label{eqn:sqh_2dof}$$ with the corresponding vector field
$$\begin{aligned}
\dot{q}_1 = & \phantom{-} \frac{\partial H_2}{\partial p_1} = \phantom{-}\lambda p_1, \\
\dot{p}_1 = & -\frac{\partial H_2}{\partial q_1} = \phantom{-}\lambda q_1,\\
\dot{q}_2 = & \phantom{-}\frac{\partial H_2}{\partial p_2} = \phantom{-}\omega_2 p_2,\\
\dot{p}_2 = & -\frac{\partial H_2}{\partial q_2} = -\omega_2 q_2\\
\end{aligned}
\label{eqn:eom_nf_2dof}$$ The equilibirum point is located at
$(0,0,0,0)$ on the zero total energy surface. It is trivial to check
that the eigenvalues of the linearized system around the equilibrium
point are $\pm \lambda$ and $\pm i \omega_2$, and hence the equilibrium
point is of saddle $\times$ center type. In this form, the Hamiltonian 
is decoupled into the “reactive” mode given by $H_r$ and the “bath” mode
given by $H_b$, hence it will be referred to as separable quadratic
Hamiltonian (SQH). This representation lends to discussing the
distribution of total energy between the two modes for the phase space
structures in uncoupled coordinates. In this form, a chemical reaction
is said to have occurred when the $q_1$ coordinate of a trajectory
changes sign and thus, an isoenergetic, $H(q_1, p_1, q_2, p_2) = h$,
dividing surface (DS) can be defined by $q_1 = 0$ hypersurface. The
constant energy defines a three dimensional surface in the four
dimensional phase space and is given by
$$\frac{\lambda}{2} (p_1^2 - q_1^2) + \frac{\omega_2}{2}(p_2^2 + q^2_2) = H_r + H_b = h > 0, \qquad H_r > 0, \quad H_b \geq 0$$
The dividing surface, $q_1 = 0$, for a constant energy is
$$\frac{\lambda}{2} p_1^2 + \frac{\omega_2}{2}(p_2^2 + q^2_2) = H_r + H_b = h > 0, \qquad H_r > 0, \quad H_b \geq 0 \label{eqn:ds_sqh_2dof}$$
which is a two dimensional surface, or precisely of geometry
$\mathbb{S}^2$, that is a 2-sphere on the three dimensional energy
surface. Thus it is codimension-1 and partitions the energy surface into
reactant $p_1 - q_1 > 0$ and product $p_1 - q_1 < 0$ regions by the
forward and backward “reaction” dividing surfaces as shown in
Ref. [@waalkens2004direct] and given by $$\begin{aligned}
p_1 & = & \pm \sqrt{\frac{2}{\lambda}\left( H_r + H_b - \frac{\omega_2}{2}(p_2^2 + q_2^2) \right)}, \quad \text{forward/backward DS} 
\end{aligned}$$ The forward and backward DS is joined at $p_1 = 0$, thus
$$\mathcal{M}(h) = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; q_1 = 0, p_1 = 0, \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) = H_b = h \geq 0 \right\} \qquad \text{NHIM}
\label{eqn:sep_quad_ham2dof_nhim}$$ which is of geometry $\mathbb{S}^1$,
that is a circle centered at the origin and of radius
$\sqrt{2h/\omega_2}$ in the $(q_2,p_2)$ plane. This is a [*normally
hyperbolic invariant manifold*]{} (NHIM) associated with the index-1
saddle and parametrized by total energy
$H(q_1, p_1, q_2, p_2) = h$ [@wiggins2013normally]. [*Invariance*]{}
follows from the vector field  since when $q_1 = p_1 =0$,
$\dot{q}_1 = \dot{p}_1 =0$. Thus $q_1$ and $p_1$ always remain zero, and
trajectories with these initial conditions remain on the NHIM, that is,
$q_1 = p_1 =0$ is invariant. It is [*normally hyperbolic*]{} since the
directions normal to the NHIM, that is, the $(q_1, p_1)$ surface, have
linear saddle like dynamics. For a two degree-of-freedom system, this
NHIM is more commonly referred in the literature as an unstable periodic
orbit.

In order to understand the relationship between the NHIM and the index 1
saddle point, we note that for $H_r =0$ and $H_b =0$ the NHIM reduces to
the point $(q_1, p_1, q_2, p_2) = (0, 0, 0, 0)$, which is the index-1
saddle point on the energy surface $H_r + H_b = 0$. Therefore, as the
total energy is increased from $0$, with $H_b$ increasing from zero, the
NHIM “grows” from the index-1 saddle point on the zero energy surface
into an invariant 1-sphere. This shows how the “influence” of the
index-1 saddle point is carried to higher energy sufaces on which the
saddle point does not exist.

The stable and unstable manifolds of the NHIM  are given by
\begin{align}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) =& \big\{  (q_1, p_1, q_2, p_2) \; | \; q_1 = p_1, \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) = H_b > 0  \big\}, \label{eqn:quad_ham2dof_umani}\\
\mathcal{W}^{\rm s}(\mathcal{M}(h)) =& \big\{  (q_1, p_1, q_2, p_2) \; | \; q_1 = -p_1, \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) = H_b > 0  \big\}, \label{eqn:quad_ham2dof_smani}
\end{align}

which is two dimensional surface and have geometry
$\mathbb{R} \times \mathbb{S}^1$ for a fixed energy. Thus, the
codimension-1 geometry of the manifolds partition the phase space into
“reactive” and “non-reactive” trajectories as shown in
Fig. \[fig:ps\_struct\_quad\_ham2dof\]. The detailed visualization is
only possible for the two degrees-of-freedom system which is also the
starting point for testing the Lagrangian descriptor based approach for
detecting the phase space structures.

![](./figures/ps_struct_quad_ham2dof.png){width="98.00000%"}

Fig. 1: Shows the phase space structures, namely the NHIM which is an
unstable periodic orbit (as black circle), its stable and unstable
manifolds (blue and red surfaces), and the energy surface (as yellow
surface) for the two degrees-of-freedom decoupled quadratic
Hamiltonian . Four example trajectories illustrate the dynamics mediated
by these phase space structures: the red and green trajectories have the
same initial configuration space coordinates, but are located outside
and inside the codimension-1 invariant manifolds which separate the
“reactive” (green) and “non-reactive” (red) trajectories. Also shown are
the two trajectories on the manifolds (magenta on the stable manifold
and cyan on the unstable manifold) that show the dynamics on the
invariant surfaces. (b) Shows the projection of the stable (blue) and
unstable (red) manifolds in the configuration space coordinates
$(q_1, q_2)$ along with the projection of the four example trajectories.
The equipotential contours are the black lines and the region of
inaccessible motion is shadded as grey.

#### Detecting the unstable periodic orbit and its manifolds

To identify NHIM and its invariant manifolds, we compute the Lagrangian
descriptor in a square domain of size 2 units around the origin and
discretize the coordinates of the two-dimensional surface. Next, we pick
a constant value for one of the two remaining coordinates, and use the
total energy equation to solve for the fourth coordinate. Due to the
form of the Hamiltonian , obtaining the coordinate from the constant
energy condition reduces simply to solving a quadratic equation.

* Isoenergetic two-dimensional surface
parametrized by $(q_1, p_1)$  On the constant energy surface,
$H(q_1, p_1, q_2, p_2) = h$, we compute Lagrangian descriptor on a
two-dimensional surface parametrized by $(q_1, p_1)$ coordinates by
defining
$$U_{q_1p_1}^+ = \left\{ (q_1, p_1, q_2, p_2) \; | \; q_2 = 0, \dot{q}_2 > 0 : p_2(q_1, p_1, q_2; h) > 0 \right\}$$
where 
$$\begin{aligned}
p_2(q_1, p_1, q_2 = 0; h) = \sqrt{\frac{2}{\omega_2}\left( h - \frac{\lambda}{2}(p_1^2 - q_1^2) \right)} 
\end{aligned}$$ 
The intersection of the two-dimensional surface $U_{q_1p_1}^+$ with the NHIM  becomes

$$\begin{aligned}
\mathcal{M}(h) \cap U_{q_1p_1}^+ = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; \right. & \left. p_1 = 0, q_1 = 0, q_2 = 0, \dot{q}_2 > 0 : \right. p_2(q_1, p_1, q_2; h) > 0 \right\}.\end{aligned}$$

Thus, the NHIM is located at the origin $(0,0)$ and marked by a red
cross in the LD plot (Fig. \[fig:LD\_SQH\_2dof\_q1p1\]). Furthermore,
the intersection of the two-dimensional surface with the unstable  and
stable manifolds  is given by

$$\begin{aligned}
\mathcal{W}^u(\mathcal{M}(h)) \cap U_{q_1p_1}^+ = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; \right. & \left. p_1 = q_1, q_2 = 0, \dot{q}_2 > 0 : \right. \nonumber \\  
& \left. p_2(q_1, p_1, q_2;h) > 0 \right\}, \\
\mathcal{W}^s(\mathcal{M}(h)) \cap U_{q_1p_1}^+ = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; \right. & \left. p_1 = -q_1, q_2 = 0, \dot{q}_2 > 0 : \right. \nonumber \\  
& \left.  p_2(q_1, p_1, q_2;h) > 0 \right\},\end{aligned}$$

which are one-dimensional for a fixed energy, and represent lines
passing through the origin shown as dashed red (unstable) and white
(stable) lines, respectively, in Fig. \[fig:LD\_SQH\_2dof\_q1p1\]. The
only points of local minima in the LD plot
(Fig. \[fig:LD\_SQH\_2dof\_q1p1\]) also lie along the lines passing
through the origin and correspond to the manifolds of the NHIM.

![\label{fig:LD_SQH_2dof_q1p1}](./figures/linear_symp_trans_2dof_M400x400_E2e-01_q1p1.png){width="33.00000%"}\ ![](./figures/linear_symp_trans_2dof_M400x400_E2e-01_q2p2.png){width="33.00000%"}\ ![](./figures/linear_symp_trans_2dof_M400x400_E2e-01_q1q2.png){width="33.00000%"}
Fig. 2: Lagrangian descriptor slice for the two degrees-of-freedom separable quadratic Hamiltonian on the isoenergetic two-dimensional surfaces (a) $U_{q_1p_1}^+$, (b) $U_{q_2p_2}^+$, (c) $U_{q_1q_2}^+$. Parameters used are $\lambda = 1.0, \omega_2 = 1.0$ for constant total energy $H_2 = 0.2$ and integration time $\tau = 10$ is fixed.



* Isoenergetic two-dimensional surface
parametrized by $(q_2, p_2)$  Next, on the fixed energy surface,
$H(q_1, p_1, q_2, p_2) = h$, we compute the Lagrangian descriptor on a
two-dimensional surface parametrized by $(q_2, p_2)$ coordinates by
defining

$$U_{q_2p_2}^+ = \left\{ (q_1, p_1, q_2, p_2) \; | \; q_1 = 0, \dot{q}_1 \geqslant 0 : p_1(q_1, q_2, p_2; h) \geqslant 0 \right\}$$

where $$\begin{aligned}
p_1(q_1 = 0, q_2, p_2; h) =& \sqrt{\frac{2}{\lambda}\left( h  - \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) \right)}\end{aligned}$$
The intersection of the two-dimensional surface $U_{q_2p_2}^+$ with the
NHIM  is given by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{q_2p_2}^+ = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; \right. & \left.   p_1 = 0, q_1 = 0, \dot{q}_1 > 0, \right. \nonumber \\ 
& \left.    \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) = h \right\}.\end{aligned}$$

which represents a circle of radius $\sqrt{2h/\omega_2}$. The radius is
$\approx 0.632$ for $h = 0.2, \omega_2 = 1.0$ and is marked as a dashdot
red circle in Fig. \[fig:LD\_SQH\_2dof\_q2p2\]. The points on the circle
are also locations of minima as shown by the one-dimensional slices at
constant $p_2$.

Next, intersection of the two-dimensional $U_{q_2p_2}^+$ with the
stable  and unstable manifolds  is given by

$$\begin{aligned}
\mathcal{W}^u(\mathcal{M}(h)) \cap U_{q_2p_2}^+ = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; \right. & \left.   p_1 = q_1, q_1 = 0, \dot{q}_1 > 0 : \right. \nonumber \\  
& \left.  \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) = h \right\}, \\
\mathcal{W}^s(\mathcal{M}(h)) \cap U_{q_2p_2}^+ = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; \right. & \left.   p_1 = -q_1, q_1 = 0, \dot{q}_1 > 0 : \right. \nonumber \\  
& \left.  \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) = h \right\},\end{aligned}$$

which are circles of same radius $\sqrt{2h/\omega_2}$. These manifolds
are identified along the boundary of the LD plot in
Fig. \[fig:LD\_SQH\_2dof\_q2p2\] along which the minima also occurs.

* Isoenergetic two-dimensinal surface
parametrized by $(q_1, q_2)$  Next, on the constant energy surface,
$H(q_1, p_1, q_2, p_2) = h$, we compute the Lagrangian descriptor on a
two-dimensional surface by defining

$$U_{q_1q_2}^+ = \left\{ (q_1, p_1, q_2, p_2) \; | \; p_1 = 0, p_2(q_1, p_1, q_2; h) > 0 \right\}$$

where $$\begin{aligned}
p_2(q_1, p_1 = 0, q_2; h) =& \sqrt{\frac{2}{\omega_2}\left( h  + \frac{\lambda}{2}q_1^2 - \frac{\omega_2}{2}q_2^2 \right)}\end{aligned}$$
The intersection of the NHIM  with the two-dimensional surface
$U_{q_1q_2}^+$ is

$$\begin{aligned}
\mathcal{M}(h) \cap U_{q_1q_2}^+ = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; \right. & \left.   p_1 = 0, q_1 = 0, p_2 > 0 : \right. \nonumber \\ 
& \left.    \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) = h \right\}.\end{aligned}$$

which represent points on the line $q_1 = 0$ which is marked by the
dashdot line in Fig. \[fig:LD\_SQH\_2dof\_q1q2\].

Next, intersection of the two-dimensional $U_{q_1q_2}^+$ with the
stable  and unstable manifolds  is given by

\begin{align}
\mathcal{W}^u(\mathcal{M}(h)) \cap U_{q_1q_2}^+ = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; \right. & \left.   p_1 = q_1, p_1 = 0, p_2 > 0 : \right. \nonumber \\  
& \left.  \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) = h \right\}, \\
\mathcal{W}^s(\mathcal{M}(h)) \cap U_{q_1q_2}^+ = \left\{ (q_1, p_1, q_2, p_2) \; \vert \; \right. & \left.   p_1 = -q_1, p_1 = 0, p_2 > 0 : \right. \nonumber \\  
& \left.  \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) = h \right\},
\end{align}

which denote lines parallel to the $q_2$ axis and marked by the dashed
red (unstable) and white (stable) lines in
Fig. \[fig:LD\_SQH\_2dof\_q1q2\]. These manifolds are identified by
minima in the Lagrangian descriptor plot in
Fig. \[fig:LD\_SQH\_2dof\_q1q2\] as also shown in the one-dimensional
slice for constant $q_2$.

Now the question is how do these phase space structures manifest when
the system is non-separable into “reactant” and “product” coordinates.
We will answer this using a linear symplectic transformation of the
separable quadratic Hamiltonian and then compute the LD on different 2D
surfaces for a fixed energy.


### Coupled quadratic Hamiltonian


To couple the coordinates in the separable quadratic Hamiltonian , we
introduce a linear transformation, $C$, such that it satisfies the
symplectic condition $$C \mathcal{J} C^T = \mathcal{J} = \begin{pmatrix}
0_N & I_N \\
-I_N & 0_N
\end{pmatrix}
\label{eqn:symp_cond}$$ where $\mathcal{J}$ is the 2N $\times$ 2N matrix
and $I_N$ denotes the $N \times N$ identity matrix. The symplectic
transformation $C$ acts on the non-separable (coupled) coordinates
$(x, y, p_x, p_y)$ to give the decoupled coordinates
$(q_1, q_2, p_1, p_2)$: $$\begin{bmatrix}
q_1 \\
q_2 \\
p_1 \\
p_2
\end{bmatrix}
= C
\begin{bmatrix}
x \\
y \\
p_x \\
p_y 
\end{bmatrix} 
\implies
\begin{bmatrix}
\dot{q_1} \\
\dot{q_2} \\
\dot{p_1} \\
\dot{p_2}
\end{bmatrix}
= C
\begin{bmatrix}
\dot{x} \\
\dot{y} \\
\dot{p_x} \\
\dot{p_y} 
\end{bmatrix},
\quad
\begin{bmatrix}
\dot{x} \\
\dot{y} \\
\dot{p_x} \\
\dot{p_y} 
\end{bmatrix} 
= -C 
\begin{bmatrix}
\dot{q_1} \\
\dot{q_2} \\
\dot{p_1} \\
\dot{p_2}
\end{bmatrix}$$ since, $C$ satisfies $C^T = C^{-1} = -C$ and some
examples are given in Appendix \[sect:examples\_C\].

Let us consider the symplectic transformation $$C =
\begin{pmatrix}
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \\
-1 & 0 & 1 & 1 \\
0 & -1 & 1 & 1
\end{pmatrix}
\label{eqn:two_dof_C}$$ The coupled (non-separable) coordinates become
$$\begin{aligned}
q_1 &= p_x \\
q_2 &= p_y \\
p_1 &= -x + p_x + p_y \\
p_2 &= -y + p_x + p_y
\end{aligned}$$ This set of transformed coordinates applied to  gives
the non-separable quadratic Hamiltonian $$\begin{aligned}
\mathcal{H}(x, p_x, y, p_y) = \frac{\lambda}{2} x^2 + \frac{\omega_2}{2} y^2 + & \frac{\omega_2}{2} \left( 2p_y^2 + p_x^2 + 2 p_x p_y - 2 y p_x - 2 y p_y \right) + \nonumber \\ 
& \frac{\lambda}{2} \left( p_y^2 + 2p_xp_y - 2xp_x - 2xp_y \right) 
\label{eqn:ham_symp_2dof}\end{aligned}$$ which gives the Hamiltonian
vector field $$\begin{aligned}
\dot{x} & = -\lambda x - \omega_2 y + \omega_2 p_x + (\lambda + \omega_2) p_y \\
\dot{p}_x & = \phantom{-}\lambda( -x + p_x + p_y )  \\
\dot{y} & = -\lambda x - \omega_2 y + (\lambda + \omega_2) p_x + (\lambda + 2\omega_2) p_y \\
\dot{p}_y & =  \phantom{-}\omega_2( -y + p_x + p_y )
\end{aligned}
\label{eqn:eom_symp_2dof}$$ where the equilibrium point is at
$(0,0,0,0)$, and its total energy is $0$. The Jacobian at this
equilibrium point has eigenvalues of type
$\lambda, -\lambda, i\omega_2, -i \omega_2$ (as shown in the
Appendix \[sect:linear\_symp\_transf\_jac\]) and is of saddle $\times$
center type, that is index-1.

In the decoupled (separable) quadratic Hamiltonian, the dividing surface
is given by $q_1 = 0$ Eqn.  which in the transformed coordinates becomes
$p_x = 0$. The dividing surface on the fixed energy surface
$\mathcal{H} = h$ is given by
\begin{equation}
{\rm DS} = \left\{ (x, p_x, y, p_y) \, | \, \frac{\lambda}{2} p_y^2 + \omega_2 p_y^2 - \lambda x p_y - \omega_2 y p_y + \frac{\lambda}{2} x^2 + \frac{\omega_2}{2} y^2 = h \right\}.
\end{equation}
In the decoupled coordinates, the NHIM is defined by $q_1 = 0, p_1 = 0$
on a fixed energy surface, which gives $-x + p_x + p_y = 0$ that is
$p_y = x$ in the coupled coordinates, and is given by
\begin{equation}
\mathcal{M}(h) = \left\{ (x, p_x, y, p_y) \, | \, p_x = 0, p_y = x, \omega_2 x^2 - \omega_2 x y + \frac{\omega_2}{2} y^2 = h \right\}.
\label{eqn:nhim_symp_2dof}
\end{equation}
Next, in the coupled coordinates, the
stable and unstable manifolds of the NHIM are given by 
\begin{align}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) &= \left\{ (x, p_x, y, p_y) \; | \; x = p_y, \frac{\omega_2}{2}\left( (-y + p_x + p_y)^2 + p_y^2 \right)  = h > 0 \right\} \label{eqn:umani_nsqh_2dof} \\
\mathcal{W}^{\rm s}(\mathcal{M}(h)) &= \left\{ (x, p_x, y, p_y) \; | \; x =  2p_x + p_y, \frac{\omega_2}{2}\left( ( -y + p_x + p_y)^2 + p_y^2 \right)  = h > 0 \right\} \label{eqn:smani_nsqh_2dof}
\end{align}

Invertible linear symplectic transformations are clearly $C^\infty$
diffeomorphisms. They not only preserve the Hamiltonian structure, but
they preserve the phase space structures relevant to reaction dynamics.
In particular, Lyapunov exponents are preserved under symplectic
diffeomorphisms which implies that normal hyperbolicity is preserved.
The no-recrossing property of the dividing surface is a result of
transversality of the Hamiltonian vector field to the dividing surface,
and such transversality properties are preserved under diffeomorphisms.
Hence the NHIM  also has a geometry $\mathbb{S}^1$ and the invariant
manifolds  and are also of geometry $\mathbb{R} \times \mathbb{S}^1$.

#### Detecting the unstable periodic orbit and its manifolds

Now we illustrate the procedure for detecting the unstable periodic
orbit, its stable and unstable manifolds by computing the Lagrangian
descriptor on isoenergetic two-dimensional surfaces parametrized by a
pair of coordinates of the separable quadratic Hamiltonian .

* Isoenergetic two-dimensional surface
parametrized by $(x,p_x)$  On a fixed energy surface,
$\mathcal{H}(x, p_x, y, p_y) = h$, we consider a two-dimensional surface
parametrized by $(x, p_x)$ coordinates by defining

$$U_{xp_x}^+ = \left\{ (x,p_x,y,p_y) \; \vert \; y = 0, \; p_y = p_y(x,p_x,y;h) : \dot{y}(x,p_x,y,p_y) > 0 \right\} , 
\label{eqn:sos_symp_2dof_xpx}$$



Combining the dividing surface condition $p_x = 0$, NHIM condition
$p_y = x$, the directionality condition $\dot{y}(x,y,p_x,p_y;h) > 0$,
gives $x > 0$. Thus, for a fixed energy, the intersection of the
one-dimensional NHIM  with the isoenergetic two-dimensional surface  is
given by

$$\begin{aligned}
\mathcal{M}(h) \cap U^{+}_{xp_x} = \left\{ (x,p_x,y,p_y) \; \vert \; \right. & \left.   y = 0, p_x = 0, p_y = x, \right. \nonumber \\ 
& \left.  \omega_2 x^2 - \omega_2 x y + \frac{\omega_2}{2} y^2 = h : x > 0 \right\} \label{eqn:nhim_x_coord} \end{aligned}$$

which is a point on the line $p_x = 0$ and $x = \sqrt{2h/\omega_2}$, and
is marked by the red cross in Fig. \[fig:LD\_NSQH\_2dof\_xpx\]. This
point is also identified by the minimum value of the Lagrangian
descriptor at this coordinate (as shown along the one-dimensional slice
at constant $p_x = 0$) and agrees with  to within the grid resolution
used for the discretization of the two-dimensional surface. Since the
NHIM is the intersection of stable and unstable manifolds, its
coordinate is also a minimum and singular point in the LD values.

Next, the intersection of the unstable and stable manifolds  with the
isoenergetic two-dimensional surface $U_{xp_x}^+$ is given by

$$\begin{aligned}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) \cap U_{xp_x}^+ =  \left\{ (x,p_x,y,p_y) \; \vert \; \right. & \left.   y = 0, x = p_y, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2}\left( ( x + p_x )^2 + x^2 \right)  = h \; : \; (\lambda + \omega_2) p_x + 2 \omega_2 x  > 0 \right\} \\
\mathcal{W}^{\rm s}(\mathcal{M}(h)) \cap U_{xp_x}^+ =  \left\{ (x,p_x,y,p_y) \; \vert \;  \right. & \left.   y = 0, x = 2p_x + p_y, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2}\left( ( x - p_x )^2 + ( x - 2p_x )^2 \right)  = h \; : \; (\lambda + 3 \omega_2) p_x - 2 \omega_2 x < 0  \right\}\end{aligned}$$

where the inequalities are derived from the directionality condition
$\dot{y} > 0$ for the surface. This is detected by the minima in the LD
plot in Fig. \[fig:LD\_NSQH\_2dof\_xpx\] and highlighted by the dashed
red (unstable) and dashed white (stable) curves. Points on these
manifolds are also picked up by the one-dimensional slices along the
constant $p_x = -0.60, 0.30$.

![\label{fig:LD_NSQH_2dof_xpx}](./figures/linear_symp_trans_2dof_M1000x1000_E2e-01_xpx.png){width="33.00000%"}\ ![\label{fig:LD_NSQH_2dof_xy}](./figures/linear_symp_trans_2dof_M1000x1000_E2e-01_xy.png){width="33.00000%"}\ ![\label{fig:LD_NSQH_2dof_ypy}](./figures/linear_symp_trans_2dof_M1000x1000_E2e-01_ypy.png){width="33.00000%"}
Fig. 3: Lagrangian descriptor plot of the non-separable quadratic Hamiltonian vector field~\eqref{eqn:eom_symp_2dof} on the isoenergetic two-dimensional surface (a) $U_{xp_x}^+$, (b) $U_{xy}^+$, (c) $U_{yp_y}^{+}$. The intersection of the NHIM and the isoenergetic two-dimensional surfaces is shown as a red cross (or dash-dot red line) and the corresponding for the manifolds is shown as dashed red (unstable) and dashed white (stable) curves.  The parameters used are $\lambda = \omega_2 = 1.0$, $h = 0.2$, and $\tau = 10$.


* Isoenergetic two-dimensional surface
parametrized by $(x,y)$  On a fixed energy surface
$\mathcal{H}(x, p_x, y, p_y) = h$, we compute the Lagrangian descriptor
on a two-dimensional surface parametrized by $(x, y)$ coordinates by
defining

$$U_{xy}^+ = \left\{ (x,p_x,y,p_y) \; | \; p_x = 0, \; p_y(x,p_x,y;h) > 0 : \dot{p_x}(x, y, p_y) \geqslant 0 \right\} , 
\label{eqn:sos_symp_2dof_xy}$$

Thus, the intersection of the NHIM  with this isoenergetic
two-dimensional surface  is given by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{xy}^+ =  \left\{ (x,p_x,y,p_y) \, | \, \right. & \left.   p_y = x, p_x = 0, p_y > 0: p_y \geqslant x,  \right. \nonumber \\  
& \left.  \omega_2 x^2 - \omega_2 x y + \frac{\omega_2}{2} y^2 = h \right\}\end{aligned}$$

which represents a one-dimensional curve and is marked as a dashdot red
line in Fig. \[fig:LD\_NSQH\_2dof\_xy\]. The NHIM is also identified by
the minima in Lagrangian descriptor values along one-dimensional
sections at constant $y$ in Fig. \[fig:LD\_NSQH\_2dof\_xy\].

Next, the intersection of the unstable and stable manifolds  with the
isoenergetic two-dimensional surface $U_{xy}^+$ is given by

$$\begin{aligned}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) \cap U_{xy}^+ =  \left\{ (x,p_x,y,p_y) \; \vert \; \right. & \left.   x = p_y, p_x = 0, p_y > 0: p_y \geqslant x, \right. \nonumber \\ 
& \left.  \omega_2 x^2 - \omega_2 x y + \frac{\omega_2}{2} y^2 = h \right\} \\
\mathcal{W}^{\rm s}(\mathcal{M}(h)) \cap U_{xy}^+ =  \left\{ (x,p_x,y,p_y) \; \vert \;  \right. & \left.   x = 2p_x + p_y, p_x = 0, p_y > 0: p_y \geqslant x, \right. \nonumber \\ 
& \left.  \omega_2 x^2 - \omega_2 x y + \frac{\omega_2}{2} y^2 = h \right\}\end{aligned}$$

where each manifold represents the same one-dimensional curve as the
NHIM and are marked by dashed red (unstable) and dashed white (stable)
lines in the Fig. \[fig:LD\_NSQH\_2dof\_xy\]. These manifolds are also
identified by points of minima in the Lagrangian descriptor values as
shown by one-dimensional sections at constant $y$ values.

* Isoenergetic two-dimensional surface
parametrized by $(y, p_y)$  On the fixed energy surface
$\mathcal{H}(x, p_x, y, p_y) = h$, we take the Lagrangian descriptor
slice by defining a two-dimensional surface

$$U_{yp_y}^+ = \left\{ (x,p_x,y,p_y) \; | \; x = 0, \, p_x(x,y,p_y;h) > 0 : \dot{x}(x, p_x, y, p_y) > 0  \right\}
\label{eqn:2dsurf_nsqh_2dof_ypy}$$

Thus, the intersection of the NHIM  with this isoenergetic
two-dimensional surface  is given by $$\begin{aligned}
\mathcal{M}(h) \cap U_{yp_y}^+ =  \left\{ (x,p_x,y,p_y) \; \vert \; \right. & \left.   x = 0, p_y = x, p_x = 0, \dot{x}(x, p_x, y, p_y) > 0 : y < 0  \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2} y^2 = h \right\}  \end{aligned}$$
which is a point $y = \sqrt{2h/\omega_2}$ on the $p_y = 0$ line and
shown as a red cross in Fig.\[fig:LD\_NSQH\_2dof\_ypy\]. This point is
also identified by the minima in the Lagrangian descriptor values as
shown along one-dimensional slice in Fig. \[fig:LD\_NSQH\_2dof\_ypy\].

Next, the intersection of the unstable  and stable manifolds  with the
isoenergetic surface $U_{yp_y}^+$ manifest as

$$\begin{aligned}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) \cap U_{yp_y}^+ =  \left\{ (x,p_x,y,p_y) \; \vert \; \right. & \left.   x = 0, x = p_y, \dot{x}(x, p_x, y, p_y) > 0 : y - p_x < 0, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2}\left( - y + p_x \right)^2  = h \right\} \\
\mathcal{W}^{\rm s}(\mathcal{M}(h)) \cap U_{yp_y}^+ =  \left\{ (x,p_x,y,p_y) \; \vert \;  \right. & \left.   x = 0, x = 2p_x + p_y, \right. \nonumber \\ 
& \left.  \dot{x}(x, p_x, y, p_y) > 0 : -\omega_2y + ( \frac{\omega_2}{2} + \lambda )p_y > 0, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2}\left( ( -y + \frac{p_y}{2} )^2 + p_y^2 \right)  = h \right\}\end{aligned}$$

which are shown as dashed red (unstable) and dashed white (stable)
curves in the Fig. \[fig:LD\_NSQH\_2dof\_ypy\]. These curves are also
identified by the points of minima in the Lagrangian descriptor values
as evident by the one-dimensional slices in
Fig. \[fig:LD\_NSQH\_2dof\_ypy\].


## Three Degrees of Freedom {#sec:normal_form_3dof}


In this section we use the Lagrangian Descriptor method to identify the
NHIM and its stable and unstable manifolds for a 3 DoF separable
quadratic Hamiltonian with an index-1 saddle point at the origin. The
advantage of using this model Hamiltonian is that we can compare the
analytical expressions for the phase space structures present with
features in LD plots.


### Decoupled quadratic Hamiltonian


We consider the three degrees of freedom benchmark example in
Ref.[@wig2016] which is a linear system (quadratic Hamiltonian) with an
equilibrium point of saddle-center-center equilibrium type at the
origin. The Hamiltonian of this system is of the form

\begin{equation}
H = \underbrace{\frac{\lambda}{2} \left(p_1^2 - q_1^2 \right)}_{H_r} + \underbrace{\frac{\omega_2}{2} \left(p_2^2 + q_2^2 \right)}_{H_{b_1}} + \underbrace{\frac{\omega_3}{2} \left(p_3^2 + q_3^2 \right)}_{H_{b_2}}, \quad \lambda, \, \omega_2, \, \omega_3 >0
\label{ham3}
\end{equation}

with the corresponding Hamiltonian vector field given by:


\begin{align}
\dot{q}_1 & = \frac{\partial H}{\partial p_1}= \lambda p_1, \\
\dot{p}_1 & = -\frac{\partial H}{\partial q_1}= \lambda q_1, \\
\dot{q}_2 & = \frac{\partial H}{\partial p_2}= \omega_2 p_2, \\
\dot{p}_2 & = -\frac{\partial H}{\partial q_2}= -\omega_2 q_2, \\
\dot{q}_3 & = \frac{\partial H}{\partial p_3}= \omega_3 p_3, \\
\dot{p}_3 & = -\frac{\partial H}{\partial q_3}= -\omega_3 q_3, 
\label{eqn:hameq3}
\end{align}


Since the total Hamiltonian decouples into the independent Hamiltonians
$H_r$, $H_{b_1}$ and $H_{b_2}$ we can analyze the phase portraits for
each separately. In the language of chemical reaction dynamics $H_r$
corresponds to the “reactive mode”, $H_{b_1}$ and $H_{b_2}$ are “bath
modes”. The equilibrium point
$(q_1, p_1, q_2, p_2, q_3, p_3)=(0, 0, 0, 0, 0, 
0)$ is an index-1 saddle point for the full three DoF system on the zero
(total) energy surface.

In this system a trajectory “reacts” when its $q_1$ coordinate changes
sign. Therefore the surface $q_1=0$ is a [*dividing surface*]{} (DS) for
trajectories, separating “reactants” from “products” (i.e. before and
after reaction). This DS is a five dimensional surface in the six
dimensional phase space. We discuss its geometrical structure, in both
the phase space and in a fixed energy surface, in more detail.

First, we consider the “energetics” of the reaction. In order for $q_1$
to change sign we must have $H_r >0$. Also, it is clear from the form of
$H_{b_1}$ and $H_{b_2}$ that $H_{b_1} \geqslant 0$ and
$H_{b_2} \geqslant 0$. Therefore, reaction requires
$H = H_r + H_{b_1} + H_{b_2} > 0$.

The energy surface in which reaction occurs is given by

$$\frac{\lambda}{2} \left(p_1^2 - q_1^2 \right) + \frac{\omega_2}{2} \left(p_2^2 + q_2^2, \right)  + \frac{\omega_3}{2} \left(p_3^2 + q_3^2, \right)= H_r + H_{b_1} + H_{b_2} = H > 0, \quad H_r > 0, \, H_{b_1}, \, H_{b_2} \geqslant 0.
\label{3DoFES}$$

The intersection of the DS, $q_1=0$, with this energy surface is given
by:

$$\frac{\lambda}{2} \, p_1^2  + \frac{\omega_2}{2} \left(p_2^2 + q_2^2, \right) + \frac{\omega_3}{2} \left(p_3^2 + q_3^2, \right)= H_r + H_{b_1} + H_{b_2} = H > 0, \quad H_r > 0, \, H_{b_1}, \, H_{b_2} \geqslant 0.
\label{3DoDS}$$

This is the isoenergetic DS. It is a 4-sphere in the six dimensional
$(q_1, p_1, q_2, p_2, q_3, p_3)$ phase space. It has two hemispheres
that are DSs for the forward and backward reactions, respectively:

$$\frac{\lambda}{2} \, p_1^2  + \frac{\omega_2}{2} \left(p_2^2 + q_2^2, \right)  + \frac{\omega_3}{2} \left(p_3^2 + q_3^2, \right)= H_r + H_{b_1} + H_{b_2} = H > 0, \quad p_1 >0, \quad \mbox{forward DS},$$

$$\frac{\lambda}{2} \, p_1^2  + \frac{\omega_2}{2} \left(p_2^2 + q_2^2, \right)  + \frac{\omega_3}{2} \left(p_3^2 + q_3^2, \right)= H_r + H_{b_1} + H_{b_2} = H > 0, \quad p_1 <0, \quad \mbox{backward DS}.$$

Thes two hemispheres “meet” at $p_1 =0$:

$$\frac{\omega_2}{2} \left(p_2^2 + q_2^2, \right) + \frac{\omega_3}{2} \left(p_3^2 + q_3^2, \right)= H_{b_1} + H_{b_2} \geqslant  0,  \quad  \mbox{NHIM},
\label{eqn:sep_quad_ham3dof_nhim}$$

which is a [*normally hyperbolic invariant 3 sphere*]{}.
[*invariance*]{} follows from since for $q_1 = p_1 =0$ then
$\dot{q}_1 = \dot{p}_1 =0$. Thus $q_1$ and $p_1$ always remain zero, and
trajectories with these initial conditions remain on , i.e.,
$q_1 = p_1 =0$ is invariant. It is [*normally hyperbolic*]{} since the
directions normal to , i.e. $q_1-p_1$, are linearized saddle like
dynamics.

In order to understand the relationship between the NHIM and the index-1
saddle point we note that for $H_r =0, \, H_{b_1} =0$ and $H_{b_2} = 0$
the NHIM reduces to the point
$(q_1, p_1, q_2, p_2, q_3, p_3) = (0, 0, 0, 0, 0 ,0)$, which is the
index-1 saddle point on the energy surface $H_r + H_{b_1} + H_{b_2} =0$.
Therefore as the total energy is increased from $0$, with $H_{b_1}$
increasing from zero and/or $H_{b_2}$ increasing from zero, we see that
the NHIM “grows” from the index-1 saddle point on the zero energy
surface into an invariant 3 sphere. This shows how the “influence” of
the index-1 saddle point is carried to higher energy sufaces on which
the saddle point does not exist.

The stable and unstable manifolds of the NHIM are given the following
equations:

\begin{align}
\mathcal{W}^u\left( \mathcal{M}(h) \right) & = & \left \{ (q_1, p_1, q_2, p_2, q_3, p_3) \, \vert \, q_1=p_1, \, \frac{\omega_2}{2} \left(p_2^2 + q_2^2 \right) + \frac{\omega_3}{2} \left(p_3^2 + q_3^2 \right) 
= H_{b_1} + H_{b_2} > 0 \right \}, \label{eqn:quad_ham3dof_umani} \\
\mathcal{W}^s\left( \mathcal{M}(h) \right) & = & \left \{ (q_1, p_1, q_2, p_2, q_3, p_3) \, \vert q_1= -p_1, \, \frac{\omega_2}{2} \left(p_2^2 + q_2^2 \right) + \frac{\omega_3}{2} \left(p_3^2 + q_3^2 \right) 
=H_{b_1} + H_{b_2} > 0  \right \}, \label{eqn:quad_ham3dof_smani} 
\end{align}

and they are four dimensional on a fixed five dimensional energy
surface. The topology of the manifolds is the product of a line
($q_1 = p_1$ or $q_1 = -p_1$) with a 3 sphere
($\frac{\omega_2}{2} \left(p_2^2 + q_2^2 \right) + \frac{\omega_3}{2} \left(p_3^2 + q_3^2 \right) = H_{b_1} + H_{b_2} >0$),
that is $\mathbb{R} \times \mathbb{S}^3$ and sometimes such geometry is
referred to as “spherical cylinders”. Since the lines $q_1 = p_1$ and
$q_1 = -p_1$ correspond to the contour $H_r = 0$, in the six dimensional
phase space, these manifolds have energy
$H = H_r + H_{b_1}  + H_{b_2} = 0 + H_{b_1}  + H_{b_2} > 0$.

To understand the dynamics governed by the invariant manifolds of the
NHIM, let us choose an initial condition
$(q_1, p_1, q_2, p_2, q_3, p_3)$ on
$\mathcal{W}^u\left(\mathcal{M}(h)\right)$. Then as
$t \rightarrow +\infty$ the $(q_2, p_2, q_3, p_3)$ components of the
trajectory with this initial condition evolve quasiperiodically and the
$(q_1, p_1)$ components grow at an exponential rate. Similarly, if we
choose an initial condition on
$\mathcal{W}^s\left(\mathcal{M}(h)\right)$, then as
$t \rightarrow +\infty$ the $(q_2, p_2, q_3, p_3)$ components of the
trajectory evolve quasiperiodically and the $(q_1, p_1)$ components
decay to zero at an exponential rate as $t \rightarrow \infty$. In other
words, trajectories starting on
$\mathcal{W}^u \left(\mathcal{M}(h)\right)$ decay at an exponential rate
to the NHIM as $t \rightarrow - \infty$ and trajectories starting on
$\mathcal{W}^s\left(\mathcal{M}(h)\right)$ decay at an exponential rate
to the NHIM as $t \rightarrow  + \infty$.

#### Detecting NHIM and its manifolds

We consider isoenergetic two-dimensional surfaces parametrized by two
coordinates and compute the Lagrangian descriptor in a square domain of
size 2 units around the origin. We discretize the coordinates of the two
dimensional surface and pick constant values for three of the four
remaining coordinates, and use the total energy equation to solve for
the sixth coordinate. Due to the form of the Hamiltonian , obtaining the
coordinate from the constant energy condition is simply solving a
quadratic equation.

* Isoenergetic two-dimensional surface
parametrized by $(q_1, p_1)$  On the constant energy surface,
$H(q_1, p_1, q_2, p_2, q_3, p_3) = h$, we compute Lagrangian descriptor
on a two-dimensional surface parametrized by $(q_1, p_1)$ coordinates by
defining $$\begin{aligned}
U_{q_1p_1}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; | \; \right. & \left.   q_2 = 0,  p_2 = 0, q_3 = 0, \dot{q}_3 > 0 : \right. \nonumber \\
& \left.   p_3(q_1, p_1, q_2, p_2, q_3; h) > 0 \right\}\end{aligned}$$
where $$\begin{aligned}
p_3(q_1, p_1, q_2 = 0, p_2 = 0, q_3 = 0; h) = \sqrt{\frac{2}{\omega_3}\left( h - \frac{\lambda}{2}\left( p_1^2 - q_1^2 \right) \right)} \end{aligned}$$
The intersection of the two-dimensional surface $U_{q_1p_1}^+$ with the
NHIM  becomes

$$\begin{aligned}
\mathcal{M}(h) \cap U_{q_1p_1}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = 0, p_1 = 0, q_2 = 0,  p_2 = 0, q_3 = 0, \dot{q}_3 > 0 : \right. \nonumber \\ 
& \left.   p_3(q_1, p_1, q_2, p_2, q_3; h) > 0 \right\}.\end{aligned}$$

Thus, the NHIM is located at the origin $(0,0)$ and marked by a red
cross in the LD plot (Fig. \[fig:LD\_SQH\_3dof\_q1p1\]).

Next, the intersection of the two-dimensional surface with the unstable 
and stable manifolds  is given by

$$\begin{aligned}
\mathcal{W}^u(\mathcal{M}(h)) \cap U_{q_1p_1}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = p_1, q_2 = 0, p_2 = 0, q_3 = 0, \dot{q_3} > 0 : \right. \nonumber \\  
& \left.  p_3(q_1, p_1, q_2, p_2, q_3; h) > 0 \right\}, \\
\mathcal{W}^s(\mathcal{M}(h)) \cap U_{q_1p_1}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = -p_1, q_2 = 0, p_2 = 0, q_3 = 0,\dot{q_3} > 0 : \right. \nonumber \\  
& \left.  p_3(q_1, p_1, q_2, p_2, q_3;h) > 0 \right\},\end{aligned}$$

which are one-dimensional for a fixed energy, and represent lines
passing through the origin shown as dashed red (unstable) and white
(stable) lines, respectively, in Fig. \[fig:LD\_SQH\_3dof\_q1p1\]. The
only points of local minima in the LD plot
(Fig. \[fig:LD\_SQH\_3dof\_q1p1\]) also lie along the lines passing
through the origin and correspond to the manifolds of the NHIM.

* Isoenergetic two-dimensional surface
parametrized by $(q_2, p_2)$  On the constant energy surface,
$H(q_1, p_1, q_2, p_2, q_3, p_3) = h$, we compute Lagrangian descriptor
on a two-dimensional surface parametrized by $(q_2, p_2)$ coordinates by
defining $$\begin{aligned}
U_{q_2p_2}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; | \; \right. & \left.   q_1 = 0, q_3 = 0, p_3 = 0, \dot{q}_1 > 0 : \right. \nonumber \\
& \left.   p_1(q_1, q_2, p_2, q_3, p_3; h) > 0 \right\}\end{aligned}$$
where $$\begin{aligned}
p_1(q_1 = 0, q_2, p_2, q_3 = 0, p_3 = 0; h) = \sqrt{\frac{2}{\lambda}\left( h - \frac{\omega_2}{2}\left( p_2^2 + q_2^2 \right) \right)} \end{aligned}$$
The intersection of the two-dimensional surface $U_{q_2p_2}^+$ with the
NHIM  becomes

$$\begin{aligned}
\mathcal{M}(h) \cap U_{q_2p_2}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = 0, p_1 = 0, p_3 = 0, q_3 = 0, \dot{q}_3 > 0 : \right. \nonumber \\ 
& \left.   \frac{\omega_2}{2}\left( p_2^2 + q_2^2  \right) = h \right\}.\end{aligned}$$

Thus, the NHIM is the circle of radius $\sqrt{2h/\omega_2}$ and marked
by a dashed line in the LD plot (Fig. \[fig:LD\_SQH\_3dof\_q2p2\]).


Next, the intersection of the two-dimensional surface with the unstable 
and stable manifolds  is given by

$$\begin{aligned}
\mathcal{W}^u(\mathcal{M}(h)) \cap U_{q_2p_2}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = p_1, q_1 = 0, q_3 = 0, p_3 = 0, \dot{q_3} > 0 : \right. \nonumber \\  
& \left.  \frac{\omega_2}{2}\left( p_2^2 + q_2^2  \right) = h \right\}, \\
\mathcal{W}^s(\mathcal{M}(h)) \cap U_{q_2p_2}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = -p_1, q_1 = 0, q_3 = 0, p_3 = 0,\dot{q_3} > 0 : \right. \nonumber \\  
& \left.  \frac{\omega_2}{2}\left( p_2^2 + q_2^2  \right) = h \right\},\end{aligned}$$

which are one-dimensional for a fixed energy, and marked by dashed red
(unstable) and white (stable) lines, respectively, in
Fig. \[fig:LD\_SQH\_3dof\_q2p2\]. The only points of local minima in the
LD plot (Fig. \[fig:LD\_SQH\_3dof\_q2p2\]) also lie along the lines
passing through the origin and correspond to the manifolds of the NHIM.


![\label{fig:LD_SQH_3dof_q1p1}](./figures/lag_des_3dof_q20-q30-p20_500x500_E2e-01_q1p1.png){width="33.00000%"}\ ![\label{fig:LD_SQH_3dof_q2p2}](./figures/lag_des_3dof_q10-p30-q30_500x500_E2e-01_q2p2.png){width="33.00000%"}\ ![\label{fig:LD_SQH_3dof_q3p3}](./figures/lag_des_3dof_p10-q20-p20_500x500_E2e-01_q3p3.png){width="33.00000%"}
| -- |
|Lagrangian descriptor plot of the separable quadratic Hamiltonian vector field~\eqref{eqn:hameq3} on the isoenergetic two-dimensional surface (a) $U_{q_1p_1}^{+}$, (b) $U_{q_2p_2}^{+}$, (c) $U_{q_3p_3}^{+}$. The parameters used are $\lambda = \omega_2 = \omega_3 = 1.0$, $h = 0.2$, and $\tau = 10$.|




* Isoenergetic two-dimensional surface
parametrized by $(q_3, p_3)$  On the constant energy surface,
$H(q_1, p_1, q_2, p_2, q_3, p_3) = h$, we compute Lagrangian descriptor
on a two-dimensional surface parametrized by $(q_3, p_3)$ coordinates by
defining

$$\begin{aligned}
U_{q_3p_3}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; | \; \right. & \left.   p_1 = 0, q_2 = 0, p_2 = 0, \dot{p_1} > 0 : \right. \nonumber \\
& \left.   q_1(p_1, q_2, p_2, q_3, p_3; h) > 0 \right\}\end{aligned}$$

where

$$\begin{aligned}
q_1(p_1 = 0, q_2 = 0, p_2 = 0, q_3, p_3; h) = \sqrt{\frac{2}{\lambda}\left( \frac{\omega_3}{2}\left( p_3^2 + q_3^2 \right) - h \right)} \end{aligned}$$

The intersection of the two-dimensional surface $U_{q_3p_3}^+$ with the
NHIM  becomes

$$\begin{aligned}
\mathcal{M}(h) \cap U_{q_3p_3}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = 0, p_1 = 0, p_2 = 0, q_2 = 0, \dot{p_1} > 0 : \right. \nonumber \\ 
& \left.   \frac{\omega_3}{2}\left( p_3^2 + q_3^2  \right) = h \right\}.\end{aligned}$$

Thus, the NHIM is the circle of radius $\sqrt{2h/\omega_3}$ and marked
by a dashed line in the LD plot (Fig. \[fig:LD\_SQH\_3dof\_q3p3\]).

Next, the intersection of the two-dimensional surface with the unstable 
and stable manifolds  is given by

$$\begin{aligned}
\mathcal{W}^u(\mathcal{M}(h)) \cap U_{q_3p_3}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = p_1, p_1 = 0, q_2 = 0, p_2 = 0, \dot{p_1} > 0 : \right. \nonumber \\  
& \left.  \frac{\omega_3}{2}\left( p_3^2 + q_3^2  \right) = h \right\}, \\
\mathcal{W}^s(\mathcal{M}(h)) \cap U_{q_3p_3}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = -p_1, p_1 = 0, q_2 = 0, p_2 = 0,\dot{p_1} > 0 : \right. \nonumber \\  
& \left.  \frac{\omega_3}{2}\left( p_3^2 + q_3^2  \right) = h \right\},\end{aligned}$$

which are one-dimensional for a fixed energy and represent circles of
radius $\sqrt{2h/\omega_3}$, and marked by dashed red (unstable) and
white (stable) lines, respectively, in Fig. \[fig:LD\_SQH\_3dof\_q3p3\].
The only points of minima and singularity in the LD plot
(Fig. \[fig:LD\_SQH\_3dof\_q3p3\]) is along the circle and thus identify
the manifolds of the NHIM.

* Isoenergetic two-dimensional surface
parametrized by $(q_1, q_2)$  On the constant energy surface,
$H(q_1, p_1, q_2, p_2, q_3, p_3) = h$, we compute Lagrangian descriptor
on a two-dimensional surface parametrized by $(q_1, q_2)$ coordinates by
defining $$\begin{aligned}
U_{q_1q_2}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; | \; \right. & \left.   p_1 = 0, p_2 = 0, q_3 = 0, \dot{q}_3 > 0 : \right. \nonumber \\
& \left.   p_3(q_1, p_1, q_2, p_2, q_3; h) > 0 \right\}\end{aligned}$$
where $$\begin{aligned}
p_3(q_1, p_1 = 0, q_2, p_2 = 0, q_3 = 0; h) = \sqrt{\frac{2}{\omega_3}\left( h - \left( \frac{\omega_2}{2} q_2^2 - \frac{\lambda}{2} q_1^2 \right)  \right)} \end{aligned}$$
The intersection of the two-dimensional surface $U_{q_1q_2}^+$ with the
NHIM  is given by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{q_1q_2}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = 0, p_1 = 0, p_2 = 0, q_3 = 0, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2} q_2^2 + \frac{\omega_3}{2} p_3^2  = h : p_3(q_1, p_1, q_2, p_2, q_3; h) > 0 \right\}.\end{aligned}$$

Thus, on the isoenergetic two-dimensional surface $U_{q_1q_2}^+$, the
NHIM is the line $q_1 = 0$ and marked by a dashdot line in the LD plot
(Fig. \[fig:LD\_SQH\_3dof\_q1q2\]).

Next, the intersection of the unstable  and stable manifolds  with
two-dimensional surface is given by

$$\begin{aligned}
\mathcal{W}^u(\mathcal{M}(h)) \cap U_{q_1q_2}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = p_1, p_1 = 0, p_2 = 0, q_3 = 0, \right. \nonumber \\  
& \left.  \frac{\omega_2}{2} q_2^2 + \frac{\omega_3}{2} p_3^2  = h :  p_3(q_1, p_1, q_2, p_2, q_3; h) > 0 \right\}, \\
\mathcal{W}^s(\mathcal{M}(h)) \cap U_{q_1q_2}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = -p_1, p_1 = 0, p_2 = 0, q_3 = 0, \right. \nonumber \\  
& \left.  \frac{\omega_2}{2} q_2^2 + \frac{\omega_3}{2} p_3^2  = h :  p_3(q_1, p_1, q_2, p_2, q_3; h) > 0 \right\},\end{aligned}$$

which are one-dimensional for a fixed energy, and marked by dashed red
(unstable) and white (stable) lines, respectively, in
Fig. \[fig:LD\_SQH\_3dof\_q1q2\]. The only points of minima in the LD
plot (Fig. \[fig:LD\_SQH\_3dof\_q1q2\]) also lie along this line at
$q_1 = 0$ and identify the manifolds of the NHIM.

* Isoenergetic two-dimensional surface
parametrized by $(q_2, q_3)$  On the constant energy surface,
$H(q_1, p_1, q_2, p_2, q_3, p_3) = h$, we compute Lagrangian descriptor
on a two-dimensional surface parametrized by $(q_2, q_3)$ coordinates by
defining $$\begin{aligned}
U_{q_2q_3}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; | \; \right. & \left.   q_1 = 0, p_2 = 0, p_3 = 0, \dot{q}_1 \geqslant 0 : \right. \nonumber \\
& \left.   p_1(q_1, q_2, p_2, q_3, p_3; h) \geqslant 0 \right\}\end{aligned}$$
where $$\begin{aligned}
p_1(q_1 = 0, q_2, p_2 = 0, q_3, p_3 = 0; h) = \sqrt{\frac{2}{\lambda}\left( h - \left( \frac{\omega_2}{2} q_2^2 + \frac{\omega_3}{2} q_3^2 \right)  \right)} \end{aligned}$$
The intersection of the two-dimensional surface $U_{q_2q_3}^+$ with the
NHIM  is given by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{q_2q_3}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = 0, p_1 = 0, p_2 = 0, p_3 = 0, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2} q_2^2 + \frac{\omega_3}{2} q_3^2  = h, p_1(q_1, q_2, p_2, q_3, p_3; h) \geqslant 0  \right\}.\end{aligned}$$

Thus, on the isoenergetic two-dimensional surface $U_{q_2q_3}^+$, the
NHIM is an ellipse of semi-major axis $\sqrt{2h/\omega_2}$ and
semi-minor axis $\sqrt{2h/\omega_3}$ if $\omega_2 < \omega_3$, and
vice-versa, otherwise. In our case, $\omega_2 = \omega_3$, the
intersection of the NHIM with the isoenergetic two-dimensional surface
becomes a circle and marked by a dashdot line in the LD plot
(Fig. \[fig:LD\_SQH\_3dof\_q2q3\]).

Next, the intersection of the unstable  and stable manifolds  with the
isoenergetic two-dimensional surface is given by

$$\begin{aligned}
\mathcal{W}^u(\mathcal{M}(h)) \cap U_{q_2q_3}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = p_1, q_1 = 0, p_2 = 0, p_3 = 0, \right. \nonumber \\  
& \left.  \frac{\omega_2}{2} q_2^2 + \frac{\omega_3}{2} q_3^2  = h, p_1(q_1, q_2, p_2, q_3, p_3; h) \geqslant 0  \right\}, \\
\mathcal{W}^s(\mathcal{M}(h)) \cap U_{q_2q_3}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = -p_1, q_1 = 0, p_2 = 0, p_3 = 0, \right. \nonumber \\  
& \left.  \frac{\omega_2}{2} q_2^2 + \frac{\omega_3}{2} q_3^2  = h, p_1(q_1, q_2, p_2, q_3, p_3; h) \geqslant 0  \right\},\end{aligned}$$

which have the same geometry as the NHIM on the isoenergetic
two-dimensional surface, $U_{q_2q_3}$, and marked by dashed red
(unstable) and white (stable) lines, respectively, in
Fig. \[fig:LD\_SQH\_3dof\_q2q3\]. The only points of minima in the LD
plot (Fig. \[fig:LD\_SQH\_3dof\_q2q3\]) also lie along this circle and
identify the manifolds of the NHIM.

* Isoenergetic two-dimensional surface
parametrized by $(q_3, q_1)$  On the constant energy surface,
$H(q_1, p_1, q_2, p_2, q_3, p_3) = h$, we compute Lagrangian descriptor
on a two-dimensional surface parametrized by $(q_3, q_1)$ coordinates by
defining

$$\begin{aligned}
U_{q_3q_1}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; | \; \right. & \left.   q_2 = 0, p_2 = 0, p_3 = 0, \dot{q_1} > 0 : \right. \nonumber \\
& \left.   p_1(q_1, q_2, p_2, q_3, p_3; h) > 0 \right\}\end{aligned}$$

where

$$\begin{aligned}
p_1(q_1, q_2 = 0, p_2 = 0, q_3, p_3 = 0; h) = \sqrt{\frac{2}{\lambda}\left( h - \left( \frac{\omega_3}{2} q_3^2 - \frac{\lambda}{2} q_1^2 \right)  \right)} \end{aligned}$$

The intersection of the two-dimensional surface $U_{q_3q_1}^+$ with the
NHIM  is given by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{q_3q_1}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = 0, p_1 = 0, q_2 = 0, p_2 = 0, p_3 = 0, \right. \nonumber \\ 
& \left.  p_1(q_1, q_2, p_2, q_3, p_3; h) > 0, \frac{\omega_3}{2} q_3^2  = h \right\}.\end{aligned}$$

which represents two points $q_3 = \pm \sqrt{2h/\omega_3}$ on the line
$q_1 = 0$, and marked by a red cross in
Fig. \[fig:LD\_SQH\_3dof\_q3q1\]. These points are also identified by
the minima in the Lagrangian descriptor values as shown by
one-dimensional slices at constant $q_1$.

Next, the intersection of the unstable  and the stable manifold  with
the isoenergetic two-dimensional surface is given by

$$\begin{aligned}
\mathcal{W}^u(\mathcal{M}(h)) \cap U_{q_3q_1}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = p_1, q_2 = 0, p_2 = 0, p_3 = 0, \right. \nonumber \\  
& \left.  p_1(q_1, q_2, p_2, q_3, p_3; h) > 0, \frac{\omega_3}{2} q_3^2  = h \right\}, \\
\mathcal{W}^s(\mathcal{M}(h)) \cap U_{q_3q_1}^+ = \left\{ (q_1, p_1, q_2, p_2, q_3, p_3) \; \vert \; \right. & \left.   q_1 = -p_1, q_2 = 0, p_2 = 0, p_3 = 0, \right. \nonumber \\  
& \left.  p_1(q_1, q_2, p_2, q_3, p_3; h) > 0, \frac{\omega_3}{2} q_3^2  = h \right\},\end{aligned}$$

where each manifold represent lines parallel to $q_1$ axis. The unstable
manifold lies on the $q_1 > 0$ plane and the stable manifold lies on the
$q_1 < 0$ plane, and are marked by dashed red (unstable) and dashed
white (stable) lines in the Fig. \[fig:LD\_SQH\_3dof\_q3q1\]. These
manifolds are again identified by points of minima in the Lagrangian
descriptor values as shown by one-dimensional slices at constant $q_1$.


### Coupled quadratic Hamiltonian


To couple the coordinates in the separable quadratic Hamiltonian , we
introduce a linear transformation such that it satisfies the symplectic
condition .

Let us consider the symplectic transformation $$C =
\begin{pmatrix}
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 \\
-1 & 0 & 0 & 1 & 1 & 1 \\
0 & -1 & 0 & 1 & 1 & 1 \\
0 & 0 & -1 & 1 & 1 & 1 \\
\end{pmatrix}
\label{eqn:three_dof_C}$$

The change of coordinates is given by $$\begin{aligned}
q_1 &= p_x \\
p_1 &= -x + p_x + p_y + p_z \\
q_2 &= p_y \\
p_2 &= -y + p_x + p_y + p_z \\
q_3 &= p_z \\
p_3 &= -z + p_x + p_y + p_z
\end{aligned}$$ This transformation applied to  gives the Hamiltonian in
the transformed coordinates 
$$\begin{aligned}
\mathcal{H}(x, p_x, y, p_y, z, p_z) = \frac{\lambda}{2} \left[ (-x + p_x + p_y + p_z)^2 -  p_x^2 \right] + 
& 
\frac{\omega_2}{2} \left[ (-y + p_x + p_y + p_z)^2 + p_y^2 \right] \\ + &\frac{\omega_3}{2} \left[ (-z + p_x + p_y + p_z)^2 + p_z^2 \right] 
\label{eqn:ham_symp_3dof}
\end{aligned}$$ 
which gives the vector field 
\begin{align}
\dot{x} & = \frac{\partial \mathcal{H}}{\partial p_x} = \frac{\lambda}{2} \left[ 2(-x + p_x + 
p_y + p_z) - 2p_x \right] +  \frac{\omega_2}{2} \left[ 2(-y + p_x + p_y + p_z) \right] + 
\frac{\omega_3}{2} \left[ 2(-z + p_x + p_y + p_z) \right] \\
& = \lambda\left( -x + p_y + p_z \right) +  \omega_2 \left(-y + p_x + p_y + p_z \right) + 
\omega_3 \left(-z + p_x + p_y + p_z \right) \\
\dot{p_x} & = - \frac{\partial \mathcal{H}}{\partial x} \\ & = \lambda \left ( -x + p_x + p_y + p_z 
\right) \\
\dot{y} & = \frac{\partial \mathcal{H}}{\partial p_y} = \frac{\lambda}{2} \left[ 2(-x + p_x + 
p_y + p_z) \right] +  \frac{\omega_2}{2} \left[ 2p_y + 2(-y + p_x + p_y + p_z) \right] + 
\frac{\omega_3}{2} \left[ 2(-z + p_x + p_y + p_z) \right] \\
& = \lambda \left(-x + p_x + p_y + p_z \right) +  \omega_2 \left (-y + p_x + 2p_y + p_z \right) + 
\omega_3 \left(-z + p_x + p_y + p_z \right) \\
\dot{p_y} & = - \frac{\partial \mathcal{H}}{\partial y} \\ & = \omega_2 \left( -y + p_x + p_y + p_z 
\right) \\
\dot{z} & = \frac{\partial \mathcal{H}}{\partial p_z} = \frac{\lambda}{2} \left[ 2(-x + p_x + 
p_y + p_z) \right] +  \frac{\omega_2}{2} \left[ 2(-y + p_x + p_y + p_z) \right] + 
\frac{\omega_3}{2} \left[ 2p_z + 2(-z + p_x + p_y + p_z) \right] \\
& = \lambda \left(-x + p_x + p_y + p_z \right) +  \omega_2 \left(-y + p_x + p_y + p_z \right) + 
\omega_3 \left( -z + p_x + p_y + 2p_z \right) \\
\dot{p_z} & = - \frac{\partial \mathcal{H}}{\partial z} \\ & = \omega_3 \left( -z + p_x + p_y + p_z 
\right)
\label{eqn:eom_symp_3dof}
\end{align}
where the equilibrium point is at
$(0,0,0,0,0,0)$ and its total energy is $0$. The Jacobian at this
equilibrium point has eigenvalues
$\lambda, -\lambda, i\omega_2, -i \omega_2, i\omega_3, -i \omega_3$ (as
shown in the Appendix \[sect:linear\_symp\_transf\_jac\]) and is saddle
$\times$ center $\times$ center, thus index-1. The energy surface is a
five dimensional surface given by  in the six dimensional phase space.

In the decoupled (separable) system , the dividing surface is defined by
$q_1 = 0$ which becomes $p_x = 0$ in the coupled (non-separable)
coordinates. The dividing surface for a fixed energy
$\mathcal{H}(x,p_x,y,p_y,z,p_z) = h$ is given by
$${\rm DS} =  \frac{\lambda}{2} \left[ \left( -x + p_y + p_z \right)^2 \right] + \frac{\omega_2}{2} \left[ p_y^2 
+ \left( -y + p_y + p_z \right)^2  \right] + \frac{\omega_3}{2} \left[ p_z^2 + \left( -z + p_y + 
p_z \right)^2  \right] = h,$$ On this dividing surface, the NHIM is
defined by $p_1 = 0$ in the separable coordinates, which gives
$-x + p_x + p_y  + p_z= 0$ that is $p_y + p_z = x$ in the non-separable
coordinates, and can be expressed as
$$\mathcal{M}(h) = \frac{\omega_2}{2} \left[ (x - p_z)^2 + \left( x - y \right)^2 \right] + \frac{\omega_3}{2} \left[ p_z^2 
+ \left( x - z \right)^2 \right] = h.
\label{eqn:nhim_symp_3dof}$$ Next, in the coupled coordinates, the
stable and unstable manifolds of the NHIM are given by 

\begin{align}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) = \left\{  (x, y, z, p_x, p_y, p_z) \; | \; \right. & \left.   x = p_y + p_z,  \right. \nonumber \\ & \left.  \;  \frac{\omega_2}{2} \left( (-y + p_x + p_y + p_z )^2 + p_y^2 \right) + \right. \nonumber \\
&\mathrel{\phantom{=}} \left.  \frac{\omega_3}{2}\left( (-z + p_x + p_y + p_z )^2 + p_z^2 \right)  = h \right\}, \label{eqn:umani_nsqh_3dof}
\end{align}

\begin{align}
\mathcal{W}^{\rm s}(\mathcal{M}(h)) = \left\{  (x, y, z, p_x, p_y, p_z) \; | \; \right. & \left.   x = 2p_x + p_y + p_z, \right. \nonumber \\ &  \left.   \;  \frac{\omega_2}{2} \left( (-y + p_x + p_y + p_z )^2 + p_y^2 \right) + \right. \nonumber \\
&\mathrel{\phantom{=}} \left.  \frac{\omega_3}{2}\left( (-z + p_x + p_y + p_z )^2 + p_z^2 \right)  = h  \right\} \label{eqn:smani_nsqh_3dof}
\end{align}


As noted earlier, invertible linear symplectic transformations preserve
the normal hyperbolicity and the transversality of the Hamiltonian
vector field to the dividing surface. Hence the NHIM  still has a
geometry $\mathbb{S}^3$ and the invariant manifolds  and are also
$\mathbb{R} \times \mathbb{S}^3$.

#### Detecting NHIM and its manifolds

Now we illustrate the procedure for detecting NHIM and its stable and
unstable manifolds using features in Lagrangian descriptor on
isoenergetic two-dimensional surfaces for the non-separable quadratic
Hamiltonian .

* Isoenergetic two-dimensional surface
parametrized by $(x, p_x)$  On the constant energy surface,
$\mathcal{H}(x, p_x, y, p_y, z, p_z) = h$ , we compute the Lagrangian
descriptor on a two-dimensional surface parametrized by $(x, p_x)$
coordinates by defining

$$\begin{aligned}
U_{xp_x}^+ = \left\{ (x, p_x, y, p_y, z, p_z)  \; \vert \; \right. & \left.   y = 0, z = 0, p_y = 0, \right. \nonumber \\ 
& \left.  p_z(x, p_x, y, p_y, z; h) > 0 \, : \, \dot{z}(x, y, z, p_x, p_y, p_z) > 0 \right\}, 
\label{eqn:sos_symp_3dof_xpx}\end{aligned}$$

Thus, on the five-dimensional energy surface, the intersection of the
three-dimensional NHIM  with the two-dimensional surface is
zero-dimensional and given by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{xp_x}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   y = 0, z = 0, p_y = 0, p_x = 0, p_y + p_z = x, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2} \left( \left( x - p_z \right)^2 + x^2 \right) + \frac{\omega_3}{2}\left( p_z^2 + x^2 \right) = h, \right. \nonumber \\
& \left.   \dot{z}(x, p_x, y, p_y, z, p_z) > 0 : x > 0\right\}. 
\label{eqn:nhim_x_coord_3dof}\end{aligned}$$

which is a point on the line $p_x = 0$, and marked by a red cross in
Fig. \[fig:NSQH\_3dof\_xpx\]. This is also identified by the location of
minima in the Lagrangian descriptor values as shown by the
one-dimensional slice in Fig. \[fig:NSQH\_3dof\_xpx\] and agrees with
$x = \sqrt{2h/\left(\omega_2 + 2\omega_3 \right)}$ to within the grid
resolution.

Next, the intersection of the unstable  and stable  manifolds with the
isoenergetic two-dimensional surface  becomes

$$\begin{aligned}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) \cap U_{xp_x}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   y = 0, z = 0, p_y = 0, x = p_y + p_z, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2} \left( x + p_x \right)^2 + \frac{\omega_3}{2}\left( ( x + p_x )^2 + x^2 \right) = h, \right. \nonumber \\ 
& \left.  \dot{z}(x, p_x, y, p_y, z, p_z) > 0 : \right. \nonumber \\ 
& \left.  (\omega_2 + 2\omega_3)x + (\lambda + \omega_2 + \omega_3)p_x > 0 \right\}, \\
\mathcal{W}^{\rm s}(\mathcal{M}(h)) \cap U_{xp_x}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   y = 0, z = 0, p_y = 0, x = 2p_x + p_y + p_z, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2} \left( x - p_x \right)^2 + \frac{\omega_3}{2}\left( ( x - p_x )^2 + (x - 2p_x)^2 \right) = h, \right. \nonumber \\ 
& \left.  \dot{z}(x, p_x, y, p_y, z, p_z) > 0 : \right. \nonumber \\ 
& \left.  (\omega_2 + 2\omega_3)x + ( -\lambda - \omega_2 -3\omega_3)p_x > 0 \right\}\end{aligned}$$

which are one-dimensional and are shown as a dashed red line (unstable)
and a dashed white line (stable) in Fig. \[fig:NSQH\_3dof\_xpx\]. These
lines are also identified by the minima in the Lagrangian descriptor
values as shown by one-dimensional slice at constant $p_x$.

![\label{fig:NSQH_3dof_xpx}](./figures/linear_symp_trans_3dof_M400x400_E2e-01_xpx.png){width="33.00000%"}\ ![\label{fig:NSQH_3dof_ypy}](./figures/linear_symp_trans_3dof_M500x500_E2e-01_ypy.png){width="33.00000%"}\ ![\label{fig:NSQH_3dof_zpz}](./figures/linear_symp_trans_3dof_M500x500_E2e-01_zpz.png){width="33.00000%"}
Fig. 4: Lagrangian descriptor plot of the non-separable quadratic Hamiltonian vector field~\eqref{eqn:eom_symp_3dof} on the (a) $U_{xp_x}^+$, (b) $U_{yp_y}^+$, (c) $U_{zp_z}^+$. The parameters used are $\lambda = \omega_2 = \omega_3 = 1.0$, integration time of $\tau = 10$, and total energy $h = 0.2$.


* Isoenergetic two-dimensional surface
parametrized by $(y,p_y)$  On a constant energy surface,
$\mathcal{H}(x, p_x, y, p_y, z, p_z) = h$, we compute Lagrangian
descriptor on a two-dimensional surface by defining

$$\begin{aligned}
U_{yp_y}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   x = 0, \; z = 0, \; p_z = 0, \right. \nonumber \\ 
& \left.   p_x(x,y,p_y,z,p_z;h) > 0 : \dot{x}(x, p_x, y, p_y, z, p_z) > 0 \right\},
\label{eqn:sos_symp_3dof_ypy}\end{aligned}$$

The intersection of the NHIM  with the two-dimensional surface is given
by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{yp_y}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   x = 0, z = 0, p_z = 0, p_x = 0, p_y + p_z = x, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2} y^2 = h, \dot{x}(x, p_x, y, p_y, z, p_z) > 0 \, : \, y < 0 \right\}. 
\label{eqn:nhim_ypy_3dof}\end{aligned}$$

which is a point on the isoenergetic two-dimensional surface , and shown
as red cross in Fig. \[fig:NSQH\_3dof\_ypy\]. This point is also
identified by the minima in the Lagrangian descriptor values as shown by
the one-dimensional slices for constant $y$.


Next, the intersection of the unstable  and stable  manifolds with the
isoenergetic two-dimensional surface  becomes

$$\begin{aligned}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) \cap U_{ypy}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   x = 0, z = 0, p_z = 0 , x = p_y + p_z, \right. \nonumber \\ 
& \left.   \frac{\omega_2}{2}\left( -y + p_x \right)^2 + \frac{\omega_3}{2} p_x^2 = h, \right. \nonumber \\
& \left.    \dot{x}(x, p_x, y, p_y, z, p_z) > 0 \, : \, -\omega_2y + \lambda p_y > 0 \right\} \\ 
\mathcal{W}^{\rm s}(\mathcal{M}(h)) \cap U_{ypy}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   x = 0, z = 0, p_z = 0, x = 2p_x + p_y + p_z, \right. \nonumber \\  
& \left.   \frac{\omega_2}{2}\left( \left( -y + \frac{p_y}{2} \right)^2 + p_y^2 \right) + \frac{\omega_3}{2} \frac{p_y^2}{4} = h : \right. \nonumber \\
& \left.  \dot{x}(x, p_x, y, p_y, z, p_z) > 0 \, : \, -2\omega_2 y + (2\lambda + \omega_2 + \omega_3)p_y > 0 \right\}\end{aligned}$$

which are lines on the two-dimensional isoenergetic surface  and are
also identified by the minima in the Lagrangian descriptor values as
shown by the one-dimensional slices at constant $y$.

* Isoenergetic two-dimensional surface
parametrized by $(z,p_z)$  On a constant energy surface,
$\mathcal{H}(x, p_x, y, p_y, z, p_z) = h$ , we compute Lagrangian
descriptor on a two-dimensional surface by defining

\begin{align}
U_{zp_z}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   x = 0, \; y = 0, \; p_y = 0, \right. \nonumber \\ & \left.   p_x(x,y,p_y,z,p_z;h) > 0 : \dot{x}(x, p_x, y, p_y, z, p_z) > 0 \right\},
\label{eqn:sos_symp_3dof_zpz}
\end{align}

The intersection of the NHIM  with the two-dimensional surface is given
by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{zp_z}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   x = 0, y = 0, p_y = 0, p_x = 0, p_y + p_z = x, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2} p_z^2 + \frac{\omega_3}{2}\left( p_z^2 + z^2 \right) = h, \right. \nonumber \\
& \left.   \dot{x}(x, p_x, y, p_y, z, p_z) > 0 \, : \, z < 0 \right\}. 
\label{eqn:nhim_zpz_3dof}\end{aligned}$$

which is a point and shown as red cross in Fig. \[fig:NSQH\_3dof\_zpz\].
This point is also identified by the minima in the Lagrangian descriptor
values as shown by the one-dimensional slices for constant $z$.

Next, the intersection of the four-dimensional unstable  and stable 
manifolds with the isoenergetic two-dimensional surface  becomes


$$\begin{aligned}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) \cap U_{zpz}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   x = 0, y = 0, p_y = 0, p_z = 0, \right. \nonumber \\ 
& \left.   \frac{\omega_2}{2}p_x^2 + \frac{\omega_3}{2}\left( (-z + p_x)^2 \right) = h :  \right. \nonumber \\ 
& \left.    \dot{x}(x, p_x, y, p_y, z, p_z) > 0 \, : \, \omega_2 p_x + \omega_3( -z + p_x ) > 0 \right\} \\ 
\mathcal{W}^{\rm s}(\mathcal{M}(h)) \cap U_{zpz}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   x = 0,  y = 0, p_y = 0, p_z + 2p_x = 0, \right. \nonumber \\  
& \left.   \frac{\omega_2}{2}\left( \frac{p_z}{2} \right)^2 + \frac{\omega_3}{2}\left( \left(-z + \frac{p_z}{2}\right) ^2 + p_z^2\right) = h :  \right. \nonumber \\ 
& \left.    \dot{x}(x, p_x, y, p_y, z, p_z) > 0 \, : \, -\omega_3 z + (\lambda + \omega_2/2 + \omega_3/2)p_z > 0 \right\} \end{aligned}$$

which are one-dimensional and shown as dashed red (unstable) and dashed
white (stable) lines in Fig. \[fig:NSQH\_3dof\_zpz\]. These invariant
manifolds are also identified by minima in the Lagrangian descriptor
values as shown by the one-dimensional slice at constant $z$.

* Isoenergetic two-dimensional surface
parametrized by $(x,p_z)$  On a constant energy surface,
$\mathcal{H}(x, p_x, y, p_y, z, p_z) = h$, we compute the Lagrangian
descriptor on a two-dimensional surface by defining

$$\begin{aligned}
U_{xp_z}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   p_x = 0, y = 0, z = 0, \right. \nonumber \\ & \left.   p_y = p_y(x,p_x,y,p_y,z;h) : \dot{y}(x, p_x, y, p_y, z, p_z) > 0 \right\},
\label{eqn:sos_symp_3dof_xpz}\end{aligned}$$

The intersection of the NHIM  with the two-dimensional surface is given
by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{xp_z}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   p_x = 0, y = 0, z = 0, p_y + p_z = x, \right. \nonumber \\ 
& \left.  \frac{\omega_2}{2}\left( (x-p_z)^2 + x^2 \right) + \frac{\omega_3}{2}\left( p_z^2 + x^2 \right) = h, \right. \nonumber \\  
& \left.   \dot{y}(x, p_x, y, p_y, z, p_z) > 0 : \omega_2 (2x - p_z) + \omega_3x > 0 \right\}. 
\label{eqn:nhim_xpz_3dof}\end{aligned}$$

which is also identified by the minima in the Lagrangian descriptor
values as shown in the one-dimensional slices at constant $p_z$.

Next, the intersection of the unstable  and stable  manifolds with the
isoenergetic two-dimensional surface  becomes

$$\begin{aligned}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) \cap U_{xp_z}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   p_x = 0, y = 0, z = 0, x = p_y + p_z, \right. \nonumber \\ 
& \left.   \frac{\omega_2}{2} \left( x^2 + (x - p_z)^2 \right)  + \frac{\omega_3}{2} \left( x^2 + p_z^2 \right)= h, \right. \nonumber \\
& \left.   \dot{y}(x, p_x, y, p_y, z, p_z) > 0 : \omega_2 (2x - p_z) + \omega_3x > 0 \right\}, \\
\mathcal{W}^{\rm s}(\mathcal{M}(h)) \cap U_{xp_z}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   p_x = 0, y = 0, z = 0, x = 2p_x + p_y + p_z, \right. \nonumber \\  
& \left.   \frac{\omega_2}{2} \left( x^2 + (x - p_z)^2 \right)  + \frac{\omega_3}{2} \left( x^2 + p_z^2 \right)= h, \right. \nonumber \\
& \left.   \dot{y}(x, p_x, y, p_y, z, p_z) > 0 : \omega_2 (2x - p_z) + \omega_3x > 0 \right\}\end{aligned}$$

which project as the same one-dimensional curve as the NHIM and
identified by the minima in the Lagrangian descriptor values as shown by
the one-dimensional slices at constant $p_z$ in
Fig. \[fig:NSQH\_3dof\_xpz\].


![\label{fig:NSQH_3dof_xpz}}](./figures/linear_symp_trans_3dof_M500x500_E2e-01_xpz.png){width="33.00000%"}\ ![\label{fig:NSQH_3dof_ypz}](./figures/linear_symp_trans_3dof_M400x400_E2e-01_ypz.png){width="33.00000%"}\ ![\label{fig:NSQH_3dof_xz}](./figures/linear_symp_trans_3dof_M500x500_E2e-01_xz.png){width="33.00000%"}
Fig. 7: LD plot of the transformed Hamiltonian vector field~\eqref{eqn:eom_symp_3dof} on the  \protect\subref{fig:NSQH_3dof_ypz} $U_{xp_z}^+$,  \protect\subref{fig:NSQH_3dof_ypz} $U_{yp_z}^+$, \protect\subref{fig:NSQH_3dof_xz} $U_{xz}^+$. The parameters used are $\lambda = \omega_2 = \omega_3 = 1.0$, integration time of $\tau = 10$, and total energy $h = 0.2$.


* Isoenergetic two-dimensional surface
parametrized by $(y,p_z)$  Next, on a constant energy surface,
$\mathcal{H}(x, p_x, y, p_y, z, p_z) = h$, we compute Lagrangian
descriptor on a two-dimensional surface by defining

$$\begin{aligned}
U_{yp_z}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   x = 0, p_x = 0, z = 0, \right. \nonumber \\
& \left.   p_y(x,p_x,y,z,p_z;h) > 0 \, : \, \dot{x}(x,p_x,y,p_y,z,p_z) > 0 \right\}, 
\label{eqn:sos_symp_3dof_ypz}\end{aligned}$$

Thus, the intersection of the NHIM  with the two-dimensional surface is
given by

$$\begin{aligned}
\mathcal{M}(h) \cap U_{yp_z}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   x = 0, p_x = 0, z = 0, p_y + p_z = x, \right. \nonumber \\ 
& \left.   \frac{(\omega_2 + \omega_3)}{2}p_z^2 + \frac{\omega_2}{2}y^2 = h, \right. \nonumber \\  
& \left.   \dot{x}(x,p_x,y,p_y,z,p_z) > 0 \, : \, y < 0 \right\}. 
\label{eqn:nhim_ypz_3dof}\end{aligned}$$

which represents a portion of the ellipse and shown as the dashdot red
curve in Fig. \[fig:NSQH\_3dof\_ypz\]. As shown by the one-dimensional
slices for constant $p_z$, Eqn.  is identified by the minima in the
Lagrangian descriptor values.

Next, the intersection of the unstable  and stable  manifolds with the
isoenergetic two-dimensional surface  becomes

$$\begin{aligned}
\mathcal{W}^{\rm u}(\mathcal{M}(h)) \cap U_{ypz}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   x = 0, p_x = 0, z = 0, x = p_y + p_z, \right. \nonumber \\ 
& \left.   \frac{\omega_2}{2} y^2  + \frac{(\omega_2 + \omega_3)}{2} p_z^2 = h, \right. \nonumber \\
& \left.   \dot{x}(x,p_x,y,p_y,z,p_z) > 0 : y < 0 \right\}, \\
\mathcal{W}^{\rm s}(\mathcal{M}(h)) \cap U_{ypz}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   x = 0, p_x = 0, z = 0, x = 2p_x + p_y + p_z, \right. \nonumber \\  
& \left.   \frac{\omega_2}{2} y^2  + \frac{(\omega_2 + \omega_3)}{2} p_z^2 = h, \right. \nonumber \\
& \left.   \dot{x}(x,p_x,y,p_y,z,p_z) > 0 : y < 0 \right\}\end{aligned}$$

which are shown as red dashed line (unstable) and white dashed line
(stable) in Fig. \[fig:NSQH\_3dof\_ypz\]. We note here that the
intersection of the manifolds with the isoenergetic two-dimensional
surface  overlays on the NHIM’s intersection . These curves are also
identified by the minima in the Lagrangian descriptor values as shown by
the one-dimensional slice for constant $p_z$.

* Isoenergetic two-dimensional surface
parametrized by $(x,z)$  Next, on a constant energy surface,
$\mathcal{H}(x, p_x, y, p_y, z, p_z) = h$ , we compute the Lagrangian
descriptor on a two-dimensional surface by defining

$$\begin{aligned}
U_{xz}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   y = 0, p_x = 0, p_z = 0, \right. \nonumber \\ & \left.   p_y = p_y(x,p_x,y,z,p_z;h) : \dot{p_z}(x, p_x, y, p_y, z, p_z) > 0 \right\},
\label{eqn:sos_symp_3dof_xz}\end{aligned}$$

The intersection of the NHIM  with the two-dimensional surface is given
by

$$\begin{aligned}
{\rm NHIM} \cap U_{xz}^+ = \left\{ (x, p_x, y, p_y, z, p_z) \; \vert \; \right. & \left.   y = 0, p_x = 0, p_z = 0, p_y + p_z = x, \right. \nonumber \\ 
& \left.  \omega_2 x^2 + \frac{\omega_3}{2}\left( x - z \right)^2 = h, \right. \nonumber \\
& \left.   \dot{p_z}(x, p_x, y, p_y, z, p_z) > 0 : \omega_3(-z + x) > 0 \right\}. 
\label{eqn:nhim_xz_3dof}\end{aligned}$$

which represents a portion of an ellipse and is identified by the minima
in the Lagrangian descriptor values as shown by the one-dimensional
slices for constant $z$ in Fig. \[fig:NSQH\_3dof\_xz\].

Next, the intersection of the unstable  and stable  manifolds with the
isoenergetic two-dimensional surface  becomes

$$\begin{aligned}
W^{\rm u}({\rm NHIM}) \cap U_{xz}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   y = 0, p_x = 0, p_z = 0, x = p_y + p_z, \right. \nonumber \\ 
& \left.   \omega_2x^2  + \frac{\omega_3}{2} \left( -z + x \right)^2 = h, \right. \nonumber \\
& \left.   \dot{p_z}(x, p_x, y, p_y, z, p_z) > 0 : \omega_3(-z + x) > 0 \right\}, \\
W^{\rm s}({\rm NHIM}) \cap U_{xz}^+ = \left\{  (x, p_x, y, p_y, z, p_z) \; | \; \right. & \left.   y = 0, p_x = 0, p_z = 0, x = 2p_x + p_y + p_z, \right. \nonumber \\  
& \left.   \omega_2x^2  + \frac{\omega_3}{2} \left( -z + x \right)^2 = h, \right. \nonumber \\
& \left.   \dot{p_z}(x, p_x, y, p_y, z, p_z) > 0 : \omega_3(-z + x) > 0 \right\}\end{aligned}$$

which are also be identified by the minima in the Lagrangian descriptor
values as shown by the one-dimensional slice for constant $z$ in
Fig. \[fig:NSQH\_3dof\_xz\].

## Summary and Outlook {#sec:summ}


In this article, we assessed Lagrangian descriptor (LD) based detection
of high dimensional phase space structures, that is normally hyperbolic
invariant manifolds, their unstable and stable manifolds, in the two and
three degrees-of-freedom quadratic Hamiltonian systems with index-1
saddle. This is done using a sytematic comparison of numerical results
with analytical expressions for a comprehensive set of coordinates; see
Appendix \[sect:LD\_add\_slices\] for more pairs of coordinates that
support the results herein.. It is to be noted that in $N$ DoF system,
there are $2N(2N-1)(N-1)$ two-dimensional sections, so the detection
strategy should include inspecting a few pairs of coordinates. Based on
our investigation, a judicious choice of coordinate for the
two-dimensional surface is to use one of the reaction coordinates in
phase space, that is $q_1$ or $p_1$, and the other one from the
remaining $2N-1$ coordinates.

The closed-form analytical expressions of the NHIM and their manifolds
are precisely identified by the minima and singular features in the
Lagrangian descriptor values. This provides further numerical evidence
for the use of LDs in detecting high dimensional phase space structures
with a simple computation that can be implemented along with trajectory
integration. So at least for the form of Hamiltonian considered here,
one can rely with certainity on using Lagrangian descriptor to detect
the NHIM and its stable and unstable manifolds. Furthermore, this
detection approach has the potential to be combined with machine
learning type methods [@feldmaier_invariant_2019]. Future work on this
approach will include specific nonlinear systems that are inspired by
applications in celestial mechanics, ship dynamics, structural
mechanics, and chemical reaction dynamics.

## Acknowledgements

We acknowledge the support of EPSRC Grant No.  EP/P021123/1 and ONR
Grant No. N00014-01-1-0769. We would like to thank Vladim[í]{}r
Kraj[ň]{}[á]{}k for useful discussions and feedback.

