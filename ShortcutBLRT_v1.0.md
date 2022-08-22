Shortcut BLRT
================
Sam McGee
2022-08-22

# Estimating Power (![\hat{G}\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BG%7D_%7Ba%7D "\hat{G}_{a}")) for LCA / LPA Models

This RMarkdown uses Tekle, Gudicha, & Vermunt’s (2016) published method
for a shortcut BLRT as a means of accomplishing a power analysis. [The
complete paper can be found
here.](https://link.springer.com/article/10.1007/s11634-016-0251-0)

In short, traditional power analyses for latent class / profile models
(LCA / LPA) rely on Monte Carlo simulation to calculate the BLRT for
repeated samples under the *H<sub>1</sub>* population, yielding a “power
by proportion of *p* values” (*PPP*). However, this method is
computationally intensive, and requires researchers to repeat the entire
procedure for each adjustment to their proposed sample size. Thus,
determining the minimum necessary sample size is very time consuming,
and may be one reason why researchers often do not perform a power
analysis using the BLRT for LCA / LPA models.

In their paper, Tekle, Gudicha, & Vermunt (2016) propose, develop, and
test a “shortcut BLRT” method that is much less computationally
demanding, allowing researchers to adjust their sample size and
recalculate a BLRT relatively quickly, thus satisfying the needs of a
power analysis for minimum sample size in a study. This shortcut method
obtains a critical value (*CV*) from the distribution under
*H<sub>1</sub>*, approximating the distribution using Monte Carlo
simulation.

For an explanation of the Shortcut BLRT method, read the following
section. Otherwise, you may skip to the Guided Example section.

## Tekle, Gudicha, & Vermunt’s (2016) “Shortcut BLRT” Method

To accomplish this, we must first estimate the *CV* under
*H<sub>0</sub>*, or *CV<sub>0</sub>*. Then, we estimate power as the
proportion of likelihood ratio (*LR*) values exceeding *CV<sub>0</sub>*
in the Monte Carlo samples generated under *H<sub>1</sub>*. The detailed
steps are as follows, adapted directly from Tekle et al. (2016);

Given nominal significance
![a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a "a"),
the
![LR](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR "LR")
rejects the null hypothesis of
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
classes in favor of
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
if the observed *LR* statistic exceeds the
![CV_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;CV_0 "CV_0"),
which is proportional to
![C_a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_a "C_a")
and is the
![(1-a)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281-a%29 "(1-a)")th
quantile of the underlying
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
distribution.

> *Reject*
> ![H\_{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_%7B0%7D "H_{0}")
> *if*
> ![LR \> C_a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%20%3E%20C_a "LR > C_a")
> <span style="float:right">*(1)*</span>

The regularity conditions are violated in these comparisons, so we must
instead generate an empirical distribution
![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0")
under the
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
distribution. From here,
![C_a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_a "C_a")
can be estimated such that;

> ![P(LR \> C_a\\ \|\\ H_0) = P(LR \> C_a\\ \|\\ F_0) = {a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P%28LR%20%3E%20C_a%5C%20%7C%5C%20H_0%29%20%3D%20P%28LR%20%3E%20C_a%5C%20%7C%5C%20F_0%29%20%3D%20%7Ba%7D "P(LR > C_a\ |\ H_0) = P(LR > C_a\ |\ F_0) = {a}")
> <span style="float:right">*(2)*</span>

In order to identify the empirical distribution
![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0"),
we use the population parameters of
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0").
However, for a power analysis, we need to estimate these parameters from
a large data set generated under the
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
distribution, known as “exemplary data.” This piece is in-line with
traditional LRT methods, except that traditional methods use this Monte
Carlo procedure to estimate a
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
value for each data set, while we are estimating
![C_a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_a "C_a")
in this shortcut method.

Let
![\boldsymbol{Y}^{b} = (y^{b}\_{1},\\ ..., y^{b}\_T)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7BY%7D%5E%7Bb%7D%20%3D%20%28y%5E%7Bb%7D_%7B1%7D%2C%5C%20...%2C%20y%5E%7Bb%7D_T%29 "\boldsymbol{Y}^{b} = (y^{b}_{1},\ ..., y^{b}_T)")
be a random sample of size
![N](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N "N")
drawn from the LCA / LPA with
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
classes
![P(Y_i, \hat\Psi\_{K})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P%28Y_i%2C%20%5Chat%5CPsi_%7BK%7D%29 "P(Y_i, \hat\Psi_{K})"),
where
![\hat\Psi\_{K}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%5CPsi_%7BK%7D "\hat\Psi_{K}")
is the *ML* estimate under
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
that was based on the exemplary data set, which uses the LCA / LPA model
with
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
classes.

The *LR* statistic is computed for the replicate sample
![b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b "b"),
for
![b = 1,\\ ... B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b%20%3D%201%2C%5C%20...%20B "b = 1,\ ... B")
samples, generating
![LR^{b}\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bb%7D_0 "LR^{b}_0")
each time. These are then rearranged in order, such that
![LR^{1}\_0 \le LR^{2}\_0 \le ... \le LR^{B}\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7B1%7D_0%20%5Cle%20LR%5E%7B2%7D_0%20%5Cle%20...%20%5Cle%20LR%5E%7BB%7D_0 "LR^{1}_0 \le LR^{2}_0 \le ... \le LR^{B}_0"),
from which we can obtain the estimate for
![C\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_%7Ba%7D "C_{a}")
as the quantile at
![\[B(1-{a})\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5BB%281-%7Ba%7D%29%5D "[B(1-{a})]"),
or the
![h](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h "h")th
quantile in the ordered statistic of the bootstrap *LR* under
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0").
In mathematical terms;

> ![\hat{C}\_{a} = Q\_{\[B(1 - {a})\]} = Q\_{\[h\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BC%7D_%7Ba%7D%20%3D%20Q_%7B%5BB%281%20-%20%7Ba%7D%29%5D%7D%20%3D%20Q_%7B%5Bh%5D%7D "\hat{C}_{a} = Q_{[B(1 - {a})]} = Q_{[h]}")
> <span style="float:right">*(3)*</span>

After estimating
![C\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_%7Ba%7D "C_{a}"),
power can be computed given the significance level
![{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7Ba%7D "{a}"),
the population values under
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1"),
and the data design factors including sample size. We estimate the power
![G\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_%7Ba%7D "G_{a}")
for the LRT as the probability that the observed value of the *LR*
exceeds the *CV* or
![C\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_%7Ba%7D "C_{a}")
given that alternative hypothesis model
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
and it’s accompanying distribution
![F_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_1 "F_1")
holds. Again, in mathematical terms;

> ![G\_{a} = P(LR \> C\_{a}\\ \|\\ F_1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_%7Ba%7D%20%3D%20P%28LR%20%3E%20C_%7Ba%7D%5C%20%7C%5C%20F_1%29 "G_{a} = P(LR > C_{a}\ |\ F_1)")
> <span style="float:right">*(4)*</span>

We use the previously estimated
![\hat{C}\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BC%7D_%7Ba%7D "\hat{C}_{a}")
and construct the empirical distribution of the *LR* statistic under
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1"),
again using Monte Carlo simulation. Specifically, given the hypothesized
LCA / LPA
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
parameter values for the
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
alternative, we generate
![M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M "M")
random samples of size
![N](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N "N")
from the population defined by the alternate hypothesis. On each drawn
sample, fit both the
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
LCA / LPA
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
and the
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
LCA / LPA
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
model, computing the *LR* statistic
![LR^{m}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bm%7D_1 "LR^{m}_1")
for each sample. This same statistic is then ordered,
![\\{LR^{1}\_1,\\ LR^{2}\_1,\\ LR^{3}\_1,\\ ...,\\ LR^{M}\_1\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7BLR%5E%7B1%7D_1%2C%5C%20LR%5E%7B2%7D_1%2C%5C%20LR%5E%7B3%7D_1%2C%5C%20...%2C%5C%20LR%5E%7BM%7D_1%5C%7D "\{LR^{1}_1,\ LR^{2}_1,\ LR^{3}_1,\ ...,\ LR^{M}_1\}"),
forming the empirical distribution of the *LR* under the alternative
hypothesis. Based on this, the power
![G\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_%7Ba%7D "G_{a}")
is computed as;

> ![\hat{G}\_{a} = \frac{1}{M} \sum\_{m=1}^{M} I\_{\[LR^{m}\_{1} \> \hat{C}\_{a}\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BG%7D_%7Ba%7D%20%3D%20%5Cfrac%7B1%7D%7BM%7D%20%5Csum_%7Bm%3D1%7D%5E%7BM%7D%20I_%7B%5BLR%5E%7Bm%7D_%7B1%7D%20%3E%20%5Chat%7BC%7D_%7Ba%7D%5D%7D "\hat{G}_{a} = \frac{1}{M} \sum_{m=1}^{M} I_{[LR^{m}_{1} > \hat{C}_{a}]}")
> <span style="float:right">*(5)*</span>

where
![I\_{\[h\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I_%7B%5Bh%5D%7D "I_{[h]}")
is an indicator function that returns 1 if
![h](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h "h"),
defined here as
![\[LR^{m}\_1 \> C\_{a}\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5BLR%5E%7Bm%7D_1%20%3E%20C_%7Ba%7D%5D "[LR^{m}_1 > C_{a}]"),
is true and returns 0 if
![h](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h "h")
is false. Power,
![\hat{G}\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BG%7D_%7Ba%7D "\hat{G}_{a}")
is directly estimated using the proportion where
![h](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h "h")
is true and the alternative model
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
is correctly identified.

The research goal is often to determine the smallest number of subjects
required to achieve sufficient power,
![G_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_0 "G_0").
We simulate our samples from the population parameters under the
alternative hypothesis using the shortcut method here. This is
determined as;

> ![\\{\min(n):\\ \hat{G}\_{a}(n) \> G\_{0} \\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B%5Cmin%28n%29%3A%5C%20%5Chat%7BG%7D_%7Ba%7D%28n%29%20%3E%20G_%7B0%7D%20%5C%7D "\{\min(n):\ \hat{G}_{a}(n) > G_{0} \}")
> <span style="float:right">*(6)*</span>

We can then evaluate
![\hat{G}\_{a}(n) = \frac{1}{M}\sum^{M}\_{m=1}I\_{\[LR^{m}\_1(n)\\ \>\\ \hat{C}\_{a}(n)\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BG%7D_%7Ba%7D%28n%29%20%3D%20%5Cfrac%7B1%7D%7BM%7D%5Csum%5E%7BM%7D_%7Bm%3D1%7DI_%7B%5BLR%5E%7Bm%7D_1%28n%29%5C%20%3E%5C%20%5Chat%7BC%7D_%7Ba%7D%28n%29%5D%7D "\hat{G}_{a}(n) = \frac{1}{M}\sum^{M}_{m=1}I_{[LR^{m}_1(n)\ >\ \hat{C}_{a}(n)]}")
where
![I\_{\[h\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I_%7B%5Bh%5D%7D "I_{[h]}")
is the same indicator function as in *(5)*,
![LR^{m}\_1(n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bm%7D_1%28n%29 "LR^{m}_1(n)")
is the *LR* statistic at the
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")th
Monte Carlo sample of size
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
from the population model under
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1"),
and
![\hat{C}\_{a}(n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BC%7D_%7Ba%7D%28n%29 "\hat{C}_{a}(n)")
is the *CV* computed based on
![B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B "B")
samples from the population model under
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0").

Here, a series of trials are performed to solve *(6)*. Using an
exemplary data set generated under the population model for
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1"),
we obtain the parameter estimate for the
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
model. The *CV* is computed based on
![B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B "B")
independent samples of size
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
drawn from
![P(\boldsymbol{Y}\_{i},\\ \hat{\psi}\_{K + 1})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P%28%5Cboldsymbol%7BY%7D_%7Bi%7D%2C%5C%20%5Chat%7B%5Cpsi%7D_%7BK%20%2B%201%7D%29 "P(\boldsymbol{Y}_{i},\ \hat{\psi}_{K + 1})")
, the alternative model. Next,
![M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M "M")
independent samples of size
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
are drawn from the same
![P(\boldsymbol{Y}\_{i},\\ \hat{\psi}\_{K + 1})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P%28%5Cboldsymbol%7BY%7D_%7Bi%7D%2C%5C%20%5Chat%7B%5Cpsi%7D_%7BK%20%2B%201%7D%29 "P(\boldsymbol{Y}_{i},\ \hat{\psi}_{K + 1})")
population model, from which we evaluate the test statistic
![LR^{m}\_{1}(n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bm%7D_%7B1%7D%28n%29 "LR^{m}_{1}(n)")
for each sample for $m = 1,\\ …,\\ M$. Tekle, Gudicha, & Vermunt (2016)
recommend a linear search algorithm, where a smaller
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
value is needed if $\hat{G}\_{a}(n) \> G\_{0} \\}$, otherwise a larger
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
is still required for sufficient power.

## Guided Example from Tekle, Gudicha, & Vermunt (2016)

We follow Tekle, Gudicha, & Vermunt’s (2016) example, using an LCA model
with
![K = 3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%3D%203 "K = 3")
classes,
![T = 10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T%20%3D%2010 "T = 10")
indicator variables, with moderate separation for class-specific
responses
![(0.8)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280.8%29 "(0.8)").
We follow this model to demonstrate that the Shortcut BLRT method
performs well under these conditions and identifies a small
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
sample size as sufficient for
![{G} = 0.8](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7BG%7D%20%3D%200.8 "{G} = 0.8")
with Type I error rate fixed to ${a} = 0.05}, both of which are
traditionally considered sufficient for research purposes.

Using the
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
class parameters established in their example, the response
probabilities are set to
![\theta\_{kt} = 0.8](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bkt%7D%20%3D%200.8 "\theta_{kt} = 0.8")
in class 1,
![\theta\_{kt} = 0.2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bkt%7D%20%3D%200.2 "\theta_{kt} = 0.2")
in class 2,
![\theta\_{kt} = 0.8](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bkt%7D%20%3D%200.8 "\theta_{kt} = 0.8")
for the first half of
![T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T "T")
indicators and
![\theta\_{kt} = 0.2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bkt%7D%20%3D%200.2 "\theta_{kt} = 0.2")
for the other half in class 3, and to
![\theta\_{kt} = 0.2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bkt%7D%20%3D%200.2 "\theta_{kt} = 0.2")
for the first half of indicators and
![\theta\_{kt} = 0.8](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bkt%7D%20%3D%200.8 "\theta_{kt} = 0.8")
for the other half in class 4. We also specify class sizes to be equal,
cycling through 1:4 repeatedly. These can be changed in the embedded
code by alternating the commented lines below the p F0.ID.

# Establishing ![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0")

``` r
# Define K and K + 1
K0 = 3
K1 = 4

# Define T
n_indicators = 10

# Create F0 data set using model parameters
F0.Data <- data.frame(matrix(nrow=1000,ncol=(n_indicators + 1))) # ncol includes T and a column for class ID

F0.ID <- rep_len(1:K1, 1000) # generate equal class sizes
# F0.ID <- sample(c(1:4), 1000, replace=TRUE, prob=c(0.4, 0.3, 0.2, 0.1)) # generate unequal class sizes

# Estimate responses based on class probabilities defined by Tekle, Gudicha, & Vermunt
for (i in 1:nrow(F0.Data)) {
  F0.Data[i,11] <- F0.ID[i]
  if(F0.Data[i,11] == 1) {
    for(j in 1:10) {
      F0.Data[i,j] <- rbinom(1, 4, 0.8)
    }
  }
  else if(F0.Data[i,11] == 2) {
    for(j in 1:10) {
      F0.Data[i,j] <- rbinom(1, 4, 0.2)
    }
  }
  else if(F0.Data[i,11] == 3) {
    for(j in 1:5) {
      F0.Data[i,j] <- rbinom(1, 4, 0.8)
    }
    for(j in 6:10) {
      F0.Data[i,j] <- rbinom(1, 4, 0.2)
    }
  }
  else if(F0.Data[i,11] == 4) {
    for(j in 1:5) {
      F0.Data[i,j] <- rbinom(1, 4, 0.2)
    }
    for(j in 6:10) {
      F0.Data[i,j] <- rbinom(1, 4, 0.8)
    }
  }
}
# Assign ID column name for clarity
colnames(F0.Data)[11] <- 'ID'
```

Now that the exemplary data set has been generated for
![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0"),
we can estimate the model parameters for both
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
and
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1").

``` r
# Fit H0 model using tidyLPA
H0.Model <- F0.Data %>%
  # select indicator variables in model
  select(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10) %>%
  # estimate number of profiles
  estimate_profiles(n_profiles = K0,
                    # adjust variances and covariances if desired; default is Model 1 for parsimony
                    models = 1)

# Get H0 estimates
H0.Estimates <- as.data.table(get_estimates(H0.Model))

# Fit H1 model
H1.Model <- F0.Data %>%
  # select indicator variables in model
  select(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10) %>%
  # estimate number of profiles
  estimate_profiles(n_profiles = K1,
                    # adjust variances and covariances if desired; default is Model 1 for parsimony
                    models = 1)

# Get H1 estimates
H1.Estimates <- as.data.table(get_estimates(H1.Model))
```

With our
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
and
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
model estimates from the sample obtained under the alternative
hypothesis, we now set up our desired sample size for this trial, along
with the length for
![B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B "B")
and
![M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M "M").
In Tekle, Gudicha, & Vermunt’s (2016) paper, they set these to 500, but
this can be adjusted here as needed. Additionally, we set up holding
vectors to retain the *-2 Log likelihood (-2LL)* fit statistic
information from each iteration, and a vector with the results of
![I\_{\[h\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I_%7B%5Bh%5D%7D "I_{[h]}")
in *(5)*.

``` r
# Set n, B, and M
n = 200
B = 10
M = 10

# Set alpha level
alpha = 0.05

# Vector with indicator function results
IndicatorM <- matrix(0, nrow=1, ncol=M)

# Holding vectors for model B and M fit results of each iteration
BFit <- rep(0, B)

MFit <- rep(0, M)
```

We generate the list
![\\{LR^{1}\_1,\\ LR^{2}\_1,\\ LR^{3}\_1,\\ ...,\\ LR^{B}\_1\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7BLR%5E%7B1%7D_1%2C%5C%20LR%5E%7B2%7D_1%2C%5C%20LR%5E%7B3%7D_1%2C%5C%20...%2C%5C%20LR%5E%7BB%7D_1%5C%7D "\{LR^{1}_1,\ LR^{2}_1,\ LR^{3}_1,\ ...,\ LR^{B}_1\}")
using the
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
parameters under the
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
population model to create an ordered
![LR^{b}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bb%7D_1 "LR^{b}_1")
statistic, in equation *(3)*. Additionally, we generate
![\\{LR^{1}\_1,\\ LR^{2}\_1,\\ LR^{3}\_1,\\ ...,\\ LR^{M}\_1\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7BLR%5E%7B1%7D_1%2C%5C%20LR%5E%7B2%7D_1%2C%5C%20LR%5E%7B3%7D_1%2C%5C%20...%2C%5C%20LR%5E%7BM%7D_1%5C%7D "\{LR^{1}_1,\ LR^{2}_1,\ LR^{3}_1,\ ...,\ LR^{M}_1\}")
from the
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
model to the same random sample drawn from
![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0"),
and save the resulting
![LR^{m}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bm%7D_1 "LR^{m}_1")
statistic, which will be used in equation *(5)*.

``` r
# Calculate the LR for H0 for each iteration B
for(i in 1:B) {
  # Draw a random sample n from the exemplary data
  Sampled.F <- sample_n(F0.Data, n)
  
  Temp.BFit <- Sampled.F %>%
    # Fit the H0 model to the sampled data
    select(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10) %>%
    estimate_profiles(n_profiles = K0, models = 1)
  
  Temp.MFit <- Sampled.F %>%
    # Fit the H1 model to the sampled data
    select(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10) %>%
    estimate_profiles(n_profiles = K1, models = 1)
  
  # Save this iteration's -2LL to BFit and MFit
  BFit[i] <- format(Temp.BFit$model_1_class_3$fit[3], scientific = F)
  MFit[i] <- format(Temp.MFit$model_1_class_4$fit[3], scientific = F)
}
```

Next, we calculate the
![\hat{C}\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BC%7D_%7Ba%7D "\hat{C}_{a}")
at the
![Q\_{\[h\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Q_%7B%5Bh%5D%7D "Q_{[h]}")th
position, using the ordered
![LR^{b}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bb%7D_1 "LR^{b}_1")
statistic.

``` r
# Sort the resulting LR vector from smallest to largest
LR0 <- as.vector(sort(BFit, na.last = NA, decreasing = TRUE))

# Identify C_a for these data
C_a <- as.numeric(LR0[round(length(LR0)*(1-alpha))])
```

Now, using our indicator function in *(5)*, we can calculate
![\hat{G}\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BG%7D_%7Ba%7D "\hat{G}_{a}")
as the proportion of
![LR^{m}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bm%7D_1 "LR^{m}_1")
from the
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
model that are preferred over
![\hat{C}\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BC%7D_%7Ba%7D "\hat{C}_{a}")
that was generated from
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0"),
all fit to our exemplary data
![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0")
which were generated with the
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
population parameters.

``` r
for(i in 1:length(MFit)){
  if(MFit[i] < C_a) { # the evaluation direction is reversed due to the nature of -2LL
    IndicatorM[i] <- 1
  }
}

Ga = (1/M)*(sum(IndicatorM))
```

Estimated power = 1

# Simulated Class Parameters & Sample

Here we can define the LCA model parameters directly and generate the
necessary
![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0")
data set to define the
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
class model, the
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
model, and determine the minimum necessary sample size
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
for sufficient power.

First, define
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K"),
which in turn defines
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1").
Additionally, we define
![T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T "T"),
the number of indicators, as well as the responses per indicator.
Following this,
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
is defined as the class separation, which in turn allows us to randomly
generate
![\boldsymbol{Y}\_k = (m^{1}\_{k},\\ ...,\\ m^{T}\_{k})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7BY%7D_k%20%3D%20%28m%5E%7B1%7D_%7Bk%7D%2C%5C%20...%2C%5C%20m%5E%7BT%7D_%7Bk%7D%29 "\boldsymbol{Y}_k = (m^{1}_{k},\ ...,\ m^{T}_{k})"),
where
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")
is an indicator. Using these class means,
![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0")
can then be generated.

``` r
# Define K
K0 = 3
K1 = K0 + 1

# Define number of indicators
n_indicators = 10
item_responses = 4

# Define class-separation using random mean generation
theta = 0.7 # {0.7 = low, 0.8 = moderate, 0.9 = high}

# Define item-response estimates by class w/ theta separation
theta.k <- data.frame(matrix(0,nrow = K1, ncol = n_indicators))
# Define the mean differences between classes based on theta
#   (Alternatively, we could manually set these values)
theta.scaled <- seq(from = 0.1,to = theta, length.out = K1)
for(i in 1:nrow(theta.k)) {
  for(j in 1:ncol(theta.k)) {
    theta.k[i,j] <- format(rbinom(1,item_responses,theta.scaled[i]), scientific = F)
  }
}

# Random F0.ID generation
F0.ID <- sample(1:K1, 1000, replace = TRUE, prob = NULL)
# Fixed F0.ID generation
# F0.ID <- rep_len(1:K1, 1000)

F0.Data <- data.frame(matrix(nrow=1000,ncol=(n_indicators+1)))
for(i in 1:nrow(F0.Data)) {
  F0.Data[i,ncol(F0.Data)] <- F0.ID[i]
  for(j in 1:K1) {
    if(F0.Data[i,ncol(F0.Data)] == j) {
      for(l in 1:n_indicators) {
        temp.mean <- (as.numeric(theta.k[j,l])/item_responses)
        F0.Data[i,l] <- rbinom(n = 1,size = item_responses,prob = temp.mean)
      }
    }
  }
}
colnames(F0.Data)[length(F0.Data)] <- 'ID'
```

To be slightly more user-friendly, we will define
![n, B, M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%2C%20B%2C%20M "n, B, M")
and
![a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a "a")
here, to be used below.

``` r
n = 50
B = 10
M = 10
alpha = 0.05
```

The Shortcut BLRT is then repeated as before. Instead of being broken
apart by section, all the previous steps are in one block of {r} code
here, but the process is identical.

``` r
# Save colnames as list of names
col_list <- colnames(F0.Data)[1:(length(F0.Data)-1)]

# Fit H0 model using tidyLPA
H0.Model <- F0.Data %>%
  # select indicator variables in model
  select(col_list[1:length(col_list)]) %>%
  # estimate number of profiles
  estimate_profiles(n_profiles = K0,
                    # adjust variances and covariances if desired; default is Model 1 for parsimony
                    models = 1)

# Get H0 estimates
H0.Estimates <- as.data.table(get_estimates(H0.Model))

# Fit H1 model
H1.Model <- F0.Data %>%
  # select indicator variables in model
  select(col_list[1:length(col_list)]) %>%
  # estimate number of profiles
  estimate_profiles(n_profiles = K1,
                    # adjust variances and covariances if desired; default is Model 1 for parsimony
                    models = 1)

# Get H1 estimates
H1.Estimates <- as.data.table(get_estimates(H1.Model))

# Vector with indicator function results
IndicatorM <- matrix(0, nrow=1, ncol=M)

# Holding vectors for model B and M fit results of each iteration
BFit <- rep(0, B)

MFit <- rep(0, M)

# Calculate the LR for H0 for each iteration B
for(i in 1:B) {
  # Draw a random sample n from the exemplary data
  Sampled.F <- sample_n(F0.Data, n)
  
  Temp.BFit <- Sampled.F %>%
    # Fit the H0 model to the sampled data
    select(col_list[1:length(col_list)]) %>%
    estimate_profiles(n_profiles = K0, models = 1)
  
  Temp.MFit <- Sampled.F %>%
    # Fit the H1 model to the sampled data
    select(col_list[1:length(col_list)]) %>%
    estimate_profiles(n_profiles = K1, models = 1)
  
  # Save this iteration's -2LL to BFit and MFit
  BFit[i] <- format(Temp.BFit$model_1_class_3$fit[3], scientific = F)
  MFit[i] <- format(Temp.MFit$model_1_class_4$fit[3], scientific = F)
}

# Sort the resulting LR vector from smallest to largest
LR0 <- as.vector(sort(BFit, na.last = NA, decreasing = TRUE))

# Identify C_a for these data
C_a <- as.numeric(LR0[round(length(LR0)*(1-alpha))])

for(i in 1:length(MFit)){
  if(MFit[i] < C_a) { # the evaluation direction is reversed due to the nature of -2LL
    IndicatorM[i] <- 1
  }
}

Ga = (1/M)*(sum(IndicatorM))
```

Estimated power = 0.8

# Recommendations and Implementation

In order to use this RMarkdown for your own data purposes, you will need
to adjust these parameters to fit the model you define. It should be
noted that while this example uses an LCA specifically, LPA models can
also be applied here, but follows the recommendations that more
indicator
![T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T "T")
variables for these models generally perform better, and the need for
higher class separation is more apparent.

My own experience using this for a real-world model involved determining
my
![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
model with an existing data set, generating the model and parameters for
the
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
and
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
models using these data, and using the
![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
parameters (including class size and item estimates) to generate
![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0").

To accomplish this, I recommend creating a tidyLPA object for your
desired model, including
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
classes. Then, create a second tidyLPA object for
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
and determine the class size probabilities. **Please note;** you should
always double-check your class sizes and probability per class. If the
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
model has a class with \<1% chance of occurring, even a large data set
will struggle to identify this class, and it is recommended to consider
whether such a low-occurring class is “real,” or if the present sample
of data is appropriate for investigating such a class. Researcher
discretion is highly advised here.

Once the tidyLPA objects for the
![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
and
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
models are identified, you can generate an
![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0")
by sampling ID assignment
1:![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1"),
using rnorm(), your model-specified class means and variances, and some
for() loops. You should wind up with an exemplary data set of
![n = 1000](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%3D%201000 "n = 1000"),
**based on the class probabilities and parameters under your**
![K + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%20%2B%201 "K + 1")
**model.**

Following this, repeat as above;

-   Sample from the exemplary data your desired size
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n").

-   Re-estimate the
    ![H_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_0 "H_0")
    and
    ![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
    model parameters using the sampled
    ![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0")
    exemplary data.

-   Store both *-2LL* estimates from each model into their respective
    ![LR^{b}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bb%7D_1 "LR^{b}_1")
    and
    ![LR^{m}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bm%7D_1 "LR^{m}_1")
    vectors.

-   Repeat the previous 3 steps
    ![M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M "M")
    number of times (typical recommendations are \~500).

-   Sort the
    ![LR^{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bb%7D "LR^{b}")
    vector, and identify your
    ![\hat{C}\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BC%7D_%7Ba%7D "\hat{C}_{a}")
    value at the
    ![Q\_{\[B(1 - {a})\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Q_%7B%5BB%281%20-%20%7Ba%7D%29%5D%7D "Q_{[B(1 - {a})]}")th
    position.

-   Using the indicator function
    ![I\_{\[h\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I_%7B%5Bh%5D%7D "I_{[h]}"),
    determine the proportion of
    ![LR^{m}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bm%7D "LR^{m}")
    that exceed your critical value,
    ![\hat{C}\_{a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BC%7D_%7Ba%7D "\hat{C}_{a}").

-   ![\hat{G}\_{a}(n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BG%7D_%7Ba%7D%28n%29 "\hat{G}_{a}(n)")
    is then estimated as the proportion of
    ![LR^{m}(n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;LR%5E%7Bm%7D%28n%29 "LR^{m}(n)")
    that correctly identify the
    ![H_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_1 "H_1")
    model from which the
    ![F_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0 "F_0")
    exemplary data were generated from.
