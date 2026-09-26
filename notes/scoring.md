There are 6 metrics which have raw values based upon the our perturbation profile and the actual experimental result. Each of the 6 raw metrics are scaled within context generating 6 scaled metrics per context which are averaged into 1 scaled metric per context. The final score is the average of the 3 scaled metrics across the contexts. 

$$
\underbrace{
\begin{pmatrix}
r_{1,1} & r_{1,2} & r_{1,3} & r_{1,4} & r_{1,5} & r_{1,6} \\
r_{2,1} & r_{2,2} & r_{2,3} & r_{2,4} & r_{2,5} & r_{2,6} \\
r_{3,1} & r_{3,2} & r_{3,3} & r_{3,4} & r_{3,5} & r_{3,6}
\end{pmatrix}
}_{3 \times 6 \text{ Raw}}
\xrightarrow{\text{scale within context}}
\underbrace{
\begin{pmatrix}
s_{1,1} & s_{1,2} & s_{1,3} & s_{1,4} & s_{1,5} & s_{1,6} \\
s_{2,1} & s_{2,2} & s_{2,3} & s_{2,4} & s_{2,5} & s_{2,6} \\
s_{3,1} & s_{3,2} & s_{3,3} & s_{3,4} & s_{3,5} & s_{3,6}
\end{pmatrix}
}_{3 \times 6 \text{ Scaled}}
\xrightarrow{\text{avg across metrics}}
\underbrace{
\begin{pmatrix}
\bar{s}_1 \\
\bar{s}_2 \\
\bar{s}_3
\end{pmatrix}
}_{3 \times 1}
\xrightarrow{\text{avg across contexts}}
\boxed{S_{\text{final}}}
$$

Where $r_{c,m}$ is the raw value of metric $m$ in context $c$, $s_{c,m}$ is the scaled value, and $\bar{s}_c = \frac{1}{6}\sum_{m=1}^{6} s_{c,m}$.

Definitons: 
- Control profile = gene expression profile of the cell without anything knocked out
- Predicted profile = our predicted expression profile of the cell with one gene knocked out
- Actual profile = experimental gene expression profile of the cell with one gene knocked out
- Effect vector (predicted) = the difference in gene expression of all genes in our predicted profile and the control profile
- Effect vector (actual) = the difference in gene expression of all genes in the actual experimental profile and the control profile
- Recall there are 18533 genes total and 300 of them are to be knocked out individually.

Scaling formula applicable to all metrics: $\frac{u - b}{r - b}$ 
- Where $u$ is our prediction raw metric; $r$ is the "replicate" experimental data raw metric, and "b" is the "baseline" data's metric. 
- The baseline data calculates the average gene expression across all cells in a construct and submits that for each of the 400 cells across all 300 perturbations.

**Perturbation discrimination (PDS):** there are 900 total experiments (300 target genes across 3 contexts). Within each context the 300 actual profiles are graphed. For each of our predicted profiles (300), 1 is added to this graph at a time. Note that profiles remove the 300 genes (of the 18533) which are the ones being simulated. The vectors in this graph have (18233 genes). Also recall that for one experiment we generated 400 different cells; we average these 400 cells together to get a single vector of 18233 genes. 
- With the 301 profiles in the graph (300 actual and 1 of our predicted) the cosine distance between our predicted profile and all the actual profiles are generated and sorted by lowest distance. 
- Amongst the 300 actual profile 1 corresponds to the experiment in which our predicted profile attempted to simulate. Ideally, we want our predicted profile to be very close to this profile as that indicates we simulated it well. 
- $PDS = 1 - \frac{Rank}{n - 1}$ where $n = 300$ and rank means how many other actual profiles are closer to our predicted than the actual profile corresponding to the perturbation experiment our predicted profile attempted to replicate. 
    - $PDS \in [0,1]$ if we guessed the same thing for all the experiments then the value would be 0.5 by expectation. **After scaling the range is technically unbounded but generally within:** $~[-1.17, 1.17]$
    - Ex. if our predicted profile was closest to the actual profile of the experiment it attempted to predict then rank would be 0 and PDS would be 1
    - Ex. if our predicted profile was farthest from the actual profile then rank would be 299 and the PDS would be 0. 
The PDS score for each of the 300 target genes is then averaged within each context to generate the raw PDS score for that context. 

**Expression accuracy (MSE)**: Again we operate intracontext and instead of dropping all 300 target genes we only drop the 1 target we are evaluating (vectors have 18532). After generating the expression errors for the 300 targets we combine those values and scale. $MSE \in [0,1]$ following **scaling**. 

Effectively, per perturbation we calculate the ratio of distance of predicted to actual and actual to control for all perturbations (300). This ratio looks like: 

$N_p = ||pred - actual||^2 - noise$ and $D_p = ||actual - control||^2 - noise$

To combine into a single metric per perturbation we find: $\displaystyle \frac{\sum_{p=1}^{300} N_p}{\sum_{p=1}^{300} D_p}$


**Differential Expression direction fidelity (FID)**: Intra-context. For each perturbation there is a differential gene expression (pred - control) in which some of the differences are statisically significant as found by Wilcoxon-rank sum test. For the actual experiment, differential gene expression (actual - control) and the statiscially significant genes can also be found. This metric compares which genes our predicted effect vector found significant and whether these genes were up/downregulated and compares it to the actual effect vector. 

$F_p = \frac{k}{max(n_{actual}, n_{pred})}$

Where $k$ is the number of significant genes which moved in the same direction as actual and $n$ refers to the number of identified genes. 

$F_p$ is then averaged across perturbations (bounded raw in $F_p \in [0,1]$) and scaled intra-context and is not bounded but generally within $FID \in [-1.85, 1]$

**Differential Expression significance overlap (JAC)**: Intra-context. For each perturbation, similar to FID we find the significant genes in our predicted effect vector and the actual effect vector. Then we calculate the proportion of genes which our prediction identified as significant relative to the sum of the number of significant genes across our predicted and actual. 

$J_p = \frac{|Pred \cap Actual|}{|Pred \cup Actual|}$ 

Average J_p for all perturbutations where the raw is bounded as $J_p \in [0,1]$ but scaled range is not clamped and is generally between $[~0, 2.85]$

Reference doesn't do that well here and leaves the most room for positive scores. 

**Differential Expression Direction reach (reach)**: Intra-context. For each perturbation,significant genes from actual. We sort these genes by significance (p-value). 

The algorithm is to iterate through the list and at each element calculate the proportion of the genes seen that have the correct direction in the perturbed state (effectively expanding window). The algorithm terminates when all genes have been seen. It returns the deepest point (most genes seen) where the proportion of correct directions is >= 0.9.

$k^* = \max\{k : P(k) \geq 0.9\}$, taking $k^* = 0$ when no prefix clears it.

$R_p = \frac{k^*}{n_{real}}$

Note: purity $P(k)$ is not monotone — it can dip below 0.9 and recover. A wrong call at position 1 doesn't necessarily mean $k^* = 0$; if later calls dilute the error back above 0.9 (e.g. 9 correct out of 10 = 0.9), the deeper prefix counts.

We then average $R_p$ across perturbations with non-empty budgets. Raw is bounded $R_p \in [0,1]$ and scaled intra-context.

**Differential Expression log fold change accuracy (NMAE)**:  Intracontext. For each perturbation, we take the genes the actual identified as significant and passed a quality gate. We calculate the log fold change (lfc) on the predicted effect vector and the actual effect vector. The mean absolute error is the absolute difference between the lfcs for each gene between these effect vectors. 

$lfc_{pred} = log_2(pred - control)$ and $lfc_{actual} = log_2(actual - control)$  
- Need to have log fold changes such that a perturbation which doubles or halves would have the same effect on the error

$\displaystyle NMAE_p = \frac{\sum_{g \in Genes}|lfc_{real} - lfc_{actual }|}{\sum_{g \in Genes} |lfc_{actual}|}$ 

We then average the $NMAE_p$ across all perturbations which has raw value of range $NMAE \in [0, +inf)$ which is scaled to $[-6,~1]$
