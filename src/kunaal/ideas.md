## Cell type partitioning

Cell types partition the gene expression space because different types of cells have distinct sets of genes which are turned on and off.
- Can run analysis to figure out the true variance in gene expression across cell types
    - Genes that are active across many cell types are likely "housekeeping" or simply highly fundamental genes
    - Can also learn which genes are uniquely active to certain cell types to learn a signature of different cell types and use that later.

**Core question: do models require delineation of cell types to effectively predict perturbations?**
- Regardless one experiment should certainly be to try our modeling approach post-clustering cells by type versus a general approach.

Arguments for yes:
- Different sets of genes will be completely off in one cell versus another meaning that the regulatory relationships will be entirely different
- Some genes can be permanently turned off and that is dependent on cell identity
- GRNs are context-specific: a gene that's silenced in one cell type has no outgoing regulatory edges to learn from in that context

Arguments for no:
- For housekeeping genes active across all cell types we can probably predict perturbations without identifying cell identity first.
- Potentially the regulatory relationships between genes are independent of being able to predict a perturbation

Seems like identification of housekeeping genes versus cell type specific genes will be important.
- **Theoretical question: do combinations of genes being on make the gene products cell type specific or is the gene product itself enough?**
- A hybrid approach could work: shared regulatory edges for housekeeping genes (conserved across contexts), context-specific edges for the rest. This reduces the amount of per-context data needed.

To distinguish between cell types can turn a single cell into a vector with each dimension as a gene and then apply distance-based clustering approaches to figure out which cells are likely to be a similar type based on distances in the gene expression space.
- If we are using data across multiple datasets it is almost certain that the gene expression space will be heterogeneous and thus we'll have to choose an algorithm which accounts for this (Kamila, latent class models, k-prototypes, UFT-k-means) or apply dimensionality reduction beforehand to create a unified space.

**Key experiment:** Compare perturbation prediction accuracy with GRNs built per-cell-type vs a single general GRN. The gap tells us how much context-specificity matters.

## Intra cell type state partitioning

Cells within a single type are not identical and differ in gene expression depending upon:
- **Replication / cell cycle**: cell cycle activities greatly differentiate what a cell is going to have active. Probably the single biggest source of intra-type variation in scRNA-seq — often the dominant signal after cell type itself, to the point that many analyses explicitly regress it out to see anything else.
    - **Is this more important than the cell type?** For overall expression variance, no (cell type determines which genes are on at all). But for perturbation response it could be *more important* for specific genes — knocking out a cell-cycle gene will have dramatically different effects on a dividing cell vs a quiescent one.
- **Internal functional state**: what exactly the cell is doing itself, most cells do many different things
    - E.g. fibroblast versus fibrocyte (quiescent and doesn't secrete ECM) though these are the same cell type
- **Epigenetic state**: DNA methylation, histone modifications, chromatin accessibility. These create stable differences between cells of the same type *without* genetic mutations. Some cells may have a gene's promoter methylated while others don't — same cell type, same genotype, but different regulatory capacity. This matters for perturbation prediction: if gene Y's promoter is closed in some cells, knocking out an upstream regulator of Y won't change Y's expression in those cells (it was already off). Creates subpopulation-specific perturbation responses that a uniform fold-change would miss.
- **Transcriptional bursting / stochastic gene expression**: transcription doesn't happen continuously — it happens in discrete bursts. Even genetically identical cells in identical conditions at the same cell cycle stage will have different expression profiles at any given snapshot because of when their last burst happened. This is a big part of why scRNA-seq data is so sparse and noisy. For perturbation prediction this means some cell-to-cell variation is *irreducible noise*, not biological signal. The fold-change approach handles this naturally (the noise is already in the control cells we're multiplying).
- **External microenvironment state**: things like whether the cell is in a place of damage, infection, or what are the surrounding cells doing, what is the nutrient profile, etc. Note: mostly lost in dissociated scRNA-seq — once cells are pulled apart for sequencing, spatial context is gone. Real biologically but hard to model from the data we have.
- **Mutation**: if a cell has mutations in certain genes that will create intra-cell-type differences
- **Technical capture efficiency**: non-biological but real in scRNA-seq. Different cells are captured with different efficiency, leading to variation in total UMI counts that doesn't reflect biology. Should be normalized away, not modeled — but important to distinguish from biological variation.

**Is this actually useful for perturbation modeling performance?**

Arguments for yes:
- What a cell is doing changes its expression profile significantly at least for a small set of genes
- Cell cycle state creates wildly different expression profiles
- Epigenetic state can create subpopulations that respond fundamentally differently to the same perturbation

Arguments for no:
- Could simply be additional noise and complexity models will be unable to learn
- Will we be able to effectively identify these intra cell states?
- Some variation (transcriptional bursting, technical noise) is irreducible and shouldn't be modeled

**Rough relevance ranking for perturbation prediction:**
1. Cell cycle — changes which regulatory programs are active
2. Epigenetic state — changes which genes CAN respond
3. Functional state — changes baseline expression levels
4. Transcriptional bursting — irreducible noise in measurements
5. Technical noise — should be normalized away, not modeled

**Practical note:** Even if intra-type states aren't modeled explicitly, this variability matters at prediction time — we need to produce 400 diverse cells per perturbation, not a single mean vector. The natural cell-to-cell variability captured in the controls is what gives us that diversity (see the distribution section below).

## Construction of gene regulatory networks

Imagine a graph with each gene as a node. Edges are directed (gene A → gene B) and weighted, with sign indicating activation (+) or repression (-).

We are building these at the level of cell types first, but will also try a singular general GRN, and the intra cell state GRNs as experiments.

### Data sources for edges (layered approach)

Options for generation of edges of the GRN

**External knowledge scaffold:** gives us a prior for directed edges on a subset of the gene interactions.
- STRING database (protein-protein interactions with confidence scores)
- KEGG/Reactome (curated signaling and metabolic pathways)
- ENCODE/TRRUST (transcription factor → target gene relationships)

**Control co-expression:** Association of the genes.
- Approach 1: generate a correlation matrix of gene expression
  - lacks directionality and A -> B -> C relationships still generate high correlations between A and C despite B as the potentially crucial moderator variable
- Approach 2: Predict each gene's gene expression using the gene expression of all other genes as features. Then run feature importance to see which of the other gene's were most useful in that prediction. The feature importance becomes the edge weights. Regression-based methods (GRNBoost2, GENIE3, Graphical LASSO).
  - GRNBoost2 and GENIE3 are specfic implementations of boosted tree relationships for gene regulatory networks
  - superior to approach 1 as it includes directionality and likely allows for A -> B -> C relationships to emphasize the direct variables
  - SHAP values will have to be used to generate the signed importance

**Embeddings:** Gene embeddings from foundation models (Geneformer, scGPT) encode latent functional relationships learned from millions of cells across many cell types. Embedding distances could be the edges of the GRN itself.

**Perturbation data:** Allows for derivation of causal relationships. 
- Use matched scRNA-seq control and perturbation data to figure out the actual new expression of cells when an arbitrary gene is knocked out. 
- Create a data matrix of cells and gene expression before knockout and another data matrix of cells and gene expression post knockout then run a statistical test to determine effect size. If the effect is significant then an edge exists with the weight as the effect size

We will use multiple data sources (2025 VCC training data, Arc Virtual Cell Atlas, public Perturb-seq datasets like Replogle 2022, Norman 2019) creating a heterogeneous gene space. 
- Batch correction via [Harmony](https://cran.r-project.org/web/packages/harmony/vignettes/Seurat.html) (integrates multiple scRNA-seq datasets into a single dataset by iteratively clustering cells and applying linear correction vectors in a PCA embedding space to align cells), also [scVI](https://docs.scvi-tools.org/en/1.0.0/tutorials/notebooks/api_overview.html)
- Dimensionality reduction into a shared gene space
- Train variational autoencoder to map cells to the same low-dimensional latent representation
- If we did reduce gene space below the VCC set then we would use a decoder for actual translations to generate the perturbed expression profile. 

An idea to merge these approaches together: 
- Start with external knowledge base to identify which edges certainly exist and what their directions might be
- Create an embeddings model which shows the distances between genes and use that to help inform where new edges may exist
- Create edges using regression models like GRNBoost2/GENIE3 on a subset of the genes which were identified through the embeddings and external knowledge base as likely to have a relationship and run feature importance to build the directed edges. 

### Open questions

**How sparse should the graph be?** 18K genes = 324M possible edges. Most should be zero. The density threshold is a hyperparameter to tune against perturbation prediction accuracy.

**How much of a GRN transfers across cell types?** If we learn edges from one cell type's perturbation data, what fraction are valid in a different cell type? 
- This is an interesting direction to explore as housekeeping and other genes which have similar expression profiles across cell lines probably can be represented by a general GRN. 

## Construction of transition functions

Transition functions represent the rule which moves one gene regulatory network to the next following a perturbation. This is both signal propagation and graph transformation, operating at different timescales:

We can turn the construction of transition functions into a supervised learning problem by using the matched control and perturb seq data. 

To learn it we can use graph neural networks (GNN) which takes in the GRN graph, features (the control seq data), and perturbation mask (the matched result) and then predicts the new expression (GEARS is an example option).

### What happens when a gene is knocked out

When CRISPRi knocks out gene X:
1. X's expression drops to ~0
2. X's outgoing edges still "exist" as cellular machinery, but carry no signal (no protein being made)
3. Every gene that X regulated loses that regulatory input and adjusts its expression
4. Those changes propagate further downstream through the network
5. If any downstream genes are themselves regulators, their edges effectively deactivate too — **the graph itself changes**
6. The cell reaches a new steady state with a modified effective regulatory network

### Iterative propagation algorithm

1. Start with control GRN (graph G) and control expression vector x
2. Set perturbation target gene k: x[k] = 0
3. **Signal propagation step:** propagate expression changes through G → get updated expression x'
4. **Graph update step:** for any gene whose expression dropped below a functional threshold, deactivate its outgoing edges → get updated graph G'
 - Effectively prune that edge rather than continuously propagating it
5. Repeat steps 3-4 with (G', x') until convergence → final steady state (G*, x*)
6. x* is the predicted mean expression profile for the perturbation

### Parameters to learn from perturbation data

- The functional threshold: what expression level = "effectively off" as a regulator?
- Propagation dynamics: dampening factor, max iterations, convergence criteria
- How much of the computed change is applied per iteration? 
- When do we stop iterating?

### Key questions


**Does the iterative propagation converge?** With feedback loops in the GRN, it might oscillate. May need a dampening factor or maximum iteration count. Worth testing on a small toy graph (5-10 genes) before scaling up.

**Once we have the new graph how do we derive the predicted profile?** The final steady-state expression vector x* from the converged propagation IS the predicted profile. This gets turned into a population of cells in the next step.

## How to turn a single prediction into a distribution of predictions?

Apply the predicted perturbation shift from the newly perturbed graph to the existing 400 control cells. Since we have a set of 400 control cells for each of the 46 guides we can average these together to get a unified set of 400 control cells.

Concretely, apply the perturbation effect as a **multiplicative fold-change** to individual control cells:

1. From the GRN propagation model, compute the predicted fold change per gene: `fc[g] = x_perturbed[g] / x_control_mean[g]`
2. Take the unified set of 400 control cells (averaged across the 46 non-targeting guides)
3. For each cell, multiply their expression by the fold change: `cell_pred[g] = cell_control[g] * fc[g]`
4. Round to integers, enforce non-negativity and per-cell sum constraint (< 1,000,000)

This preserves natural cell-to-cell variability (cell cycle, functional state, stochastic heterogeneity) while layering the perturbation effect on top.

**Open question:** Is a uniform fold-change across all cells realistic? Or do some cells buffer the perturbation better than others based on their state? A more sophisticated version could learn state-dependent fold changes.