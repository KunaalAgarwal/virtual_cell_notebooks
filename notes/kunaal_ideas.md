# Axes

Representation of cells

Gene regulatory network (GRN; knowledge graph): a representation with genes as nodes and edges as relationships

Transition function (f(basal cell, perturb) = perturbed cell): 
- Linear propagation approaches
- Perturbation embeddings

## Representation of cells

Individual cell types have different gene spaces (i.e. distinct sets of genes which are off/on) and potentially have different gene regulatory networks. This means models will need some way of differentiating cell types and cell state. Cell state is differences in expressions of cells of the same type due to cell cycle activity, internal functional state, epigenetics, transcriptional bursting, etc.

Approaches on representation: 
- Embedding cells: 
    - Pretrained foundation models (ex. Geneformer, scGPT) can output an embedding vector in that learned latent space.
    - Multi-layer perceptrons were used in TxPert
- Creating regulatory networks within cell contexts

It is possible that cell type and state are not extremely relevant to predictions and representing cells in these ways may not be necessary.

## Gene regulatory networks (GRN)

Graph representing baseline knowledge of the interactions of genes prior to any modeling. Genes are nodes with weighted signed edges (positive is activation and negative is repression). 

The construction of the GRN depends on the reprsentation of cells. 

### Methods for construction of GRN: 

Use external knowledge bases to create gene-gene interaction information: 
- STRING database (protein-protein interactions with confidence scores)
    - TxPert uses this
- GO (gene ontology) Database
    - TxPert and GEARS paper uses this to create their knowledge base; it creates edges between genes if there is an annotation shared in the GO knowledge base
- KEGG/Reactome (curated signaling and metabolic pathways)
- ENCODE/TRRUST (transcription factor → target gene relationships)
- Brief aside is that TxPert additionally used PxMap/TxMap which are empirically derived graphs in addition to STRING and GO. 
    - PxMap was generated through associating morphological features of cells after genetic/chemical perturbations. I.e. if cells had similar morphology after two different perturbations those perturbations were associated with one another. 
    - TxMap was generated through associating transcriptomes of cells after genetic/chemical perturbations
    - **We can try to recreate these maps**
- Another aside is that if we combine approaches TxPert showed Exphormer-MG is best for singular perturbations. 
    - It takes the union of all the graphs combining it into a single graph and runs a graph transformer over it

Create a gene association map from the vcc control data:
- Approach 1: generate a correlation matrix of gene expression
  - lacks directionality and A -> B -> C relationships still generate high correlations between A and C despite B as the potentially crucial moderator variable
- Approach 2: Predict each gene's gene expression using the gene expression of all other genes as features. Then run feature importance to see which of the other gene's were most useful in that prediction. The feature importance becomes the edge weights. Regression-based methods (GRNBoost2, GENIE3, Graphical LASSO).
  - GRNBoost2 and GENIE3 are specfic implementations of boosted tree relationships for gene regulatory networks
  - superior to approach 1 as it includes directionality and likely allows for A -> B -> C relationships to emphasize the direct variables
  - SHAP values will have to be used to generate the signed importance

### Combine knowledge graphs

Feature engineering: derive features from these graphs and use those in any kind of model (e.g. linear regression)

Gene-gene similarity matrix: develop a 18533x18533 matrix where each entry is how related each gene is to the other. 

Graph propagation

Graph neural networks: TxPert used this


## Transition functions

The construction here depends on the cell and GRN representation. Something like TxPert uses the GRN to generate a perturbation embedding which shifts the embedded cell 

For the non-GNN approaches conversion into features allows for supervised learning with other datasets and then transfer of the model to the vcc cell lines for prediction.