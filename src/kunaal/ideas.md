First step is for models to reliably learn distinguish contexts. 
- Will likely need more data besides just 3 cell lines

Then intra-context we need to train models to recognize cells in different life-cycles/functions (e.g. actively dividing, quiescent, apoptosing)

We need to ensure that these intra-context states are biologically relevant (i.e. add in knowledge base)
- Give models understanding of what the genes do and how a combination of certain genes being on means a certain function is likely occurring
- Give models understanding of what certain cells do in general (i.e. neurons are effectively quiescent and fibroblasts secrete extracellular matrix.)

Once we have reliable intra-context states we need to create a hierarchical graph of the relationships between genes. 
- We will teach models the weights of the edges
- With that a model should be able to calculate downstream effects of perturbations. 