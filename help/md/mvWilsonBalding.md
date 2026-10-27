## name
mvWilsonBalding
## title
Wilson-Balding subtree prune-and-regraft move
## description
Prunes a subtree with its parent node and regrafts it elsewhere in a time tree, at a new age.
## details
`mvWilsonBalding` picks a node that is neither the root, a child of the root, nor a sampled-ancestor tip, and removes it together with its parent node. It then picks uniformly among the places the subtree can go: any branch older than the subtree, where the new parent age is drawn uniformly along the branch, or any older fossil tip, which becomes a sampled ancestor. Unlike `mvFNPR`, it changes the attachment age and can create or remove sampled ancestors, so it moves a fossil between tip and ancestor in another part of the tree in one step.

On budding trees the parent's species flag passes to the sibling left behind, and the new parent takes the flag of the branch it splits. A sampled ancestor continues its species, and on a branch a coin decides which child continues. The root is never pruned or regrafted above.
## authors
## see_also
mvFNPR
mvCollapseExpandFossilBranch
## example
    moves.append( mvWilsonBalding(tree, weight=taxa.size()) )
## references
- citation: Wilson IJ, Balding DJ (1998). Genealogical inference from microsatellite data. Genetics 150:499-510.
- citation: Gavryushkina A, Welch D, Stadler T, Drummond AJ (2014). Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration. PLoS Computational Biology 10:e1003919.
