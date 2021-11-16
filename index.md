---
title: "Functional groups"
output:
  html_document:
    toc: true
    toc_float: true
    keep_md: yes
runtime: shiny
---





## Sequence pre-processing

Summary statistics for each of the functional groups declared in the app.

The app includes the P1 and P11 naive datasets and the P4 non-naive dataset.

For P1 and P11 the following filtration criteria were applied:

* Functional sequence, no stop codons or frame shifts.
* Sequences which start from position 1 of the V gene.
* Sequences which didn't have gaps open (-) and didn't include any N's
* After changing into group annotations, sequences which had more than a single assignment in naive repertoires were remove.

The groups were created with similarity of 95% based on complete linkage and functional sequences and up to position 318.

## Projects sequence depth

![](index_files/figure-html/unnamed-chunk-2-1.png)<!-- -->
