# Change log

5th October 2026

## Serogroup

No serogroup marker changes.

## Virulence genes

| Classification | Change      | Details                                                                                                       |
|----------------|-------------|---------------------------------------------------------------------------------------------------------------|
| Added          | `als`       | Added as an individual virulence gene.                                                                        |
| Added          | `aphA`      | Added as an individual virulence gene.                                                                        |
| Added          | `aphB`      | Added as an individual virulence gene.                                                                        |
| Added          | `bap1`      | Added as an individual virulence gene.                                                                        |
| Added          | `frhA`      | Added as an individual virulence gene.                                                                        |
| Added          | `gbpA`      | Added as an individual virulence gene.                                                                        |
| Added          | `hapR`      | Added as an individual virulence gene.                                                                        |
| Added          | `hfq`       | Added as an individual virulence gene.                                                                        |
| Added          | `hlyU`      | Added as an individual virulence gene.                                                                        |
| Added          | `luxS`      | Moved out of the Lux cluster and added as an individual virulence gene.                                       |
| Added          | `mam7`      | Added as an individual virulence gene.                                                                        |
| Added          | `mfrhA`     | Added as an individual virulence gene.                                                                        |
| Added          | `ompT`      | Added as an individual virulence gene.                                                                        |
| Added          | `prtV`      | Added as an individual virulence gene.                                                                        |
| Added          | `rpoN`      | Added as an individual virulence gene.                                                                        |
| Added          | `tdh (dth)` | Added as an individual virulence gene.                                                                        |
| Added          | `tlh`       | Added as an individual virulence gene.                                                                        |
| Added          | `vqmA`      | Added as an individual virulence gene.                                                                        |
| Added          | `wbfZ`      | Added as an individual virulence gene.                                                                        |
| Removed        | `ace`       | Removed from the individual virulence-gene library; retained in the CTX prophage cluster.                     |
| Removed        | `acfA`      | Removed from the individual virulence-gene library; retained in the VPI-1 genomic-island cluster.             |
| Removed        | `acfB`      | Removed from the individual virulence-gene library; retained in the VPI-1 genomic-island cluster.             |
| Removed        | `acfC`      | Removed from the individual virulence-gene library; retained in the VPI-1 genomic-island cluster.             |
| Removed        | `acfD`      | Removed from the individual virulence-gene library; retained in the VPI-1 genomic-island cluster.             |
| Removed        | `tagA`      | Removed from the individual virulence-gene library; retained in the VPI-1 genomic-island cluster.             |
| Removed        | `vasX`      | Removed from the individual virulence-gene library; retained in type VI secretion system auxiliary cluster 2. |
| Removed        | `zot`       | Removed from the individual virulence-gene library; retained in the CTX prophage cluster.                     |

## Virulence clusters

| Classification | Change                                        | Details                                                                                                                     |
|----------------|-----------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| Added          | CTX prophage region                           | Added the 13-gene CTX prophage cluster, including `ctxA` and `ctxB`.                                                        |
| Added          | Cqs quorum sensing cluster                    | Added `cqsS` and `cqsA`.                                                                                                    |
| Added          | Type III secretion system T3SS2-alpha cluster | Added the 15-gene alpha cluster.                                                                                            |
| Added          | Type III secretion system T3SS2-beta cluster  | Added the 38-gene beta cluster.                                                                                             |
| Added          | Type VI secretion system auxiliary cluster 1  | Added `hcp1`, `vgrG1`, `tseL`, and `tsiV1`.                                                                                 |
| Added          | Type VI secretion system auxiliary cluster 2  | Added `hcp2`, `vgrG2`, `vasX`, and `tsiV2`.                                                                                 |
| Added          | Type VI secretion system auxiliary cluster 3  | Added `tseH` and `tsiH`.                                                                                                    |
| Added          | Type VI secretion system main cluster         | Added `vasH`, `vgrG3`, and `tsiV3`.                                                                                         |
| Added          | Flagellar cluster 1                           | Added `motA` and `motB`.                                                                                                    |
| Added          | Flagellar cluster 2                           | Added `fliA` and `flhA`.                                                                                                    |
| Added          | Flagellar cluster 3                           | Added `flhB`, `fliE`, `flrB`, `flrC`, `flrA`, `flaB`, and `flaD`.                                                           |
| Added          | Flagellar cluster 4                           | Added `flaC`, `flaA`, `flgK`, `flgJ`, `flgD`, `flgM`, `flgN`, and `flgP`.                                                   |
| Added          | Vibrio polysaccharide (VPS) cluster           | Added the 21-gene VPS cluster.                                                                                              |
| Added          | O-antigen locus                               | Added the nine-gene O-antigen locus.                                                                                        |
| Added          | Type II secretion system cluster              | Added the 12-gene type II secretion-system cluster.                                                                         |
| Removed        | TCP cluster                                   | Superseded by the larger VPI-1 genomic-island cluster.                                                                      |
| Added          | VPI-1 genomic island                          | Expanded the former TCP definition to 24 genes, including `acfA`–`acfD`, `aldA`, `tagA`, `tagD`, `tagE`, `mop`, and `tcpP`. |
| Renaming       | `tcpN` to `toxT`                              | Switched to the more commmonly used `toxT` name.                                                                            |
| Removed        | Lux Operon                                    | Replaced the five-gene Lux operon cluster.                                                                                  |
| Added          | Lux operon cluster 1                          | Added `luxO` and `luxU`.                                                                                                    |
| Added          | Lux operon cluster 2                          | Added `luxQ` and `luxP`.                                                                                                    |
| Removed        | `luxS` from Lux cluster membership            | `luxS` is now searched as an individual virulence gene instead.                                                             |
| Fix            | MSHA pilus membership                         | Added `MshL` and aligned member order with tab B.                                                                           |
| Fix            | RTX toxin operon member order                 | Reordered members to match tab B: `rtxD`, `rtxB`, `rtxC`, `rtxA`.                                                           |

----

*Older history is currently missing.*