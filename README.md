# ZHConnect
Connectivity analyses within Canton ZH.

## (Sub)project workflow
### Preliminary work (Philipp):
- The canton is dissected into (ideally) homogeneous patches of similar and reasonable size to be considered potential future protected areas (pFPAs).
- For the habitat type of interest, a resistance map is created, based on SDMs.

### Connectivity analysis
Samples are created by adding random subsets of the patches to the existing PAs. For each sample _i_, a Δ-analysis is performed in Graphab: For each patch _j_ within the sample, a connectivity metric (we used the Harary index) is calculated for the entire network. Subsequently, the same metric is calculated after removing patch _j_. The difference Δ is attributed to patch _j_. The results are summarised in a raster file which contains values for every patch present in sample _i_. In a second step, an average Δ metric is calculated across all samples.

```mermaid
flowchart TD
    subgraph E ["For sample i; i ∈ ℕ, i < N"]
        A["Existing PAs"] --> C[Sample]
        B["Potential habitat patches"] -- 1: Random subset --> C
        C -- 2: Graphab analysis --> D[ΔH for each patch in the sample]
    end
    subgraph K [General workflow]
        F[Sample 1] --> J[Averaged ΔH raster]
        G[Sample 2] --> J
        H["..."] --> J
        I[Sample N] --> J
    end
    E --> H
```

#### Considerations:
- The number of samples N should be large enough to get stable estimates of the connectivity metric.
- However, computation time t(N) is a linear function with a positive slope. To limit computation effort, N must be set to a reasonable value. Preliminary analyses showed relatively stable estimates when each pFPA was included in about 30 samples.
- Connectivity metrics often use parameters. For Harary index, the relevant parameter is a maximum distance _d_ below which two patches are considered connected. A reasonable estimate of _d_ must be provided. For single-species analyses, the species-specific value may be obtained from the literature. For multi-species analyses, a small (but no too small) value can be used to ensure most species of interest are able to migrate across that distance.
- Since we are using resistance maps, distance _d_ must be provided in the unit of cumulative cost across a path. Euclidean distance is converted to cumulative cost by 1) creating a graph with a high threshold value, 2) establishing a relationship between euclidean distance and cumulative cost within a relevant range, and 3) using the parameters to translate the selected euclidean distance value into a cumulative cost value.

<p align="center">
  <img src="./fig/Distance_conversion.svg" alt="Relationship between euclidean distance and cumulative cost" width="300"/>
</p>

#### ad 1: Random subset
A random subset is created from the existing PAs in combination with samples from pFPAs.

The following flowchart describes the sampling algorithm:
`NSAMPLES`= N
`CUMAREA` = Target area (total additional protected area for a future scenario)

Patches are provided as polygons with an ID and an area.
```mermaid
flowchart TD
  A[Start] --> B[Initialize variables<br>Polygon IDs: pol_ids<br>Polygon areas: pol_areas<br>Cumulative area: cum_area = 0]
  B --> C{i < NSAMPLES?}
  C -- Yes --> D[Initialize sampling]
  D --> E{cum_area < CUMAREA?}
  E -- No --> Z[Save sample as raster.<br>Increase i by 1.]
  Z --> B
  E -- Yes --> F[Select random patch from the remaining patches]
  F --> G{Is pol_ids empty?}
  G -- Yes --> I[Initialize pol_ids and pol_areas]
  G -- No --> H{Is the sampled patch already in the sample=pol_ids?}
  H -- No --> J[Add the patch ID to pol_ids and its area to pol_areas]
  H -- Yes --> K[Skip]
  J --> L[Update cum_area]
  K --> L
  L --> M[Remove sampled patch from remaining patches]
  M --> N{Is remaining patches empty?}
  N -- Yes --> O[Reset remaining to full set of polygons]
  N -- No --> E
  O --> E
  C -- No --> P[End]
```
This algorithm is implemented in /rsc/sample_target_area.R

#### ad 2: Graphab analysis
The Graphab analysis consists of the following steps
1. Create a project.
   - a) Select the respective sample raster as input landscape.
   - b) Habitat has the value 1, background the value 0.
   - c) Connexity is set to 8 (Queen's case).
2. Create a complete linkset.
   - a) Provide the resistance map of the respective habitat type as cost raster.
   - b) Set the maximum link distance to _d_.
4. Create a graph.
5. Calculate the Δ metric for each patch within the sample.
