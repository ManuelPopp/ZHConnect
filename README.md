# ZHConnect
Connectivity analyses within Canton ZH

```mermaid
flowchart TD
    subgraph E ["Workflow sample i; i ∈ ℕ, i < N"]
        A["Existing PAs"] --> C[Sample]
        B["Potential habitat patches"] -- 1: Random subset --> C
        C -- 2: Graphab analysis --> D[ΔH for each patch in the sample]
    end
    subgraph K [General workflow]
        F[Sample 1] --> J[Averaged raster]
        G[Sample 2] --> J
        H["..."] --> J
        I[Sample N] --> J
    end
    E --> H
```
