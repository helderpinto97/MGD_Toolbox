# MIR_Greedy_Decomposition (MGD) Toolbox

Mutual Information Rate (MIR) is a key method for quantifying the dynamical coupling between two processes in a network. However, MIR calculations can vary significantly depending on the conditions, as high-order dependencies among system processes play a significant role. 

In this work, we present a methodology to identify multiplets of variables that either maximize or minimize dynamic coupling. This approach decomposes the maximal MIR into unique, redundant, and synergistic components, enabling a quantification of the relative importance of high-order effects compared to dyadic interactions.

## Repository Contents

- **`functions/`**: Contains the toolbox for identifying multiplets of variables that maximize or minimize dynamic coupling. These tools facilitate the analysis of high-order effects and their significance relative to dyadic effects.
  
- **`Toy_VAR_Simulation`**: A script demonstrating the usage of the toolbox through a toy Vector AutoRegressive (VAR) simulation example.

## References

H. Pinto, Y. Antonacci, V. R. Vergara, L. Faes, and A. P. Rocha, "Assessing Redundancy and Synergy in Brain-Heart Interactions: A Conditioning Approach," 2024 13th Conference of the European Study Group on Cardiovascular Oscillations (ESGCO), Zaragoza, Spain, 2024, pp. 1-2, doi: 10.1109/ESGCO63003.2024.10766979.
