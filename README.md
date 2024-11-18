# `neodowns`: Simulating Candidate Competition in Group-Based Elections

`neodowns` Performs computer simulations of candidate competition under different electoral systems based on the spatial model proposed in Atsusaka and Landsman (2024). Working Paper. \## Installation

`neodowns` can be installed using the following code:

``` r
remotes::install_github(
  "YukiAtsusaka/neodowns",
  dependencies = TRUE
)
```

# Example

Our workhorse function `sim_data` simulates spatial positions of voters and candidates under a set of conditions defined by model parameters. `neodowns` implements the algorithm based on the neodownsian model of candidate competition.

## Simulate candidate and voter positions

``` r
# Generate simulated data
sim <- sim_data(N_voters = 1000,
                N_groups = 3,
                dist = "Normal")

```

## Performs simulations based on the neodownsian model

``` r

# Implement the algorithm
out <- neodowns(sim, n_iter = 1000)                

head(out)
# # A tibble: 6 × 43
#    iter  delta epsilon delta_md delta_ex delta_rcv epsilon_rcv_1
#   <dbl>  <dbl>   <dbl>    <dbl>    <dbl>     <dbl>         <dbl>
# 1     1 NA       NA       NA      NA        NA            NA    
# 2     2  1.01     1.00     1.03    0.976     1.00          1.00 
# 3     3  0.996    1.00     1.04    0.956     0.981         1.00 
# 4     4  0.987    1.00     1.04    0.936     0.995         0.999
# 5     5  0.962    1.00     1.01    0.913     0.993         1.00 
# 6     6  0.973    1.00     1.05    0.890     0.980         1.00 
```
