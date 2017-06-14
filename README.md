
# A hierchical bayesian implementation

The Hierarchical Bayesian Model is set up this way:

α ~Exponential(λ)
β ~Dirichlet(1)
θp ~Dirichlet(αβ)
S|θp ~Multinomial(θp)
Ri|S ~Multinomial(Si)

There are four levels of abstraction in the model: α  and β are population-level hyperparameters: β is a vector representing the distribution of strategies across a population (monkeys, kids, tsimane); and α controls how uniform the strategies are for each participant relative to the population. I set λ to 0.1 representing a slight preference for non-uniformity in strategies. θpθp is a vector representing a distribution over strategies for each participant. Each participant is then a multinomial over strategies given their θ. Each participant’s responses RiRi are draws from each of the strategeies Si∈S.