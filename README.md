+A hierarchical bayesian model implementation.
 +
 +There are 4 levels of abstraction in this model: 
 +
 +
 +Theta_Participant ~ Dirichlet($\alpha$ $\Beta$)
 +Strategies | Theta_Participant ~ Multinomial(Theta_Participant)
 +Responses | Strategies ~ Multinomial(Strategies)
 +
 +
 +at the bottom are participant's responses, which are
 +multinomial draws from strategies. strategies (e.g. OOMM) give
 +a dirichlet distribution over responses. 
A hierarchical Bayesian implementation

# A hierchical bayesian implementation

There are four levels of abstraction in the model: α  and β are population-level hyperparameters: β is a vector representing the distribution of strategies across a population (monkeys, kids, tsimane); and α controls how uniform the strategies are for each participant relative to the population. I set λ to 0.1 representing a slight preference for non-uniformity in strategies. θpθp is a vector representing a distribution over strategies for each participant. Each participant is then a multinomial over strategies given their θ. Each participant’s responses RiRi are draws from each of the strategeies Si∈S.