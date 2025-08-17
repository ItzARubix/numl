# numl

Implementations of various numerical algorithms with an emphasis on accuracy.

# What numl is and is not

While numl is not currently a full numerical analysis or scientific computation library (and perhaps will never be one), it provides accurate implementations of certain algorithms for use in a wide range of projects.

I am not writing numl in an attempt to compete with [bacon-sci](https://crates.io/crates/bacon-sci). Rather, I am writing numl because I want to implement some of the algorithms that they implement in a slightly different way which I hope to be more suitable for my needs. 

As an example, bacon-sci's implementation of the numerical derivative uses five function evaluations and does not have any safeguards against overly small step values, which is something that all numerical analysis textbooks discuss extensively. The implementation in numl uses only three function evaluations while picking a step size algorithmically in a way that minimizes the possibility of problematic round-off errors.

My hope is that numl eventually becomes a fully fledged numerical analysis library. If it turns out that the resulting crate is suitable for mathematics or scientific computing, so be it. 

I do not wish for numl to become yet another in the graveyard of scientific computing libraries that haven't been updated in years on crates.io; I understand that they may be dormant because there are no improvements left to be made for them, but the lack of upkeep in all of those crates makes it hard to tell which ones are quality. 

# Contributing

Right now, numl is extremely young and I am inexperienced. This is my first crate and I am not yet good at writing documentation and writing APIs in a way that is most useful to the end user. I would be ecstatic if anyone was willing to help, not just with implementing algorithms but also with making the APIs cleaner, helping with comments and docs, and providing any advice in general.
