# APHIS-simulation
The ALB spread simulation takes a set of trees with an initial ALB infestation, and models insect spread for a desired length of time.

# Input
The input is a dataframe of trees that are potential ALB victims. All are assumed to be Acers. The dataframe should have the following fields, with these exact case-sensitive names (there can be other fields, they'll be ignored):

- "x", "y": location coordinates. Point of origin doesn't matter, this will only be used to calculate distances between points. SPECIFY THE LINEAR UNITS BELOW (feet, meters).
- "mean_noforestdist": the mean distance to the forest edges in each of the 8 cardinal and intercardinal directions, in m. If the tree isn't in a forest landcover, this is 0.
- "dbh": tree's DBH in cm
 
OPTIONAL:
- "infested": integer, 0 being uninfested, 1 being infested. If this field is included and there is at least 1 tree with a 1 status, the trees marked as 1 will be the seed of the outbreak. If either this field is not included or all the values are 0, a select number of trees will be randomly assigned to be the seed of the outbreak (number settable below).

The "mean_noforestdist" field can be calculated with the script "Get non forest landcover in 8 directions.R."

# Output
A dataframe that is a copy of the input dataframe, with additional fields for dates of infestation and removal ("year_infested" and "year_removed").

# Risk model: determining which trees are infested
Uninfested trees have a certain probability of being infested each timestep, as a function of number of and distance to post-emergence infested source trees and various tree characteristics (see below).

Source trees are trees which are infested, post-emergence, and have not been removed. A tree which is infested is pre-emergence until the lag period has passed (see biology settings of the model). Pre-emergence trees do not infest other trees. The lag period could be set to 0 to remove the pre- and post-emergence distinction.

A tree's risk is calculated according to the assigned risk model (currently only Long Island's model ("LI") is supported).

## LI risk model
The Long Island model assesses tree risk as a function of source pressure, DBH of tree, number of nearby Acer neighbors, and landcover. 


The probability that a tree will be infested is:

*prob = &mu; * source.term * dbh.term * density.term * distance.term*

where *&mu;* is a parameter representing the upper limit of infestation probability.

*source.term = 1/(1 + (sp / &alpha;)<sup>&gamma;</sup>)*

where *sp* is a measure of source pressure and *&alpha;* and *&gamma;* are parameters. 

*sp = &Sigma;<sub>i=1</sub><sup>N</sup> &beta; * exp((dsn(distance<sub>i</sub> <sup>&delta;</sup>)))*

for N source trees within the maximum distance of 5280 feet, where *dsn*, *&beta;*, and *&delta;* are parameters and *distance<sub>i</sub>* is the distance from the source tree to the ith neighor.

*dbh.term = 1/(1 + (DBH / b<sub>1</sub>)<sup>b<sub>2</sub></sup>)*

where DBH is tree's DBH, in cm, and *b<sub>1</sub>* and *b<sub>2</sub>* are parameters.

*density.term = 1/(1 + (N<sub>acer</sub> / c<sub>1</sub>)<sup>c<sub>2</sub></sup>)*

where $N_{acer}$ is the number of Acer neighbors within 30 m of the target (will be calculated and updated by the script) and *c<sub>1</sub>* and *c<sub>2</sub>* are parameters.

*distance.term = 1/(1 + (dist<sub>lc</sub> / e<sub>1</sub>)<sup>e<sub>2</sub></sup>)*
where $dist_lc$ is the mean distance to no-forest landcover and *e<sub>1</sub>* and *e<sub>2</sub>* are parameters.

# Tree growth
In the case of a long simulation, growth may change a tree's risk over time. This script can grow trees using a linear function of DBH:

*DBH<sub>t+1</sub> = DBH<sub>t</sub> + (m*DBH<sub>t</sub> + b)*
where this year's DBH is *DBH<sub>t+1</sub>*, last year's DBH is *DBH<sub>t</sub>*, *m* is the slope of growth, and *b* is the intercept. Either *m* or *b* can be set to 0.

# Model flow
1. Calculate this year's probability of infestation for potential hosts, based on the chosen risk model.
2. Choose the trees that get infested this year. For each tree, compare a random number to its probability to make the choice. 
3. For those trees that get infested, record the year they get infested to both track time until emergence, and generally track spread in output.
4. Apply management. As parameters, 4 different probabilities are defined for detection and removal: post-emergence tree in surveyed area, pre-emergence tree in surveyed area, post-emergence tree NOT in surveyed area, pre-emergence tree NOT in surveyed area.
5. Decide what's being surveyed: identify trees that were detected in the previous year. Define a circle of radius 1.5 miles around each detected tree: all areas so marked are surveyed areas. (This is limited to 500 trees, for computer memory reasons.) If no trees were detected in the previous year, no surveys will be done.
6. Identify all post-emergence trees within surveyed areas, and use a random number against the probability of detection of these trees to choose which are detected and removed.
7. If the removal budget has not been met, identify all pre-emergence infested trees within surveyed areas. Use a random number on each, compared to the detection probability, to choose which got detected and removed.
8. Repeat steps b and c, with non-surveyed areas and post- and pre-emergence trees. All detected trees will be removed, up to the budget.
9. Apply tree growth to all live trees.
