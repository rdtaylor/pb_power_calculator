Plackett-Burman design power calculator
=======================================

Dependencies:
-------------

* FrF2
* foreach
* doMC


Usage:
------

In order to run the either execute the script directly from the command-line
after providing execution privileges:

    $ chmod +x power.r
    $ ./power.r

or source from within the R command-line interface:

    > source('power.r')

To modify the number of runs, number of simulations to run, etc.,
see the getConfigurations() function at the top of the file.


General approach:
-----------------

The script randomly simulates normally distributed data and offsets the
observations where the "effect_column" is high by "effect_size."  The program
runs ANOVA and reports the proportion of simulations wherein the null
hypothesis is correctly rejected.  Lastly, the script reports a confidence
interval, assuming a normal distribution of power:

$$\bar{y} \pm Z * SE = \bar{y} \pm Z * \sqrt{\hat{p}(1 - \hat{p}) / n}$$

where Z is dependent upon the alpha setting, $\hat{p}$ is the mean power,
and n is the number of simulations.
