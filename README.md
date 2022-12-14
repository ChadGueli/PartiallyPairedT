# Partially Paired t-Tests

## Purpose
This library is intended for comparing means with partially matched observations.
Data involving partially matched observations features two sets of observations,
$x$ and $y$, from different groups.
The data is partially matched because some but not all of the $x$ observations are
paired with a member of $y$, and the same holds true for $y$.
This pairing creates difficulties for testing for a difference means because it
invalidates the assumptions of independence between observations and groups.
Normally, with totally paired data, the solution is simply to take the difference of
observations within pairs and test the differences.
Implementing this strategy with partially paired data is not possible because
there are no observations to compare the unpaired data against.

This package provides five methods for analyzing group differences with partially
paired data; Liptak's Weighted Z-test, Kim et al.'s modified t-based statistic,
Looney and Jones’s corrected Z-test, and both MLE-based test's Lin and Stivers’s
MLE-based test under heteroscedasticity, and Ekbohm’s MLE-based test under homoscedasticity.
For example, given data on tree growth, we can look at the growth between 3 and
5 years of age. 
```{r}
x <- Loblolly[Loblolly$age == 3, "height"]
y <- Loblolly[Loblolly$age == 5, "height"]
```
Since this data is complete, we force NAs.
```{r}
x[1:5] <- NA
y[11:14] <- NA
```
Notice that there are still at least 2 values in $x$ and $y$, paired by index, that are both complete, and that there are at least two NAs in $x$ with complete compliments
in $y$, and vice-versa.
Now, we can run any of our tests on the data.
```{r}
liptak.Z(x, y)
```
We may also pass the data as a matrix.
```{r}
X <- cbind(x, y)
kim.mod.t(X)
```
Notably, Lin and Stivers’s MLE-based test is the default, but Ekbohm's will be
performed when `var.equal=TRUE`.
```{r}
mle.test(X)
mle.test(X, var.equal=TRUE)
```
Finally, the form of Looney and Jones’s corrected Z-test means that it is not
necessary for both groups to be incomplete. For example, if we replace the NAs
in $y$ with their true values, we can still run the test.
```{r}
y <- Loblolly[Loblolly$age == 5, "height"]
LJ.Zcorr(x, y)
