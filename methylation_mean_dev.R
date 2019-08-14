"""
#!/usr/bin/Rscript
## Date Created: 18-Jun-2019
## Author: N. S. Araujo

## Test whether the mean methylation in a group of transcripts is significantly greater than the general mean.
## Having the mean and SD from all the transcripts (population) I can test if the mean of the population subset (my transcript group) significantly differ from the the general mean using a z-test.
"""

install.packages("BSDA")
library("BSDA")

install.packages("TeachingDemos")
library(TeachingDemos)

#########################
### For B. terrestris

## DET mC mean is greater than the general mean
## Bombus general mean (mu, stdev) vs DET mean (x = DET mean, n = number of DET)
z.test(x=0.78, mu = 0.66,  stdev = 1.29, alternative = c("greater"), n=1203 )

# test output
"""
	One Sample z-test

data:  0.78
z = 3.2264, n = 1.2030e+03, Std. Dev. = 1.2900e+00, Std. Dev. of the sample
mean = 3.7193e-02, p-value = 0.0006267
alternative hypothesis: true mean is greater than 0.66
95 percent confidence interval:
0.7188236       Inf
sample estimates:
mean of 0.78 
0.78 
"""

## Transcripts highly expressed in nurses mC mean is greater than the general mean
## Bombus general mean (mu, stdev) vs nurses mean (x = nurses mean, n = number of DEG)
z.test(x=0.95, mu = 0.66,  stdev = 1.29, alternative = c("greater"), n=436 )

# test output
"""
One Sample z-test

data:  0.95
z = 4.6941, n = 436.00000, Std. Dev. = 1.29000, Std. Dev. of the sample
mean = 0.06178, p-value = 1.339e-06
alternative hypothesis: true mean is greater than 0.66
95 percent confidence interval:
0.8483813       Inf
sample estimates:
mean of 0.95 
0.95 
"""

## Transcripts highly expressed in foragers mC mean is greater than the general mean
## Bombus general mean (mu, stdev) vs foragers mean (x = foragers mean, n = number of DEG)
z.test(x=0.73, mu = 0.66,  stdev = 1.29, alternative = c("greater"), n=767 )

# test output
"""
One Sample z-test

data:  0.73
z = 1.5028, n = 767.000000, Std. Dev. = 1.290000, Std. Dev. of the sample
mean = 0.046579, p-value = 0.06644
alternative hypothesis: true mean is greater than 0.66
95 percent confidence interval:
0.653384      Inf
sample estimates:
mean of 0.73 
0.73 
"""



#########################
### For T. angustula

## DET mC mean is greater than the general mean
z.test(x=1.30, mu = 1.24,  stdev = 2.74, alternative = c("greater"), n=241 )

# test output
"""
One Sample z-test

data:  1.3
z = 0.33995, n = 241.0000, Std. Dev. = 2.7400, Std. Dev. of the sample mean
= 0.1765, p-value = 0.3669
alternative hypothesis: true mean is greater than 1.24
95 percent confidence interval:
1.009685      Inf
sample estimates:
mean of 1.3 
1.3 
"""

## Transcripts highly expressed in nurses mC mean is greater than the general mean
z.test(x=1.57, mu = 1.24,  stdev = 2.74, alternative = c("greater"), n=179 )

# test output
"""
One Sample z-test

data:  1.57
z = 1.6114, n = 179.0000, Std. Dev. = 2.7400, Std. Dev. of the sample mean
= 0.2048, p-value = 0.05355
alternative hypothesis: true mean is greater than 1.24
95 percent confidence interval:
1.233139      Inf
sample estimates:
mean of 1.57 
1.57 
"""

## Transcripts highly expressed in foragers mC mean is greater than the general mean
z.test(x=1.25, mu = 1.24,  stdev = 2.74, alternative = c("greater"), n=62 )

# test output
"""
One Sample z-test

data:  1.25
z = 0.028737, n = 62.00000, Std. Dev. = 2.74000, Std. Dev. of the sample
mean = 0.34798, p-value = 0.4885
alternative hypothesis: true mean is greater than 1.24
95 percent confidence interval:
0.6776233       Inf
sample estimates:
mean of 1.25 
1.25 
"""
