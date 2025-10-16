\name{mccimm}
\alias{mccimm_modsem}
\alias{JN_plot}
\alias{mccimm_modsem_fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Simulated Confidence Intervals for Moderated Mediation with modsem Object
}
\description{
Percentile and bias-corrected confidence intervals for moderating mediation created by 5,000,000 Monte Carlo Simulated samples. Can also create confidence intervals for a function.
}
\usage{
mccimm_modsem(object = est_lms, Z="NA", W="NA",
                   A1="NA", Z1="NA", W1="NA", ZW1="NA",
                   A2="NA", Z2="NA", W2="NA", ZW2="NA",
                   A3="NA", Z3="NA", W3="NA", ZW3="NA",
                   A4="NA", Z4="NA", W4="NA", ZW4="NA",
                   R=5)
\usage{
JN_plot(mccimmObject, ci="bc",
                     min_z = -3, max_z = 3, detail = 300,
                     lower.quantile = 0.025, upper.quantile = 0.975,
                     alpha = 0.2, sd.line = 2)
\usage{
mccimm_modsem_fun(object = est_lms, Z="NA", W="NA",
                   A1="NA", Z1="NA", W1="NA", ZW1="NA",
                   A2="NA", Z2="NA", W2="NA", ZW2="NA",
                   A3="NA", Z3="NA", W3="NA", ZW3="NA",
                   A4="NA", Z4="NA", W4="NA", ZW4="NA",
                   Sfunction="NULL", R=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
modsem object
}
  \item{Z}{
moderator Z
}
  \item{W}{
moderator W
}
  \item{A1}{
A1 path
}
  \item{Z1}{
Z1 path
}
  \item{W1}{
W1 path
}
  \item{ZW1}{
ZW1 path
}
  \item{A2}{
A2 path
}
  \item{Z2}{
Z2 path
}
  \item{W2}{
W2 path
}
  \item{ZW2}{
ZW2 path
}
  \item{A3}{
A3 path
}
  \item{Z3}{
Z3 path
}
  \item{W3}{
W3 path
}
  \item{ZW3}{
ZW3 path
}
  \item{A4}{
A4 path
}
  \item{Z4}{
Z4 path
}
  \item{W4}{
W4 path
}
  \item{ZW4}{
ZW4 path
}
  \item{R}{
Number of Monte Carlo Simulated Samples R x 1 million
}
  \item{mccimmObject}{
Object created from mccimm
}
  \item{ci}{
ci = "bc"(default): Bias-corrected confidence intervals for the figure;
     "percent": percentile confidence intervals
}
  \item{Sfunction}{
The defined function for simulation
}

}
\details{
	By specifying the path labels for the direct paths (a1, a2, a3, and a4); the moderating effects from moderator Z to the direct paths (z1, z2, z3, and z4); the moderating effects from moderator W to the direct paths (w1, w2, w3, and w4), the interaction effects between W and Z on the direct paths (three-way interactions zw1, zw2, zw3, and zw4), as well as the variable name of the moderator Z and/or moderator W, confidence intervals and p-values of the mediating effects or moderated mediating effects can be estimated by Monte Carlo simulation.
	The above general model has covered many models in Hayes' templates.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

One Moderator Z
Model 1: A1, Z1, Z
Model 4: A1, A2
Model 5: A1, A2 (mediating effect) or A1, Z1, Z (moderated direct effect)
Model 6: A1, A2, A3
Model 7: A1, Z1, A2, Z
Model 8: A1, Z1, A2, Z (moderated mediating effect) or A1, Z1, Z (moderated direct effect)
Model 14: A1, A2, Z2, Z
Model 15: A1, A2, Z2, Z (moderated mediating effect) or A1, Z1, Z (moderated direct effect)
Model 58: A1, Z1, A2, Z2, Z
Model 59: A1, Z1, A2, Z2, Z (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)

Two Moderators Z & W
Model 2: A1, Z1, W1, Z, W
Model 3: A1, Z1, W1, ZW1, Z, W
Model 9: A1, Z1, W1, A2, Z, W
Model 10: A1, Z1, W1, A2, Z, W (moderated mediating effect) or A1, Z1, W1, Z, W (moderated direct effect)
Model 11: A1, Z1, W1, ZW1, A2, Z, W
Model 12: A1, Z1, W1, ZW1, A2, Z, W (moderated mediating effect) or A1, Z1, W1, ZW1, Z, W (moderated direct effect)
Model 13: A1, Z1, W1, ZW1, A2, Z, W (moderated mediating effect) or A1, Z1, Z (moderated direct effect)
Model 16: A1, A2, Z2, W2, Z, W
Model 17: A1, A2, Z2, W2, Z, W (moderated mediating effect) or A1, Z1, W1, Z, W (moderated direct effect)
Model 18: A1, A2, Z2, W2, ZW2, Z, W
Model 19: A1, A2, Z2, W2, ZW2, Z, W (moderated mediating effect) or A1, Z1, W1, ZW1, Z, W (moderated direct effect)
Model 20: A1, A2, Z2, W2, ZW2, Z, W (moderated mediating effect) or A1, Z1, Z (moderated direct effect)
Model 21: A1, Z1, A2, W2, Z, W
Model 22: A1, Z1, A2, W2, Z, W (moderated mediating effect) or A1, Z1, Z (moderated direct effect)
Model 28: A1, Z1, A2, W2, Z, W (moderated mediating effect) or A1, Z1, Z (moderated direct effect)
Model 29: A1, Z1, A2; W2, Z, W (moderated mediating effect) or A1, Z1, W1, Z, W (moderated direct effect)
Model 60: A1, Z1, W1, A2, W2, Z, W
Model 61: A1, Z1, W1, A2, W2, Z, W (moderated mediating effect) or A1, W1, W (moderated direct effect)
Model 62: A1, Z1, W1, A2, W2, Z, W (moderated mediating effect) or A1, Z1, Z (moderated direct effect)
Model 63: A1, Z1, W1, A2, W2, Z, W (moderated mediating effect) or A1, W1, W (moderated direct effect)
Model 64: A1, Z1, A2, Z2, W2, Z, W
Model 65: A1, Z1, A2, Z2, W2, Z, W (moderated mediating effect) or A1, Z1, Z (moderated direct effect)
Model 66: A1, Z1, A2, Z2, W2, Z, W (moderated mediating effect) or A1, W1, W (moderated direct effect)
Model 67: A1, Z1, A2, Z2, W2, Z, W (moderated mediating effect) or A1, Z1, W1, Z, W (moderated direct effect)
Model 68: A1, Z1, W1, ZW1, A2, Z2, Z, W
Model 69: A1, Z1, W1, ZW1, A2, Z2, Z, W (moderated mediating effect) or A1, Z1, W1, ZW1, Z, W (moderated direct effect)
Model 70: A1, Z1, A2, Z2, W2, ZW2, Z, W
Model 71: A1, Z1, A2, Z2, W2, ZW2, Z, W (moderated mediating effect) or A1, Z1, W1, ZW1, Z, W (moderated direct effect)
Model 72: A1, Z1, W1, ZW1, A2, Z2, W2, ZW2, Z, W
Model 73: A1, Z1, W1, ZW1, A2, Z2, W2, ZW2, Z, W (moderated mediating effect) or A1, Z1, W1, ZW1, Z, W (moderated direct effect)
Model 75: A1, Z1, W1, A2, Z2, W2, Z, W
Model 76: A1, Z1, W1, A2, Z2, W2, Z, W (moderated mediating effect) or A1, Z1, W1, Z, W (moderated direct effect)

For single-level moderated mediation models, A1, A2, A3, and A4 are the direct effects, whereas Z1, Z2, Z3, Z4 and W1, W2, W3, W4 are two-way interaction effects between the moderator and the predictor on a criterion.  Three-way interaction effects between the moderators on the effects from the predictor on a criterion are represented by ZW1, ZW2, ZW3, and ZW4.
When there are three-way interactions, the difference in simple slope tests are provided based on Dawson (2014).
For cross-level moderated mediation models, A1, A2, A3, and A4 are the mean of random slopes where Z1, Z2, Z3, Z4 and W1, W2, W3, W4 are the effects of the Level 2 moderators on the random slopes.
The general model actually goes beyond the above models because more than two mediating effects (maximum four) can be specified.  Please note that when there are two or more moderating effects, Index of Moderated Mediation is not defined (same for PROCESS).

EXAMPLE A: (Model 7 with 7 parameters in New/Additional Parameters)
mcObject <- mccimm_modsem(est_lms, A1="a1", A2="a2", Z1="z1", Z="Autonomy")

EXAMPLE B: Generate Confidence Intervals with a function
mccimm_modsem_fun(est_lms, A1="a1", Z1="z1", A2="a2", Sfunction="(a1+z1)*a2")

EXAMPLE C: Plot Johnson_Neyman Figure
JNN <- JN_plot(mcObject, ci="bc") +
  ggplot2::ggtitle("Johnson-Neyman Figure of Moderated Mediating Effects") +
  ggplot2::labs(x = "Moderator", y = paste("Mediating Effect"))

ggplot2::ggsave("Johnson-Neyman Figure.png", plot=JNN, width = 22, height = 12, units = "cm")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
