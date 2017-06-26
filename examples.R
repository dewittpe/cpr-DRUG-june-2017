# file:   examples.R
# author: Peter DeWitt
#
# Examples of manipulating and building calls.

library(cpr)
library(magrittr)
library(lme4)

# match.fun --------------------------------------------------------------------
#
args(cpr:::cp.formula)

# For building the initial control polygon, the end use can specify the
# regression method.  Default is lm from the stats package.

cp_default <- cp(log10(pdg) ~ bsplines(day),
                 data = spdg,
                 keep_fit = TRUE)

# Just as easily we can build a mixed effects model via lme4::lmer
cp_lmer <- cp(log10(pdg) ~ bsplines(day) + (1 | id),
              data = spdg,
              method = lmer,
              keep_fit = TRUE)

# Check the class and summary of the fits to show they are as expected.
class(cp_default$fit)
class(cp_lmer$fit)

summary(cp_default$fit)
summary(cp_lmer$fit)

# How does the argument `method` allow for the two different regression types?
?match.fun

# Important: "... avoiding undesired matching of objects of other types."

is.function(lm)
lm <- sqrt(3)   # DON'T DO THIS, FOR ILLUSTRATIVE PURPOSES ONLY
lm

mylm <- match.fun(lm)
is.function(lm)
is.function(mylm)

identical(stats::lm, mylm)

fit1 <- stats::lm(mpg ~ wt, data = mtcars)
fit2 <- mylm(mpg ~ wt, data = mtcars)
identical(coefficients(fit1), coefficients(fit2))

# clear the workspace
rm(list = ls())

# match.call -------------------------------------------------------------------
# before moving onto using match.fun in your own function, you'll likely need to
# know about match.call
?match.call

foo <- function(x, ...) {
  match.call()  # return the call
}

foo()

foo(4, arg3 = 5)

# A call can be coerced to a list, and therefore, can be modified.
as.list(foo(4, arg3 = 5))

# We will use this in our next example
rm(list = ls())



# EXAMPLE ----------------------------------------------------------------------
# build a function myregression to use any regression method

#' @param formula a regression formula
#' @param data the data.frame containing the variables used in the formula
#' @param method the regression method
#' @param ... other arguements passed to the regression \code{method}

# We'll build this function in steps to illustrate the process

# STEP 1.  Let's just look at modifying call
myregression <- function(formula, data, method, ...) {
  cl <- as.list(match.call())
  cl <- cl[-c(1, which(names(cl) %in% "method"))]
  cl
}

# Omitting the data argument
myregression(mpg ~ wt + hp, method = lm)

# including the data argument
myregression(mpg ~ wt + hp, mtcars, method = lm)

# arguments in a different order,
myregression(method = lm, data = mtcars, formula = mpg ~ wt + hp)

# Additional arguments (subset)
myregression(mpg ~ wt + hp, data = mtcars, method = lm, subset = cyl == 4)

# STEP 2: using match.fun and do.call to fit the regression models
myregression <- function(formula, data, method, ...) {
  regression <- match.fun(method)

  cl <- as.list(match.call())
  cl <- cl[-c(1, which(names(cl) %in% "method"))]

  fit <- do.call(regression, cl)

  if (isS4(fit)) {
    fit@call <- match.call()
  } else {
    fit$call <- match.call()
  }
  fit
}

lm(mpg ~ wt + hp, data = mtcars)
myregression(mpg ~ wt + hp, data = mtcars, method = lm)

glm(I(mpg > 20) ~ wt, data = mtcars, family = binomial())
myregression(I(mpg > 20) ~ wt, data = mtcars, method = glm, family = binomial())

lmer(mpg ~ wt + hp + (1 | am), data = mtcars)
myregression(mpg ~ wt + hp + (1 | am), data = mtcars, method = lmer)

# Update -----------------------------------------------------------------------
# stats::update is a very useful function
?stats::update
?stats::update.formula

# Simple example
fit <- lm(mpg ~ wt, data = mtcars)
fit
update(fit, data = subset(mtcars, hp < 180))
update(fit, subset = hp < 180)

update(fit, formula = . ~ hp)
update(fit, formula = . ~ . + hp)

# Now, let's look at the use of update with a cpr_cp object

init_cp <- cp(log10(pdg) ~ bsplines(day, df = 14) + (1 | id),
              data = spdg)
summary(init_cp$fit)

# Forgot to specify the correct method
init_cp <- update(init_cp, method = lmer, keep_fit = TRUE)
summary(init_cp$fit)

# change the formula, move back to using stats::lm for speed in these examples
# and we'll have fixed effects for age too.
init_cp <- update(init_cp, formula = . ~ . - (1 | id) + age, method = lm)

init_cp$call

# Now, let's try to change the arguments passed to the bsplines call.   Let's
# use order = 3 instead of the default order = 4

# this gives the wrong df and we've lost age
new_cp <- update(init_cp, formula = . ~ bsplines(day, order = 3))
new_cp$call

# This will error because there would be two bsplines() calls
new_cp <- update(init_cp, formula = . ~ . + bsplines(day, order = 3))

# Need to remove the old bsplines call.  It must be written exactly as it was in
# in the init_cp call.

# This errors
new_cp <- update(init_cp, formula = . ~ . - bsplines(df = 14, day) + bsplines(day, order = 3))


# this works:
new_cp <- update(init_cp, formula = . ~ . - bsplines(day, df = 14) + bsplines(day, df = 14, order = 3))
new_cp$call

# Well, this sucks.  What if we didn't know for sure all the arguments or what
# about...
init_cp <- cp(log10(pdg) ~ bsplines(day, iknots = c(0.12, 0.15)), data = spdg)

# and add one knot.
new_cp <- update(init_cp, formula = . ~ . - bsplines(day, iknots = c(0.12, 0.15)) + bsplines(day, iknots = c(-0.8, 0.12, 0.15)))

# Recall: for the CPR algorithm we may need to start with 50, or even 100
# internal knots.

# How do we modify the value of arguments to a function within a formula which
# is part of another call?

# Let's start by looking at the formula.  We can always pull of just the formula
# from the call.
f <- log10(pdg) ~ age + bsplines(day, df = 14)
f
is.list(f)
as.list(f)
is.recursive(f)

f %>% as.list %>% lapply(as.list)

# So, formula are recursive lists. We can use that.  Let's build a generic
# function update_func_args to update arguments (...) of function fun within the
# formula f.
update_func_args <- function(f, fun, ...) {

  find_func <- function(x, fun = fun, dots = dots) {

    if (is.call(x) && grepl(fun, deparse(x[[1]]))) {
      for (d in names(dots)) {
        x[[d]] <- dots[[d]]
      }
      x
    } else if (is.recursive(x)) {
      as.call(lapply(as.list(x), find_func, fun, dots))
    } else {
      x
    }
  }

  dots <- as.list(match.call(expand.dots = FALSE))$...
  out  <- lapply(as.list(f), find_func, fun, dots)
  out  <- eval(as.call(out))
  environment(out) <- environment(f)
  out
}

update_func_args(f, "bsplines", df = 16, iknots = c(0.2, 0.8))

update_func_args(f, "bsplines", df = NULL, iknots = c(0.2, 0.8))

# To use with init_cp
init_cp <- cp(log10(pdg) ~ bsplines(day, iknots = c(0.12, 0.15)), data = spdg)
init_cp$call

new_cp <-
  update(init_cp,
         formula = update_func_args(formula(init_cp), "bsplines", iknots = NULL, df = 14))
new_cp$call
plot(init_cp, new_cp, color = TRUE, show_spline = TRUE)


# Here is an example for why the is.call(x) is needed.  Consider the diamonds
# data set from ggplot2.  price as function of the cut, color, and carat.  We'll
# bin the carats with five breaks.  Then update to seven breaks and included the
# lowest value.

data("diamonds", package = "ggplot2")
fit <- lm(price ~ cut + color + cut(carat, breaks = 5),
          data = diamonds)
summary(fit)

updated_fit <-
  update(fit,
         formula = update_func_args(formula(fit),
                                    fun = "cut",
                                    breaks = 7,
                                    right = FALSE))
summary(updated_fit)

# The function update_bsplines in the cpr package is provided so that the end
# user need not write out all of the above.  Also, it will ignore silly
# arguments, or rather, arguments that mean nothing to bsplines()
init_cp$call
new_cp <- update_bsplines(init_cp, df = 10, order = 3, iknots = NULL, silly = "yep")
new_cp$call
plot(init_cp, new_cp, color = TRUE, show_spline = TRUE)

# end of file ------------------------------------------------------------------
