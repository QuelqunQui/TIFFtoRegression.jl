# TIFFtoRegression

[![Build Status](https://travis-ci.org/QuelqunQui/TIFFtoRegression.jl.svg?branch=master)](https://travis-ci.org/QuelqunQui/TIFFtoRegression.jl)

[![Coverage Status](https://coveralls.io/repos/QuelqunQui/TIFFtoRegression.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/QuelqunQui/TIFFtoRegression.jl?branch=master)

[![codecov.io](http://codecov.io/github/QuelqunQui/TIFFtoRegression.jl/coverage.svg?branch=master)](http://codecov.io/github/QuelqunQui/TIFFtoRegression.jl?branch=master)

This Package is the result of personnal needs to gain access to mathematical expression from capture images from pdf files (that were often from some scan of old papers).  It might need to be split completely but can basicly load a TIFF into a Matrix, extract line from that matrix after boolean transformation and then perform a linear regression implementing ANOVA.
The Line detection is imited to continuous lines (i. e. dashed excluded) that are on a "plain" background (e. i. no-grid).  I'd like to add grid detection (if times permit) to be able to get rid of it.

## Detection functions

### FieldValues(::Int64,::Array{Int16},::Int64,[::Int64])

`FieldValues` is a utilitary fonction to read TIFF tagged values or their pointers, only small endian is implemented.

### AquiTIFF(::String)

Loads a TIFF file from a given location and reads it.  Not all tags have been implemented.  minimum tag requierements for TIFF grey scale have all been implemented according to TIFF v6.  Specific Dantec Dynamics tag have also been implemented.  Only small endian implemented.
Output is a matrix giving the level of gray of each pixel.

### AquiSerieD(::String, ::String, [::Int64])

Allows to acquiere an entire Serie of TIFF using Dantec Dynamics standard numerotation, given the folder where they are, the name of the first image and the number of images there is.

### AquiSerieW(::String,::String, [::Int64])

Allows to acquiere a Serie of TIFF using Windows standard numérotation, given the folder where they are, the name of the first image and the number there are.

## Detection of Lines in matrix
### Line
`::Line` type is created, allowing to stock values in `Line.YValues` as a `::Vector` and the starting point of the line in `Line.XStart` as a `::Int64`
### DetectLine(::Matrix)

The algorithme takes a Matrix and follow the lines it detect, given them back as `Line` types, all measured in pixels (or Matrix elements as you might call them from here on).
Algorithme works as :
 0. from bottom to top
 1. First column -> detect LinePx
   * if no LinePx go to the right
 2. continue climbing until Px!=LinePx
 3. from there moove to the right following the line
   * check same Px and Px+1 and Px-1
   * if Px+-1 check also Px+-2, etc
   * if no Px: line stop
   * repeat 3 until line stop
 4. remoove selected Pixels
 5. repeat operation until top right

### PlotLine(::Array{Any},[::Int64])

Is a tool to plot `::Line` directly based on `PyPlot`'s plot, that take as input an Array of `::Line` (i. e. the output of `DetectLine`).

## Multiple Regression Functions

### NextFun(::Vector,[::String, ::Int64])

Create a function of the numerical input given, ready to be used in a multiple regression comparaison.  Implemented type are Polynomes, Trigonometric function (as standardize cos), Fractions, Exponentials, Logarithms, a Mix-option and a User friendly choice-at-each-step also including inverse-Trigonometric functions.

### StatReg : ANOVA Implemented

The function calculate the associated Statistical functions to the regression.  The p-value corresponding to each tested members of the model is given, calculated by a Student-test and a p-value for the whole model to asses its concordancy with the data is calculated according to a Fisher-test.

### MultiRegOpt(::Vector, ::Vector,[::String,::Int64,::String, ::Float])

Implementation of the regression of Y(x).  The first vector being x and the second y.  The possible choice about the regression are the type of function used.  The choices are:
* "Poly" for polynoms
* "Custom" for stepwise choice entered by the User
* "Trigo" for standardize cos
* "Frac" for Fractions
* "Expo" for Exponentials
* "Loga" for Logarithms
* "Mixte" for a mix of the above
